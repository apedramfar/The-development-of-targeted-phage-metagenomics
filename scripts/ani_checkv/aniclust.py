import gzip
import platform
import resource
import time
from pathlib import Path

import click

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def _open_text(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_seqs(path):
    with _open_text(path) as handle:
        seq_id = next(handle).split()[0][1:]
        seq = ""
        for line in handle:
            if line[0] == ">":
                yield seq_id, seq
                seq_id = line.split()[0][1:]
                seq = ""
            else:
                seq += line.rstrip()
        yield seq_id, seq


def log_time(start):
    program_time = round(time.time() - start, 2)
    peak_ram = round(max_mem_usage(), 2)
    print("time: %s seconds, peak RAM: %s GB" % (program_time, peak_ram))


def max_mem_usage():
    max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    divisor = 1e6 if platform.system() == "Linux" else 1e9
    return round((max_mem_self + max_mem_child) / divisor, 2)


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "--fna",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    required=True,
    metavar="PATH",
    help="Path to nucleotide sequences",
)
@click.option(
    "--ani",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    required=True,
    metavar="PATH",
    help="Path to ANI file",
)
@click.option(
    "--out",
    type=click.Path(dir_okay=False, writable=True, path_type=str),
    required=True,
    metavar="PATH",
    help="Path to output file",
)
@click.option(
    "--exclude",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    metavar="PATH",
    help="Path to sequence IDs to exclude",
)
@click.option(
    "--keep",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    metavar="PATH",
    help="Path to sequence IDs to keep",
)
@click.option(
    "--min_ani",
    type=click.FloatRange(min=0.0, max=100.0),
    default=95,
    show_default=True,
    metavar="FLOAT",
    help="Minimum average nucleotide identity",
)
@click.option(
    "--min_qcov",
    type=click.FloatRange(min=0.0, max=100.0),
    default=0,
    show_default=True,
    metavar="FLOAT",
    help="Minimum alignment coverage of longer sequence",
)
@click.option(
    "--min_tcov",
    type=click.FloatRange(min=0.0, max=100.0),
    default=85,
    show_default=True,
    metavar="FLOAT",
    help="Minimum alignment coverage of shorter sequence",
)
@click.option(
    "--min_length",
    type=click.IntRange(min=1),
    default=1,
    show_default=True,
    metavar="INT",
    help="Minimum sequence length",
)
def cli(fna, ani, out, exclude, keep, min_ani, min_qcov, min_tcov, min_length):
    """Centroid based sequence clustering."""
    start = time.time()

    print("\nreading sequences...")
    exclude_ids = set(Path(exclude).read_text().splitlines()) if exclude else None
    keep_ids = set(Path(keep).read_text().splitlines()) if keep else None
    seqs = {}
    for seq_id, seq in parse_seqs(fna):
        if len(seq) < min_length:
            continue
        if exclude_ids and seq_id in exclude_ids:
            continue
        if keep_ids and seq_id not in keep_ids:
            continue
        seqs[seq_id] = len(seq)
    seqs = [
        seq_id
        for seq_id, _ in sorted(
            seqs.items(), key=lambda item: item[1], reverse=True
        )
    ]
    print("%s sequences retained from fna" % len(seqs))
    log_time(start)

    print("\nstoring edges...")
    num_edges = 0
    edges = dict((seq_id, []) for seq_id in seqs)
    with _open_text(ani) as handle:
        for line in handle:
            qname, tname, _num_alns, ani_value, qcov, tcov = line.split()
            if qname == tname or qname not in edges or tname not in edges:
                continue
            if (
                float(qcov) < min_qcov
                or float(tcov) < min_tcov
                or float(ani_value) < min_ani
            ):
                continue
            edges[qname].append(tname)
            num_edges += 1
    print("%s edges retained from blastani" % num_edges)
    print("%s edges currently stored" % sum(len(value) for value in edges.values()))
    log_time(start)

    print("\nclustering...")
    clust_to_seqs = {}
    seq_to_clust = {}
    for seq_id in seqs:
        if seq_id in seq_to_clust:
            continue
        clust_to_seqs[seq_id] = [seq_id]
        seq_to_clust[seq_id] = seq_id
        for mem_id in edges[seq_id]:
            if mem_id not in seq_to_clust:
                clust_to_seqs[seq_id].append(mem_id)
                seq_to_clust[mem_id] = seq_id
    print("%s total clusters" % len(clust_to_seqs))
    log_time(start)

    print("\nwriting clusters...")
    with open(out, "w") as handle:
        for seq_id, mem_ids in clust_to_seqs.items():
            handle.write(seq_id + "\t" + ",".join(mem_ids) + "\n")
    log_time(start)


if __name__ == "__main__":
    cli()

