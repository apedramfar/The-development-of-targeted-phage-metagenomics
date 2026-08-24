import gzip

import click

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def _open_text(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def parse_blast(handle):
    for line in handle:
        r = line.split()
        yield {
            "qname": r[0],
            "tname": r[1],
            "pid": float(r[2]),
            "len": float(r[3]),
            "qcoords": sorted([int(r[6]), int(r[7])]),
            "tcoords": sorted([int(r[8]), int(r[9])]),
            "qlen": float(r[-2]),
            "tlen": float(r[-1]),
            "evalue": float(r[-4]),
        }


def yield_alignment_blocks(handle):
    key, alns = None, None
    for aln in parse_blast(handle):
        key = (aln["qname"], aln["tname"])
        alns = [aln]
        break
    for aln in parse_blast(handle):
        if (aln["qname"], aln["tname"]) == key:
            alns.append(aln)
        else:
            yield alns
            key = (aln["qname"], aln["tname"])
            alns = [aln]
    if alns:
        yield alns


def prune_alns(alns, min_length=0, min_evalue=1e-3):
    keep = []
    cur_aln = 0
    qry_len = alns[0]["qlen"]
    for aln in alns:
        qcoords = aln["qcoords"]
        aln_len = max(qcoords) - min(qcoords) + 1
        if aln_len < min_length or aln["evalue"] > min_evalue:
            continue
        if cur_aln >= qry_len or aln_len + cur_aln >= 1.10 * qry_len:
            break
        keep.append(aln)
        cur_aln += aln_len
    return keep


def compute_ani(alns):
    return round(
        sum(a["len"] * a["pid"] for a in alns) / sum(a["len"] for a in alns), 2
    )


def compute_cov(alns):
    coords = sorted([a["qcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        else:
            nr_coords.append([start, stop])
    qcov = round(
        100.0
        * sum(stop - start + 1 for start, stop in nr_coords)
        / alns[0]["qlen"],
        2,
    )

    coords = sorted([a["tcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        else:
            nr_coords.append([start, stop])
    tcov = round(
        100.0
        * sum(stop - start + 1 for start, stop in nr_coords)
        / alns[0]["tlen"],
        2,
    )

    return qcov, tcov


@click.command(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "-i",
    "input",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    required=True,
    metavar="PATH",
    help="Path to blastn input file (format: 'std 6 qlen slen')",
)
@click.option(
    "-o",
    "output",
    type=click.Path(dir_okay=False, writable=True, path_type=str),
    required=True,
    metavar="PATH",
    help="Path to ANI output file",
)
@click.option(
    "-l",
    "length",
    type=click.IntRange(min=0),
    default=0,
    show_default=True,
    metavar="INT",
    help="Minimum alignment length to keep",
)
def cli(input, output, length):
    with _open_text(output, "wt") as out, _open_text(input) as handle:
        fields = ["qname", "tname", "num_alns", "pid", "qcov", "tcov"]
        out.write("\t".join(fields) + "\n")
        for alns in yield_alignment_blocks(handle):
            alns = prune_alns(alns, min_length=length)
            if not alns:
                continue
            qname, tname = alns[0]["qname"], alns[0]["tname"]
            row = [qname, tname, len(alns), compute_ani(alns), *compute_cov(alns)]
            out.write("\t".join(str(value) for value in row) + "\n")


if __name__ == "__main__":
    cli()

