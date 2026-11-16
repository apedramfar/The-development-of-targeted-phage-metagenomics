# ============================================================
# INPHARED2 comparison and ANI analysis of minicell-enriched vOTUs
# ============================================================

# Software:
# Mash v2.3
# BLAST+ v2.17.0
# CheckV supporting script: anicalc.py
#
# INPHARED2 database:
# 10 August 2026; 36,350 phage genomes


# ------------------------------------------------------------
# 1. Sketch individual vOTUs
# ------------------------------------------------------------

mash sketch -i -k 21 -s 1000 \
    -o enriched_vOTUs \
    enriched_vOTUs.fasta


# ------------------------------------------------------------
# 2. Compare vOTUs against INPHARED2
# ------------------------------------------------------------

mash dist \
    INPHARED2_genomes.fa.msh \
    enriched_vOTUs.msh \
    > mash_all_hits.tsv


# ------------------------------------------------------------
# 3. Retain the closest reference for each vOTU
# ------------------------------------------------------------

awk -F'\t' 'BEGIN{OFS="\t"}
{
    q=$2
    d=$3+0

    if (!(q in best) || d < best[q]) {
        best[q]=d
        line[q]=$0
    }
}
END {
    for (q in line) print line[q]
}' mash_all_hits.tsv \
| sort -t$'\t' -k2,2 \
> mash_best_hit_per_vOTU.tsv


# ------------------------------------------------------------
# 4. Identify vOTUs sharing >=10/1,000 hashes with their
#    closest INPHARED2 reference
# ------------------------------------------------------------

awk -F'\t' '{
    split($5,a,"/")
    if (a[1] >= 10) print
}' mash_best_hit_per_vOTU.tsv \
> mash_hits_ge10.tsv


# ------------------------------------------------------------
# 5. ANI analysis of selected vOTU-reference pairs
# ------------------------------------------------------------
# FASTA sequences of the selected vOTUs and their corresponding
# INPHARED2 reference genomes were combined before BLASTn.

cat selected_vOTUs.fasta \
    selected_INPHARED2_references.fasta \
    > vOTUs_plus_references.fasta


# Build nucleotide BLAST database

makeblastdb \
    -in vOTUs_plus_references.fasta \
    -dbtype nucl \
    -out ANI_db


# Perform nucleotide alignments in the format required by anicalc.py

blastn \
    -query vOTUs_plus_references.fasta \
    -db ANI_db \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -out ANI_blast.tsv


# ------------------------------------------------------------
# 6. Calculate ANI and alignment coverage
# ------------------------------------------------------------

python3 anicalc.py \
    -i ANI_blast.tsv \
    -o ANI_results.tsv
