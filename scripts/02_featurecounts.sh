#!/usr/bin/env bash
set -euo pipefail
config=$1
bamdir=$2
outdir=$3
threads=$(awk -F '\t' '$1 == "THREADS" {print $2; exit}' "$config")
annotation=$(awk -F '\t' '$1 == "ANNOTATION_GTF" {print $2; exit}' "$config")
strandedness=$(awk -F '\t' '$1 == "STRANDEDNESS" {print $2; exit}' "$config")
mkdir -p "$outdir"
shopt -s nullglob
bams=("$bamdir"/*_Aligned.sortedByCoord.out.bam)
[[ "${#bams[@]}" -gt 0 ]] || { echo "No STAR BAM files found" >&2; exit 1; }
featureCounts -a "$annotation" -o "$outdir/gene_counts.txt" -p -B -C -T "$threads" -s "$strandedness" \
  -t exon -g gene_id "${bams[@]}" > "$outdir/featureCounts.log" 2>&1
