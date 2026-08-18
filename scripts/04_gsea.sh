#!/usr/bin/env bash
set -euo pipefail
input=${1:?Provide DESeq2 result table}
gmt=${2:?Provide GMT gene-set file}
outdir=${3:-results/gsea}
label=${4:-rnaseq}
mkdir -p "$outdir"
rank="$outdir/${label}.rnk"
awk -F ',' 'NR > 1 && $1 != "" && $3 != "NA" {g=$1; sub(/\..*/, "", g); print g "\t" $3}' "$input" | sort -k2,2gr -u > "$rank"
[[ -s "$rank" ]] || { echo "Rank file is empty; expected gene ID and log2FoldChange columns" >&2; exit 1; }
command -v gsea-cli.sh >/dev/null 2>&1 || { echo "gsea-cli.sh is required for this stage" >&2; exit 1; }
gsea-cli.sh GSEAPreranked -rnk "$rank" -gmx "$gmt" -rpt_label "$label" -out "$outdir" -nperm 1000 -set_min 15 -set_max 500 -collapse No_Collapse -zip_report false
