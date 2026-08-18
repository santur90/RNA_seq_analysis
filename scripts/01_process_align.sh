#!/usr/bin/env bash
set -euo pipefail

samples=$1
config=$2
outdir=$3
threads=$4
get_config() { awk -F '\t' -v key="$1" '$1 == key {print $2; exit}' "$config"; }
trim_threads=$(get_config TRIM_THREADS)
rrna_index=$(get_config RRNA_INDEX)
star_index=$(get_config STAR_INDEX)
min_length=$(get_config MIN_LENGTH)
quality_cutoff=$(get_config QUALITY_CUTOFF)
mkdir -p "$outdir/qc" "$outdir/trimmed" "$outdir/rrna" "$outdir/bam" "$outdir/logs"
while IFS=$'\t' read -r sample condition batch r1 r2; do
  [[ "$sample" == "sample" || -z "$sample" ]] && continue
  [[ -f "$r1" && -f "$r2" ]] || { echo "FASTQ files missing for $sample" >&2; exit 1; }
  trim_galore --paired --gzip --fastqc --cores "$trim_threads" --quality "$quality_cutoff" --length "$min_length" \
    -o "$outdir/trimmed" "$r1" "$r2" > "$outdir/logs/${sample}.trim.log" 2>&1
  r1t="$outdir/trimmed/${sample}_R1_val_1.fq.gz"
  r2t="$outdir/trimmed/${sample}_R2_val_2.fq.gz"
  bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p "$threads" -x "$rrna_index" \
    -1 "$r1t" -2 "$r2t" --un-conc-gz "$outdir/rrna/${sample}_rRNA_removed.fq.gz" \
    2> "$outdir/logs/${sample}.rrna.log" | samtools view -@ "$threads" -b -o "$outdir/rrna/${sample}_rrna.bam" -
  STAR --runThreadN "$threads" --genomeDir "$star_index" \
    --readFilesIn "$outdir/rrna/${sample}_rRNA_removed.fq.1.gz" "$outdir/rrna/${sample}_rRNA_removed.fq.2.gz" \
    --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --outSAMattributes Standard --outFileNamePrefix "$outdir/bam/${sample}_" \
    > "$outdir/logs/${sample}.STAR.log" 2>&1
  bam="$outdir/bam/${sample}_Aligned.sortedByCoord.out.bam"
  samtools index -@ "$threads" "$bam"
  samtools flagstat -@ "$threads" "$bam" > "$outdir/bam/${sample}.flagstat.txt"
  samtools idxstats "$bam" > "$outdir/bam/${sample}.idxstats.txt"
done < "$samples"
