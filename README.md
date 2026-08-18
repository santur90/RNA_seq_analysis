# RNA_seq_analysis

A reproducible bulk RNA-seq workflow for adapter trimming, rRNA depletion, STAR alignment, gene-level quantification, DESeq2 differential expression, and preranked GSEA.

## Workflow

`FASTQ -> Trim Galore/FastQC -> rRNA filtering -> STAR -> featureCounts -> DESeq2 -> GSEA`

## Quick start

```bash
mamba env create -f environment.yml
mamba activate rnaseq
./run_rnaseq.sh --config config.tsv --samples samples.tsv --threads 12
```

Use `--dry-run` to inspect the commands without processing FASTQ files.

## Inputs

`samples.tsv` is tab-separated:

| sample | condition | batch | fastq_r1 | fastq_r2 |
|---|---|---|---|---|
| sample_1 | MDS | batch1 | data/sample_1_R1.fastq.gz | data/sample_1_R2.fastq.gz |
| sample_2 | Ctrl | batch1 | data/sample_2_R1.fastq.gz | data/sample_2_R2.fastq.gz |

The DESeq2 metadata table must contain `sampleName`, `condition`, and `batch` columns. The count matrix uses gene IDs as rows and sample names as columns.

## Design

The default design is `~ batch + condition`. The pipeline removes rRNA reads before genome alignment, produces coordinate-sorted indexed BAM files, and quantifies exon-level reads with featureCounts. DESeq2 exports all results, significant genes, normalized counts, and QC plots. GSEA consumes a ranked table and a user-provided GMT gene-set file.

## Outputs

- `results/qc/`: FastQC, trimming, rRNA, STAR, and MultiQC reports
- `results/bam/`: sorted/indexed genome alignments and alignment statistics
- `results/counts/`: featureCounts matrix and assignment summary
- `results/deseq2/`: differential expression tables and QC plots
- `results/gsea/`: preranked GSEA reports
- `results/logs/`: stage logs
