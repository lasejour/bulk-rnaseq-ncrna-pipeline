# bulk-rnaseq-ncrna-pipeline

An end-to-end Nextflow pipeline for bulk RNA-Seq quality assessment, adapter trimming,
genome mapping and quantification.

---

## Pipeline Overview
[Pipeline DAG](dag.png)

The pipeline handles both single-end and paired-end reads and performs:

- Initial and post-trimming QC (FastQC + MultiQC)
- Adapter trimming (Cutadapt)
- rRNA and ribosomal spacer depletion QC (STAR)
- Splice-junction-aware two-pass mapping (STAR)
- Transcript-level quantification with auto-strandedness detection (Salmon)
- Gene-level quantification across exon and gene features, stranded and unstranded (featureCounts)
- A merged summary table across all MultiQC steps

---

## Requirements
- Nextflow ≥ 25.10
- Conda (for the local profile)
- All tool versions are managed per-process in nextflow.config — no manual installs needed
- To create an envirnoment with all tools, users can use the envirnoment.yaml file 

---

## Installation
```
git clone https://github.com/lasejour/bulk-rnaseq-ncrna-pipeline.git
cd bulk-rnaseq-ncrna-pipeline
```

---

## Usage

```
nextflow run main.nf \
    --fastq_dir 'data/' \
    --pattern '*_{1,2}.fastq.gz' \
    --output my_results \
    -profile local \
    -params-file params.yml
```

or

```
nextflow run main.nf \
    --sample_file 'sample_file.csv' \
    --output my_results \
    -profile local \
    -params-file params.yml
```

For help:
```
nextflow run main.nf --help
```

---

## Input

The pipeline accepts reads in two ways:

Option 1 - Directory + pattern
```
--fastq_dir data/
--pattern '*_{1,2}.fastq.gz'
```

Option 2 - Sample sheet (CSV):
```
SampleName,Read1,Read2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

---

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--fastq_dir` | Path to input FASTQ directory | — |
| `--pattern` | Regex pattern to match FASTQ files | — |
| `--sample_file` | Path to CSV sample sheet (alternative to fastq_dir) | — |
| `--output` | Path to results directory | `results` |
| `--prep` | Read orientation: `single-end` or `paired-end` | — |
| `--genome_index` | Path to pre-built STAR index directory | — |
| `--salmon_index` | Path to pre-built Salmon index directory | — |
| `--rRNA_index` | Path to STAR rRNA genome index | — |
| `--rRNA_Spacers_index` | Path to STAR ribosomal spacers index | — |
| `--genome_annot_file` | Path to genome GTF annotation file | — |
| `--genome_fasta_file` | Path to primary assembly FASTA (for index construction) | — |
| `--transcriptome_file` | Path to transcriptome FASTA (for Salmon index construction) | — |
| `--genome_id_file` | Path to Ensembl transcript/gene ID table | — |
| `--splice_junctions` | Enable two-pass splice-junction-aware mapping | `false` |
| `--min_unique_reads` | Min uniquely mapping reads for splice junction filtering | `3` |
| `-params-file` | Path to params YAML file | `params.yml` |
| `-profile` | Execution profile (e.g. `local`) | — |

---

## Outputs

All results are written to `--output` (default: `results/`):

| Directory | Contents |
|---|---|
| `initialQC/` | FastQC and MultiQC reports on raw reads |
| `cleanedQC/` | FastQC and MultiQC reports post-trimming |
| `cutadapt/` | Trimmed FASTQ files and trimming logs |
| `salmonResults/` | Salmon quantification output and MultiQC summary |
| `genomeMultiQC/` | STAR alignment MultiQC summary |
| `riboMULTIQC/` | rRNA depletion mapping MultiQC summary |
| `spacerMULTIQC/` | Ribosomal spacer mapping MultiQC summary |
| `featureCounts/` | Raw counts tables (exon/gene × stranded/unstranded) |
| `finalMultiQC/` | Combined stranded and unstranded final MultiQC reports |
| `mergedSummary/` | Single merged summary table across all MultiQC steps |
| `pipeline_info/` | Execution timeline, report, trace, and DAG |

---

## Tools & Versions

| Tool | Version |
|---|---|
| FastQC | 0.12.1 |
| Cutadapt | 5.2 |
| STAR | 2.7.11b |
| Salmon | 1.10.2 |
| featureCounts (Subread) | 2.0.6 |
| MultiQC | 1.32 |

---


## Citation

If you use this pipeline in your work, please cite:

> Sejour L. bulk-rnaseq-ncrna-pipeline (v2.0.0). Precision RNA Medicine Core,
> Beth Israel Deaconess Medical Center. https://github.com/lasejour/bulk-rnaseq-ncrna-pipeline

---

## Author

**Leinal Sejour**
Bioinformatician I, Precision RNA Medicine Core
Beth Israel Deaconess Medical Center / Harvard Medical School
[GitHub](https://github.com/lasejour)
