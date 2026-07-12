# FASTQ Genomic Processing Pipeline & 23andMe Visualizer

A lightweight, local genomic data analysis suite for processing raw sequencer files (FASTQ) and mapping personal genotyping records (23andMe) into an interactive health-impact dashboard.

---

## 🧬 What this Project Does

1. **DNA Sequence Alignment & Variant Calling (FASTQ)**:
   Takes raw, paired-end short reads (`FASTQ` format) from a sequencing machine, trims adapters/quality controls via `fastp`, aligns them to a reference genome via `bowtie2`, processes them using `samtools`, and identifies SNPs and Indels using `bcftools`, generating a final Variant Call Format (`VCF`) file.
2. **Personal Genetic Health Mapping (23andMe)**:
   Parses a raw 23andMe microarray genotyping text file, resolves genomic positions (GRCh38) and consequence categories (e.g., *missense variant*, *stop gained*) in batches via the **Ensembl VEP REST API**, and maps variants to an interactive visual HTML dashboard showing category-level genetic loads.

---

## 🛠️ Toolchain Installation (macOS)

Ensure you have Homebrew installed, then run:
```bash
brew install fastp bowtie2 samtools bcftools wget
```

---

## 🚀 How to Run the Code

### 1. FASTQ Sequence Processing Pipeline
We use a self-contained, synthetic *E. coli* sandbox dataset to demonstrate 100% variant-calling accuracy:

```bash
# Step A: Generate synthetic reference FASTA and paired-end FASTQ reads
python3 scratch/generate_toy_data.py

# Step B: Run the automated pipeline (Quality Control -> Indexing -> Mapping -> Call Variants)
bash run_pipeline.sh
```
*Output variant results are saved to `results/variants.vcf` and HTML quality control reports are saved to `results/fastp_report.html`.*

### 2. 23andMe Gene-to-Feature Visualizer
Downloads a sample public, consented 23andMe file from the Personal Genome Project (PGP) and compiles a dashboard:

```bash
# Step A: Fetch dataset, query Ensembl REST API, and generate the CSV matrices
python3 scratch/generic_gene_feature_matrix.py

# Step B: Inject JSON records into the self-contained HTML visual dashboard
python3 scratch/generate_html_report.py
```
*Output data matrices are saved to `results/generic_gene_features.csv` and the interactive dashboard is saved to `results/genomic_report.html`.*

---

## 📊 What the HTML Dashboard Shows

Open `results/genomic_report.html` in any web browser to view:
* **Interactive Radar Chart**: Displays the individual's average genetic risk load/density across 8 primary health categories (Methylation, Brain, Cardiovascular, Dietary, Vitamins, Detox, Inflammation, and Physical Traits).
* **Circular Progress Gauges**: Circular risk meters mapping the average risk score (0.0 to 1.0) of 24 distinct genes.
* **Real-time Filter & Search**: Instant filter tabs by category and reactive text searches.
* **Consequence Modal Details**: Clicking any card shows precise SNP annotations, including genomic positions (GRCh38), Ensembl consequence tags, and patient-specific clinical interpretations.
