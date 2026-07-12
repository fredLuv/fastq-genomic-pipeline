#!/bin/bash
# FASTQ Genomic Processing Pipeline
# Exit immediately if any command fails
set -e

# Setup directories
WORKDIR="/Users/fred/Library/CloudStorage/GoogleDrive-fredlu0521@gmail.com/My Drive/Ethan-Resume/fastq_pipeline"
REF_DIR="$WORKDIR/ref"
READS_DIR="$WORKDIR/reads"
RESULTS_DIR="$WORKDIR/results"

mkdir -p "$RESULTS_DIR"

echo "=========================================================="
echo " Starting FASTQ Genomic Data Processing Pipeline"
echo "=========================================================="
echo "Working directory: $WORKDIR"

# Step 1: Quality Control and Adapter Trimming via fastp
echo ""
echo "=== Step 1: Quality Control and Trimming (fastp) ==="
fastp \
  -i "$READS_DIR/read1.fastq" \
  -I "$READS_DIR/read2.fastq" \
  -o "$RESULTS_DIR/clean1.fastq" \
  -O "$RESULTS_DIR/clean2.fastq" \
  -h "$RESULTS_DIR/fastp_report.html" \
  -j "$RESULTS_DIR/fastp_report.json" \
  --detect_adapter_for_pe \
  --thread 4

echo "Trimming completed. Reports saved to:"
echo " - HTML report: $RESULTS_DIR/fastp_report.html"
echo " - JSON report: $RESULTS_DIR/fastp_report.json"

# Step 2: Build Reference Genome Index via Bowtie2
echo ""
echo "=== Step 2: Building Index for Reference Genome ==="
bowtie2-build \
  --quiet \
  "$REF_DIR/ref.fa" \
  "$REF_DIR/ecoli_index"

echo "Index built successfully at $REF_DIR/ecoli_index"

# Step 3: Align Paired-End Reads via Bowtie2
echo ""
echo "=== Step 3: Aligning Reads to Reference (Bowtie2) ==="
bowtie2 \
  -x "$REF_DIR/ecoli_index" \
  -1 "$RESULTS_DIR/clean1.fastq" \
  -2 "$RESULTS_DIR/clean2.fastq" \
  -S "$RESULTS_DIR/aligned.sam" \
  --threads 4

echo "Alignment completed. SAM output saved to $RESULTS_DIR/aligned.sam"

# Step 4: Convert SAM to sorted BAM and create index via Samtools
echo ""
echo "=== Step 4: Processing Alignments (Samtools) ==="
echo "Converting SAM to BAM..."
samtools view -S -b "$RESULTS_DIR/aligned.sam" > "$RESULTS_DIR/aligned.bam"

echo "Sorting BAM file..."
samtools sort "$RESULTS_DIR/aligned.bam" -o "$RESULTS_DIR/aligned.sorted.bam"

echo "Indexing BAM file..."
samtools index "$RESULTS_DIR/aligned.sorted.bam"

# Clean up intermediate SAM and unsorted BAM files to save space
rm -f "$RESULTS_DIR/aligned.sam" "$RESULTS_DIR/aligned.bam"

echo "Alignment processing completed. Sorted BAM: $RESULTS_DIR/aligned.sorted.bam"

# Step 5: Pileup and Variant Calling via bcftools
echo ""
echo "=== Step 5: Variant Calling (bcftools) ==="
# Generate genotype likelihoods (mpileup) and call variants
bcftools mpileup \
  -Ou \
  -f "$REF_DIR/ref.fa" \
  "$RESULTS_DIR/aligned.sorted.bam" | \
bcftools call \
  -mv \
  -Ob \
  -o "$RESULTS_DIR/variants.bcf"

# Convert BCF to human-readable VCF format
bcftools view \
  "$RESULTS_DIR/variants.bcf" > "$RESULTS_DIR/variants.vcf"

echo "Variant calling completed. VCF results: $RESULTS_DIR/variants.vcf"

# Step 6: Summary of detected variants
echo ""
echo "=========================================================="
echo " Pipeline Complete. Summary of Detected Variants:"
echo "=========================================================="
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$RESULTS_DIR/variants.vcf"
echo "=========================================================="
