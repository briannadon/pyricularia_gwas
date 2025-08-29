#!/bin/bash

# Script to answer questions from the GWAS pipeline results
# Questions:
# Q1: What is the modal GC content of the samples?
# Q2: What is the total number of mismatches?
# Q3: What percent of reads uniquely map to the reference genome?
# Q4: How many total filtered variants were found?
# Q5: How many variant sites have at least 1 alternate allele call?

# Set the directories
QUALITY_STATS_DIR="submission/workflow/results/quality_stats"
ALIGNMENT_STATS_DIR="submission/workflow/results/alignment_stats"
ALIGNMENTS_DIR="submission/workflow/results/alignments"
VARIANTS_FILE="submission/workflow/results/variants/all_samples.vcf"

# Function to check if required tools are available
check_dependencies() {
    echo "Checking conda environments for required tools..."
    
    # Check if conda is available
    if ! command -v conda &> /dev/null; then
        echo "Error: conda is not available. Please ensure conda is installed and in PATH"
        exit 1
    fi
    
    # Check if required environments exist
    source /opt/conda/etc/profile.d/conda.sh
    local environments=("samtools")
    for env in "${environments[@]}"; do
        if ! conda env list | grep -q "^$env "; then
            echo "Error: conda environment '$env' not found"
            exit 1
        fi
    done
    
    echo "All required conda environments are available"
}

# Check dependencies before running
check_dependencies

echo "GWAS Pipeline Question Answers"
echo "=============================="
echo ""

# Question 1: Poor Quality Sequences
echo "Question 1: How many total sequences were flagged for being poor quality?"
echo "------------------------------------------------------------------------"

# Check if directory exists
if [ ! -d "$QUALITY_STATS_DIR" ]; then
    echo "Error: Directory $QUALITY_STATS_DIR not found"
    exit 1
fi

total_poor_quality=0
sample_count=0

echo "Analyzing FastQC results for poor quality sequences..."

# Find all HTML files and extract poor quality information
for html_file in "$QUALITY_STATS_DIR"/*.html; do
    if [ -f "$html_file" ]; then
        sample_name=$(basename "$html_file" .html)
        echo "Processing $sample_name..."

        # Extract the number of sequences flagged as poor quality from the FastQC HTML
        poor_quality=$(grep -o 'Sequences flagged as poor quality</td><td>[0-9]*' "$html_file" | sed 's/Sequences flagged as poor quality<\/td><td//;s/>//')
        
        # If not found, set to 0
        if [ -z "$poor_quality" ]; then
            poor_quality=0
        fi

        if [ "$poor_quality" -gt 0 ]; then
            echo "  $sample_name: $poor_quality poor quality sequences"
            total_poor_quality=$((total_poor_quality + poor_quality))
            sample_count=$((sample_count + 1))
        else
            echo "  $sample_name: No poor quality sequences detected"
        fi
    fi
done

echo "Total sequences flagged for poor quality: $total_poor_quality"
echo ""

# Question 2: Total reads with quality score less than 20
echo "Question 2: What is the total number of reads with quality score less than 20?"
echo "---------------------------------------------------------------------------"

if [ ! -d "$ALIGNMENTS_DIR" ]; then
    echo "Error: Directory $ALIGNMENTS_DIR not found"
    exit 1
fi

total_low_quality_reads=0
sample_count=0

echo "Calculating total reads with quality score less than 20..."

# Activate samtools environment
source /opt/conda/etc/profile.d/conda.sh
conda activate samtools

for bam_file in "$ALIGNMENTS_DIR"/*.bam; do
    if [ -f "$bam_file" ] && [[ "$bam_file" != *.bai ]]; then
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing $sample_name..."
        
        # Count reads with mapping quality less than 20
        low_quality_reads=$(samtools view "$bam_file" | awk -F'\t' '$5 < 20 {count++} END {print count}')
        
        if [ ! -z "$low_quality_reads" ] && [ "$low_quality_reads" -gt 0 ]; then
            echo "  $sample_name: $low_quality_reads reads with quality < 20"
            total_low_quality_reads=$((total_low_quality_reads + low_quality_reads))
            sample_count=$((sample_count + 1))
        fi
    fi
done

conda deactivate

echo "Total reads with quality score less than 20: $total_low_quality_reads"
echo ""

# Question 3: Percentage of reads uniquely mapped
echo "Question 3: What percent of reads uniquely map to the reference genome?"
echo "----------------------------------------------------------------------"

total_reads=0
uniquely_mapped_reads=0

echo "Calculating unique mapping percentages..."

# Activate samtools environment
source /opt/conda/etc/profile.d/conda.sh
conda activate samtools

for bam_file in "$ALIGNMENTS_DIR"/*.bam; do
    if [ -f "$bam_file" ] && [[ "$bam_file" != *.bai ]]; then
        sample_name=$(basename "$bam_file" .bam)
        echo "Processing $sample_name..."
        
        # Get total reads and uniquely mapped reads
        sample_stats=$(samtools flagstat "$bam_file")
        total_sample_reads=$(echo "$sample_stats" | grep "in total" | awk '{print $1}')
        mapped_reads=$(echo "$sample_stats" | grep "mapped" | head -1 | awk '{print $1}')
        
        # Count uniquely mapped reads (mapping quality > 0)
        unique_reads=$(samtools view -q 1 "$bam_file" | wc -l)
        
        echo "  $sample_name: $unique_reads uniquely mapped out of $total_sample_reads total reads"
        
        total_reads=$((total_reads + total_sample_reads))
        uniquely_mapped_reads=$((uniquely_mapped_reads + unique_reads))
    fi
done

conda deactivate

if [ $total_reads -gt 0 ]; then
    unique_percentage=$(awk "BEGIN {printf \"%.2f\", $uniquely_mapped_reads * 100 / $total_reads}")
    echo "Total uniquely mapped reads: $uniquely_mapped_reads out of $total_reads total reads"
    echo "Percentage of reads uniquely mapped: ${unique_percentage}%"
else
    echo "No reads found in BAM files."
fi
echo ""

# Question 4: Total variants
echo "Question 4: How many total variants were found?"
echo "-----------------------------------------------"

if [ ! -f "$VARIANTS_FILE" ]; then
    echo "Error: VCF file $VARIANTS_FILE not found"
    exit 1
fi

echo "Counting total variants..."

# Count all variants (excluding header lines)
total_variants=$(grep -v "^##" "$VARIANTS_FILE" | wc -l)

echo "Total variants: $total_variants"
echo ""

# Question 5: Variant sites with at least 1 alternate allele call
echo "Question 5: How many variant sites have at least 1 alternate allele call?"
echo "------------------------------------------------------------------------"

echo "Counting variant sites with alternate alleles..."

# Count variants that have alternate alleles (not just reference calls)
# Look for variants where at least one sample has a non-reference genotype
variant_sites_with_alt=$(grep -v "^##" "$VARIANTS_FILE" | awk '
{
    # Check if any sample has a non-reference genotype
    has_alt = 0
    for (i=10; i<=NF; i++) {
        # Skip if field is missing
        if ($i != "." && $i != "0/0" && $i != "0|0") {
            has_alt = 1
            break
        }
    }
    if (has_alt) count++
}
END {print count}
')

echo "Variant sites with at least 1 alternate allele call: $variant_sites_with_alt"
echo ""

echo "=============================="
echo "Summary of Answers:"
echo "Q1. Total sequences flagged for poor quality: $total_poor_quality"
echo "Q2. Total reads with quality < 20: $total_low_quality_reads"
echo "Q3. Percentage uniquely mapped: ${unique_percentage}%"
echo "Q4. Total variants: $total_variants"
echo "Q5. Variant sites with alt alleles: $variant_sites_with_alt"
echo "=============================="
