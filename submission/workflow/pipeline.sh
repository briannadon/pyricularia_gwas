#!/bin/bash

# GWAS Pipeline Script for Pyricularia oryzae
# This script performs quality assessment, alignment, and statistics for paired-end reads

set -e  # Exit on any error

# Configuration
DATA_DIR="../data"
REFERENCE="GCA_002368485.1_ASM236848v1_genomic.fna"
OUTPUT_DIR="results"
LOG_DIR="logs"
THREADS=8
SORT_THREADS=4

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$OUTPUT_DIR/quality_stats"
mkdir -p "$OUTPUT_DIR/alignments"
mkdir -p "$OUTPUT_DIR/alignment_stats"
mkdir -p "$OUTPUT_DIR/variants"

# Function to log messages with timestamps
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/pipeline.log"
}

# Function to check if BWA index exists
check_bwa_index() {
    local ref_base="${DATA_DIR}/${REFERENCE}"
    log "Checking for BWA index files in: $ref_base"
    log "Current working directory: $(pwd)"
    log "Checking for files:"
    log "  ${ref_base}.amb - $(if [[ -f "${ref_base}.amb" ]]; then echo "EXISTS"; else echo "MISSING"; fi)"
    log "  ${ref_base}.ann - $(if [[ -f "${ref_base}.ann" ]]; then echo "EXISTS"; else echo "MISSING"; fi)"
    log "  ${ref_base}.bwt - $(if [[ -f "${ref_base}.bwt" ]]; then echo "EXISTS"; else echo "MISSING"; fi)"
    log "  ${ref_base}.pac - $(if [[ -f "${ref_base}.pac" ]]; then echo "EXISTS"; else echo "MISSING"; fi)"
    log "  ${ref_base}.sa - $(if [[ -f "${ref_base}.sa" ]]; then echo "EXISTS"; else echo "MISSING"; fi)"
    
    if [[ -f "${ref_base}.amb" && -f "${ref_base}.ann" && -f "${ref_base}.bwt" && -f "${ref_base}.pac" && -f "${ref_base}.sa" ]]; then
        log "BWA index found for $REFERENCE"
        return 0
    else
        log "BWA index not found for $REFERENCE. Creating index..."
        return 1
    fi
}

# Function to create BWA index
create_bwa_index() {
    log "Creating BWA index for $REFERENCE"
    cd "$DATA_DIR"
    source /opt/conda/etc/profile.d/conda.sh
    conda activate bwa
    bwa index "$REFERENCE"
    conda deactivate
    cd - > /dev/null
    log "BWA index created successfully"
}

# Function to get quality statistics using FastQC
run_quality_stats() {
    local sample=$1
    local r1="$OUTPUT_DIR/alignments/${sample}_1_downsampled.fastq"
    local r2="$OUTPUT_DIR/alignments/${sample}_2_downsampled.fastq"
    local fastqc_dir="$OUTPUT_DIR/quality_stats"
    local summary_file="$fastqc_dir/${sample}_summary.txt"
    
    # Check if quality stats already exist
    if [[ -f "$summary_file" ]]; then
        log "Quality statistics already exist for $sample, skipping..."
        return 0
    fi
    
    log "Running quality statistics for sample $sample"
    
    # Run FastQC with conda environment activation
    source /opt/conda/etc/profile.d/conda.sh
    conda activate fastqc
    fastqc -o "$fastqc_dir" -t "$THREADS" "$r1" "$r2"
    conda deactivate
    
    # Generate summary statistics
    echo "Sample: $sample" > "$summary_file"
    echo "Read 1: $(wc -l < "$r1" | awk '{print $1/4}') reads" >> "$summary_file"
    echo "Read 2: $(wc -l < "$r2" | awk '{print $1/4}') reads" >> "$summary_file"
    
    log "Quality statistics completed for $sample"
}

# Function to downsample reads using seqtk
downsample_reads() {
    local sample=$1
    local r1="${DATA_DIR}/${sample}_1.fastq.gz"
    local r2="${DATA_DIR}/${sample}_2.fastq.gz"
    local downsample_r1="$OUTPUT_DIR/alignments/${sample}_1_downsampled.fastq"
    local downsample_r2="$OUTPUT_DIR/alignments/${sample}_2_downsampled.fastq"
    
    # Check if downsampled files already exist
    if [[ -f "$downsample_r1" && -f "$downsample_r2" ]]; then
        log "Downsampled files already exist for $sample, skipping..."
        return 0
    fi
    
    log "Downsampling reads for sample $sample (seed=100, 10% of reads)"
    
    # Activate seqtk environment
    source /opt/conda/etc/profile.d/conda.sh
    conda activate seqtk
    
    # Downsample to 10% of reads with seed=100 (no gzipping)
    seqtk sample -s 100 "$r1" 0.1 > "$downsample_r1"
    seqtk sample -s 100 "$r2" 0.1 > "$downsample_r2"
    
    conda deactivate
    
    log "Downsampling completed for $sample"
}

# Function to align reads using BWA
run_alignment() {
    local sample=$1
    local downsample_r1="$OUTPUT_DIR/alignments/${sample}_1_downsampled.fastq"
    local downsample_r2="$OUTPUT_DIR/alignments/${sample}_2_downsampled.fastq"
    local ref_path="${DATA_DIR}/${REFERENCE}"
    local bam_output="$OUTPUT_DIR/alignments/${sample}.bam"
    local bai_output="$OUTPUT_DIR/alignments/${sample}.bam.bai"
    
    # Check if BAM file and index already exist
    if [[ -f "$bam_output" && -f "$bai_output" ]]; then
        log "BAM file and index already exist for $sample, skipping alignment..."
        return 0
    fi
    
    log "Aligning downsampled reads for sample $sample to reference"
    
    # Activate conda environments for BWA and samtools
    source /opt/conda/etc/profile.d/conda.sh
    
    # Run BWA mem alignment and save to temporary SAM file
    conda activate bwa
    bwa mem -t "$THREADS" "$ref_path" "$downsample_r1" "$downsample_r2" > "$OUTPUT_DIR/alignments/${sample}_temp.sam"
    conda deactivate
    
    # Process with samtools
    conda activate samtools
    samtools view -@ "$THREADS" -bS "$OUTPUT_DIR/alignments/${sample}_temp.sam" | \
    samtools sort -@ "$SORT_THREADS" -m 2G -T "$OUTPUT_DIR/alignments/tmp_${sample}" -o "$bam_output" -
    
    # Index the BAM file
    samtools index "$bam_output"
    
    # Clean up temporary SAM file
    rm "$OUTPUT_DIR/alignments/${sample}_temp.sam"
    conda deactivate
    
    log "Alignment completed for $sample"
}

# Function to get alignment statistics
run_alignment_stats() {
    local sample=$1
    local bam_file="$OUTPUT_DIR/alignments/${sample}.bam"
    local stats_file="$OUTPUT_DIR/alignment_stats/${sample}_stats.txt"
    
    # Check if stats already exist
    if [[ -f "$stats_file" ]]; then
        log "Alignment statistics already exist for $sample, skipping..."
        return 0
    fi
    
    log "Generating alignment statistics for $sample"
    
    # Activate samtools environment
    source /opt/conda/etc/profile.d/conda.sh
    conda activate samtools
    
    # Get basic alignment statistics
    samtools flagstat "$bam_file" > "$stats_file"
    
    # Get basic coverage information
    echo "" >> "$stats_file"
    echo "=== Basic Statistics ===" >> "$stats_file"
    echo "Total reads: $(samtools view -c "$bam_file")" >> "$stats_file"
    echo "Mapped reads: $(samtools view -c -F 4 "$bam_file")" >> "$stats_file"
    echo "Unmapped reads: $(samtools view -c -f 4 "$bam_file")" >> "$stats_file"
    
    # Get mapping quality distribution
    echo "" >> "$stats_file"
    echo "=== Mapping Quality Distribution ===" >> "$stats_file"
    samtools view "$bam_file" | cut -f5 | sort -n | uniq -c | sort -k2 -n >> "$stats_file"
    
    conda deactivate
    log "Alignment statistics completed for $sample"
}

# Function to call variants using samtools mpileup on all samples
run_variant_calling() {
    local vcf_output="$OUTPUT_DIR/variants/all_samples.vcf"
    local ref_path="${DATA_DIR}/${REFERENCE}"
    
    # Check if VCF file already exists
    if [[ -f "$vcf_output" ]]; then
        log "VCF file already exists, skipping variant calling..."
        return 0
    fi
    
    log "Calling variants for all samples"
    
    # Get list of all BAM files
    local bam_files=""
    for sample in $samples; do
        local bam_file="$OUTPUT_DIR/alignments/${sample}.bam"
        if [[ -f "$bam_file" ]]; then
            bam_files="$bam_files $bam_file"
        else
            log "Warning: BAM file not found for sample $sample, skipping..."
        fi
    done
    
    if [[ -z "$bam_files" ]]; then
        log "Error: No BAM files found for variant calling"
        return 1
    fi
    
    log "Found $(echo $bam_files | wc -w) BAM files for variant calling"
    
    # Activate samtools environment
    source /opt/conda/etc/profile.d/conda.sh
    conda activate samtools
    
    # Call variants using samtools mpileup on all BAM files
    samtools mpileup -f "$ref_path" -v -u $bam_files > "$vcf_output"
    
    conda deactivate
    log "Variant calling completed for all samples"
}

# Main pipeline execution
main() {
    log "Starting GWAS pipeline"
    
    # Check and create BWA index if needed
    if ! check_bwa_index; then
        create_bwa_index
    fi
    
    # Get list of all paired samples
    cd "$DATA_DIR"
    samples=$(ls *_1.fastq.gz | sed 's/_1\.fastq\.gz$//' | sort)
    cd - > /dev/null
    
    log "Found $(echo "$samples" | wc -l) paired samples"
    
    # Process each sample
    for sample in $samples; do
        log "Processing sample: $sample"
        
        # Downsample reads
        downsample_reads "$sample"

        # Run quality statistics
        run_quality_stats "$sample"   
        
        # Run alignment on downsampled reads
        run_alignment "$sample"
        
        # Run alignment statistics
        run_alignment_stats "$sample"
        
        log "Completed processing sample: $sample"
    done
    
    # Run variant calling on all samples
    run_variant_calling
    
    # Generate summary report
    log "Generating summary report"
    echo "GWAS Pipeline Summary Report" > "$OUTPUT_DIR/pipeline_summary.txt"
    echo "Generated on: $(date)" >> "$OUTPUT_DIR/pipeline_summary.txt"
    echo "Total samples processed: $(echo "$samples" | wc -l)" >> "$OUTPUT_DIR/pipeline_summary.txt"
    echo "" >> "$OUTPUT_DIR/pipeline_summary.txt"
    
    # Add sample-specific summaries
    for sample in $samples; do
        echo "=== $sample ===" >> "$OUTPUT_DIR/pipeline_summary.txt"
        if [[ -f "$OUTPUT_DIR/quality_stats/${sample}_summary.txt" ]]; then
            cat "$OUTPUT_DIR/quality_stats/${sample}_summary.txt" >> "$OUTPUT_DIR/pipeline_summary.txt"
        fi
        if [[ -f "$OUTPUT_DIR/alignment_stats/${sample}_stats.txt" ]]; then
            echo "Alignment stats available in: ${sample}_stats.txt" >> "$OUTPUT_DIR/pipeline_summary.txt"
        fi
    # Add variant calling summary
    if [[ -f "$OUTPUT_DIR/variants/all_samples.vcf" ]]; then
        echo "" >> "$OUTPUT_DIR/pipeline_summary.txt"
        echo "=== Variant Calling ===" >> "$OUTPUT_DIR/pipeline_summary.txt"
        echo "Combined variants available in: all_samples.vcf" >> "$OUTPUT_DIR/pipeline_summary.txt"
        echo "Total variants: $(grep -v '^#' "$OUTPUT_DIR/variants/all_samples.vcf" | wc -l)" >> "$OUTPUT_DIR/pipeline_summary.txt"
    fi
        echo "" >> "$OUTPUT_DIR/pipeline_summary.txt"
    done
    
    log "Pipeline completed successfully!"
    log "Results available in: $OUTPUT_DIR"
    log "Logs available in: $LOG_DIR"
}

# Check if required tools are available
check_dependencies() {
    log "Checking conda environments for required tools..."
    
    # Check if conda is available
    if ! command -v conda &> /dev/null; then
        echo "Error: conda is not available. Please ensure conda is installed and in PATH"
        exit 1
    fi
    
    # Check if required environments exist
    source /opt/conda/etc/profile.d/conda.sh
    local environments=("bwa" "samtools" "fastqc" "seqtk")
    for env in "${environments[@]}"; do
        if ! conda env list | grep -q "^$env "; then
            echo "Error: conda environment '$env' not found"
            exit 1
        fi
    done
    
    log "All required conda environments are available"
}

# Run dependency check and main pipeline
check_dependencies
main
