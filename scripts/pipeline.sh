#!/bin/bash
set -e  # Stop execution if any command fails

(1) Create necessary directories
mkdir -p data log log/cutadapt out/merged out/trimmed out/star res

(2) Download sample FASTQ files
echo "Downloading sample FASTQs..."
while read -r url; do
    file=$(basename "$url")
    if [ -f data/"$file" ]; then
        echo "$file already exists, skipping"
    else
        wget -P data "$url"
        # Optional: MD5 check
        if [ -f "$url.md5" ]; then
            md5sum -c "$url.md5"
        fi
    fi
done < data/urls

(3) Download contaminants database
contaminants_url="https://bioinformatics.cnio.es/data/courses/decont/contaminants.fasta.gz"
if [ ! -f res/contaminants.fasta.gz ]; then
    echo "Downloading contaminants database..."
    wget -P res "$contaminants_url"
fi

(4) Remove small nuclear RNAs from contaminants (snRNA)
if [ ! -f res/contaminants.fasta ]; then
    gunzip -c res/contaminants.fasta.gz | \
    seqkit grep -v -p "small nuclear" > res/contaminants.fasta
fi

(5) Index contaminants database using STAR
bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

(6) Detect sample IDs automatically
SAMPLES=$(ls data/*.fastq.gz | sed 's/.*\///' | cut -d- -f1 | sort | uniq)
echo "Detected samples: $SAMPLES"

(7) Merge technical replicates
for sid in $SAMPLES; do
    bash scripts/merge_fastqs.sh data out/merged "$sid"
done

(8) Initialize pipeline log
pipeline_log="log/pipeline.log"
echo "Pipeline run: $(date)" > "$pipeline_log"

(9) Trim adapters with cutadapt
for sid in $SAMPLES; do
    trimmed_file="out/trimmed/${sid}.trimmed.fastq.gz"
    if [ -f "$trimmed_file" ]; then
        echo "$sid already trimmed, skipping" >> "$pipeline_log"
        continue
    fi

    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o "$trimmed_file" out/merged/${sid}.fastq.gz \
        > log/cutadapt/${sid}.log

    echo "Cutadapt summary for $sid" >> "$pipeline_log"
    grep "Total reads processed" log/cutadapt/${sid}.log >> "$pipeline_log"
    grep "Reads with adapters" log/cutadapt/${sid}.log >> "$pipeline_log"
done

(10) STAR alignment to remove contaminant reads
for trimmed_file in out/trimmed/*.fastq.gz; do
    sid=$(basename "$trimmed_file" .trimmed.fastq.gz)
    star_out="out/star/$sid"
    mkdir -p "$star_out"

    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
        --outReadsUnmapped Fastx --readFilesIn "$trimmed_file" \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix "$star_out/" \
        > "$star_out/STAR.log"

    echo "STAR summary for $sid" >> "$pipeline_log"
    grep "Uniquely mapped reads %" "$star_out/Log.final.out" >> "$pipeline_log"
    grep "Reads mapped to multiple loci %" "$star_out/Log.final.out" >> "$pipeline_log"
done

(11) Pipeline finished
echo "Pipeline finished successfully!" >> "$pipeline_log"
