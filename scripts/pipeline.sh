
set -e  

mkdir -p data log/cutadapt out/merged out/trimmed out/star

# (1) Download sample FASTQ if not present
while read -r url; do
    file=$(basename "$url")
    if [ -f data/"$file" ]; then
        echo "$file already exists, skipping"
    else
       wget -P data "$url"
    fi
done < data/urls

# (2) Download contaminants fasta if not present
contaminants_url="https://masterbioinformatica.com/decont/contaminants.fasta.gz"
if [ ! -f res/contaminants.fasta.gz ]; then
    wget -P res "$contaminant_url"
fi

# (3) Index the contaminnats
bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

# (4) Detect unique sample IDs
SAMPLES=$(ls data/*.fastq.gz | sed 's/.*\///' | cut -d- -f1 | sort | uniq)
echo "Detected samples: $SAMPLES"

for sid in $SAMPLES; do
    bash scripts/merge_fastqs.sh data out/merged "$sid"
done

# (5) Execute cutadapt
logfile=log/pipeline.log
echo "Pipeline run: $(date)" > "$logfile"

for sid in $SAMPLES; do
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o out/trimmed/${sid}.trimmed.fastq.gz out/merged/${sid}.fastq.gz \
        > log/cutadapt/${sid}.log
    echo "Cutadapt summary for $sid" >> "$logfile"
    grep "Total reads processed" log/cutadapt/${sid}.log >> "$logfile"
    grep "Reads with adapters" log/cutadapt/${sid}.log >> "$logfile"
done

# (6) Execute STAR aligment for eliminating contaminants
for fname in out/trimmed/*.fastq.gz; do
    sid=$(basename "$file" .trimmed.fastq.gz)
    mkdir -p out/star/"$sid"
    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
        --outReadsUnmapped Fastx --readFilesIn "$fname" \
        --readFilesCommand gunzip -c \ --outFileNamePrefix out/star/"$sid"/ \
        > out/star/$sid/STAR.log
    echo "STAR summary for $sid" >> "$logfile"
    grep ""Uniquely mapped reads %" out/star/$sid/Log.final.out >> "$logfile"
    grep "Reads mapped to multiple loci %" out/star/$sid/Log.final.out >> "$logfile"
done

echo "Pipeline finished successfully!"
