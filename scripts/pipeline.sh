#!/bin/bash
set -e  # para que falle si hay algún error

# 1️⃣ Descargar todos los FASTQs
while read -r url; do
    bash scripts/download.sh "$url" data
done < data/urls

# 2️⃣ Descargar el FASTA de contaminantes y filtrarlo si es necesario
CONTAMINANTS_URL="https://masterbioinformatica.com/decont/contaminants.fasta.gz"
bash scripts/download.sh "$CONTAMINANTS_URL" res yes

# 3️⃣ Indexar el FASTA de contaminantes
bash scripts/index.sh res/contaminants.fasta res/contaminants_idx

# 4️⃣ Obtener IDs de las muestras y mergear los FASTQs técnicos
SAMPLES=$(ls data/*.fastq.gz | sed 's/.*\///' | cut -d- -f1 | sort | uniq)

for sid in $SAMPLES; do
    bash scripts/merge_fastqs.sh data out/merged "$sid"
done

# 5️⃣ Ejecutar cutadapt para recortar adaptadores
mkdir -p out/trimmed log/cutadapt

for file in out/merged/*.fastq.gz; do
    sid=$(basename "$file" .fastq.gz)
    cutadapt -m 18 -a TGGAATTCTCGGGTGCCAAGG --discard-untrimmed \
        -o out/trimmed/"$sid".trimmed.fastq.gz "$file" > log/cutadapt/"$sid".log
done

# 6️⃣ Ejecutar STAR para eliminar contaminantes
for file in out/trimmed/*.fastq.gz; do
    sid=$(basename "$file" .trimmed.fastq.gz)
    mkdir -p out/star/"$sid"
    STAR --runThreadN 4 --genomeDir res/contaminants_idx \
        --outReadsUnmapped Fastx --readFilesIn "$file" \
        --readFilesCommand zcat --outFileNamePrefix out/star/"$sid"/
done

# 7️⃣ Crear log final
mkdir -p log
> log/pipeline.log
for sid in $SAMPLES; do
    echo "Sample: $sid" >> log/pipeline.log
    grep "Reads with adapters" log/cutadapt/"$sid".log >> log/pipeline.log || true
    grep "Uniquely mapped reads %" out/star/"$sid"/Log.final.out >> log/pipeline.log || true
    echo "" >> log/pipeline.log
done

echo "Pipeline finished successfully!"
