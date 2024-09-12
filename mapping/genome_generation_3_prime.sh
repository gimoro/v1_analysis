##!/bin/bash
## Remapping v1 experiments with STARSolo
## Izaskun Mallona & Giulia Moro
##
## December 2023
## 

NTHREADS=30
export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

## genome generation for the 44 nt capture 

ID=mouse_LTRs_egfp_3_prime

WD=/home/gmoro/star_solo_v1_experiments

mkdir -p "$WD"/"$ID"

cd "$WD"/"$ID"

NAMES=(mouse alien)

## ln the files for the genome

ln -s /home/gmoro/indices/v1_mouse/GRCm38.p6.genome.fa genome.fa
ln -s /home/gmoro/indices/v1_mouse/gencode.vM25.annotation.gtf genome.gtf
ln -s /home/gmoro/indices/v1_mouse/alien_WPRE.fa alien.fa
ln -s /home/gmoro/indices/v1_mouse/3_prime_egfp_LTRs_egfp_capture.gtf alien.gtf

declare -A FAS
FAS=([mouse]=genome.fa [alien]=alien.fa)

declare -A GTFS
GTFS=([mouse]=genome.gtf [alien]=alien.gtf)

# check files exist/their sizes

for ref in ${NAMES[@]}
do
    echo "$ref"

    echo ${GTFS[$ref]}
    ls -lh ${GTFS[$ref]}

    echo ${FAS[$ref]}
    ls -lh ${FAS[$ref]}

    echo "-----------------"
done

## prepend the species/assembly to each scaffold/chromosome

for ref in ${NAMES[@]}
do
    echo "$ref"
    sed "s/^>/>${ref}_/g" ${FAS[$ref]} > $(basename ${FAS[$ref]}.prepend)

    grep -v "#" ${GTFS[$ref]} | sed "s/^/${ref}_/g" > $(basename ${GTFS[$ref]}.prepend)
done

cat *fa.prepend > combined_correct.fa
cat *gtf.prepend > combined_correct.gtf


STAR --runMode genomeGenerate \
      --runThreadN "$NTHREADS" \
      --genomeDir "$ID" \
      --genomeFastaFiles combined_correct.fa



