#!/bin/bash
##
## Generating .bw files with umi and cb filt 
##
## Izaskun Mallona & Giulia Moro
## December 2023

WD_24=/home/gmoro/star_solo_v1_experiments/star_24
NTHREADS=30

cd "$WD_24"

for item in first_unmod_1 first_mod_2 first_unmod_3 first_mod_4 second_unmod second_mod_24 third_unmod third_mod
do

  echo $item
  bam="$WD_24"/"$item"_subset.bam

  samtools index -@ NTHREADS "$bam"

  TAB="$(printf '\t')" #defining what tab is
  valid_rgs="$WD_24"/"$item"_Solo.out/Gene/filtered/barcodes.tsv #these are the barcodes defined by starsolo, need to change every time

  # extracting unique alignments
  samtools view -H "$bam" > first_header.txt
  samtools view "$bam" | grep -P '^@|NH:i:1\b' | cat first_header.txt - | samtools view -Sb -@ $NTHREADS > unique_filtered.bam

  # extracting barcodes 

  samtools view -D CB:"$valid_rgs" unique_filtered.bam -h -b -o "$item"_filtered.bam
  
  # extracting header

  samtools view -H "$item"_filtered.bam > curr_header.txt
  
  samtools view "$item"_filtered.bam |  sort -k27 -k28 -k 20 -u | cat curr_header.txt - | \
           samtools view -Sb -@ $NTHREADS | samtools sort -@ $NTHREADS > "$item"_cb_ub_filtered.bam
  
  samtools index -@ "$NTRHREADS" "$item"_cb_ub_filtered.bam 

  bamCoverage -b "$item"_cb_ub_filtered.bam -o "$item"_coverage.bw \
  --binSize 1 \
  --numberOfProcessors "$NTHREADS"

  rm curr_header.txt
  rm first_header.txt
  rm unique_filtered.bam
  rm "$item"_filtered.bam
  rm "$item"_cb_ub_filtered.bam
  rm "$item"_cb_ub_filtered.bam.bai

done

WD_44=/home/gmoro/star_solo_v1_experiments/star_44

cd "$WD_44"

for item in second_mod_44
do

  echo $item
  bam="$WD_44"/"$item"_subset.bam

  samtools index -@ NTHREADS "$bam"

  TAB="$(printf '\t')" #defining what tab is
  valid_rgs="$WD_44"/"$item"_Solo.out/Gene/filtered/barcodes.tsv #these are the barcodes defined by starsolo, need to change every time

  # extracting unique alignments

  samtools view -H "$bam" > first_header.txt
  samtools view "$bam" | grep -P '^@|NH:i:1\b' | cat first_header.txt - | samtools view -Sb -@ $NTHREADS > unique_filtered.bam

  # extracting barcodes 

  samtools view -D CB:"$valid_rgs" unique_filtered.bam -h -b -o "$item"_filtered.bam  

  # extracting header

  samtools view -H "$item"_filtered.bam > curr_header.txt
  
  samtools view "$item"_filtered.bam |  sort -k27 -k28 -k 20 -u | cat curr_header.txt - | \
           samtools view -Sb -@ $NTHREADS | samtools sort -@ $NTHREADS > "$item"_cb_ub_filtered.bam
  
  samtools index -@ "$NTRHREADS" "$item"_cb_ub_filtered.bam 

  bamCoverage -b "$item"_cb_ub_filtered.bam -o "$item"_coverage.bw \
  --binSize 1 \
  --numberOfProcessors "$NTHREADS"

  rm curr_header.txt
  rm unique_filtered.bam
  rm first_header.txt
  rm "$item"_filtered.bam
  rm "$item"_cb_ub_filtered.bam
  rm "$item"_cb_ub_filtered.bam.bai

done

WD_3_prime=/home/gmoro/star_solo_v1_experiments/star_3_prime

cd "$WD_3_prime"

for item in second_mod_3_prime
do

  echo $item
  bam="$WD_3_prime"/"$item"_subset.bam

  samtools index -@ NTHREADS "$bam"

  TAB="$(printf '\t')" #defining what tab is
  valid_rgs="$WD_3_prime"/"$item"_Solo.out/Gene/filtered/barcodes.tsv #these are the barcodes defined by starsolo, need to change every time

  # extracting unique alignments

  samtools view -H "$bam" > first_header.txt
  samtools view "$bam" | grep -P '^@|NH:i:1\b' | cat first_header.txt - | samtools view -Sb -@ $NTHREADS > unique_filtered.bam

  # extracting barcodes 

  samtools view -D CB:"$valid_rgs" unique_filtered.bam -h -b -o "$item"_filtered.bam  

  # extracting header

  samtools view -H "$item"_filtered.bam > curr_header.txt
  
  samtools view "$item"_filtered.bam |  sort -k27 -k28 -k 20 -u | cat curr_header.txt - | \
           samtools view -Sb -@ $NTHREADS | samtools sort -@ $NTHREADS > "$item"_cb_ub_filtered.bam
  
  samtools index -@ "$NTRHREADS" "$item"_cb_ub_filtered.bam 

  bamCoverage -b "$item"_cb_ub_filtered.bam -o "$item"_coverage.bw \
  --binSize 1 \
  --numberOfProcessors "$NTHREADS"

  rm curr_header.txt
  rm first_header.txt
  rm unique_filtered.bam
  rm "$item"_filtered.bam
  rm "$item"_cb_ub_filtered.bam
  rm "$item"_cb_ub_filtered.bam.bai

done

