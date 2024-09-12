#!/bin/bash
##
## Rerunning v1 experiments with 24 nt genome 
##
## Izaskun Mallona & Giulia Moro
## December 2023

FASTQPATH_1=/home/gmoro/fasqs_third_scRNAseq/
FASTQPATH_2=/home/gmoro/fastqs_fourth_scRNAseq/
FASTQPATH_3=/home/gmoro/fastqs_fith_scRNAseq/

GTF=/home/gmoro/star_solo_v1_experiments/mouse_LTRs_egfp_24_nt/combined_correct.gtf
GENOME_FA=/home/gmoro/star_solo_v1_experiments/mouse_LTRs_egfp_24_nt/combined_correct.fa
IDX=/home/gmoro/star_solo_v1_experiments/mouse_LTRs_egfp_24_nt/mouse_LTRs_egfp_24_nt
SRC=/home/imallona/src/ebrunner_spectral
WHITES="$SRC"/07_barcodes_translation_sbg/data/
NTHREADS=5

WD=/home/gmoro/star_solo_v1_experiments

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

mkdir -p $WD/data_24 ; cd $_

# only taking the files where the genome is with the 24 nt capture

ln -s "$FASTQPATH_1"/20220126.B-o27286_1_1-Unmod1_R1.fastq.gz first_unmod_1_r1.fq.gz
ln -s "$FASTQPATH_1"/20220126.B-o27286_1_1-Unmod1_R2.fastq.gz first_unmod_1_r2.fq.gz

ln -s "$FASTQPATH_1"/20220126.B-o27286_1_2-Mod2_R1.fastq.gz first_mod_2_r1.fq.gz
ln -s "$FASTQPATH_1"/20220126.B-o27286_1_2-Mod2_R2.fastq.gz first_mod_2_r2.fq.gz

ln -s "$FASTQPATH_1"/20220126.B-o27286_1_3-Unmod3_R1.fastq.gz first_unmod_3_r1.fq.gz
ln -s "$FASTQPATH_1"/20220126.B-o27286_1_3-Unmod3_R2.fastq.gz first_unmod_3_r2.fq.gz

ln -s "$FASTQPATH_1"/20220126.B-o27286_1_4-Mod4_R1.fastq.gz first_mod_4_r1.fq.gz
ln -s "$FASTQPATH_1"/20220126.B-o27286_1_4-Mod4_R2.fastq.gz first_mod_4_r2.fq.gz

ln -s "$FASTQPATH_2"/20220421.A-o28127_1_4-SPECTRAL_unmod_R1.fastq.gz second_unmod_r1.fq.gz
ln -s "$FASTQPATH_2"/20220421.A-o28127_1_4-SPECTRAL_unmod_R2.fastq.gz second_unmod_r2.fq.gz

ln -s "$FASTQPATH_2"/20220421.A-o28127_1_1-SPECTRAL_24_R1.fastq.gz second_mod_24_r1.fq.gz
ln -s "$FASTQPATH_2"/20220421.A-o28127_1_1-SPECTRAL_24_R2.fastq.gz second_mod_24_r2.fq.gz

ln -s "$FASTQPATH_3"/20220510.B-o28306_1_1-Unmod_R1.fastq.gz third_unmod_r1.fq.gz
ln -s "$FASTQPATH_3"/20220510.B-o28306_1_1-Unmod_R2.fastq.gz third_unmod_r2.fq.gz

ln -s "$FASTQPATH_3"/20220510.B-o28306_1_2-Mod_R1.fastq.gz third_mod_r1.fq.gz
ln -s "$FASTQPATH_3"/20220510.B-o28306_1_2-Mod_R2.fastq.gz third_mod_r2.fq.gz

mkdir -p "$WD"/star_24 ; cd $_

for item in first_unmod_1 first_mod_2 first_unmod_3 first_mod_4 second_unmod second_mod_24 third_unmod third_mod
do
    echo $item
    r2="$WD"/data_24/"$item"_r2.fq.gz
    r1="$WD"/data_24/"$item"_r1.fq.gz
    
    nice -n 19 STAR --runThreadN "$NTHREADS" \
         --genomeDir "$IDX" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$item"_ / \
         --readFilesIn "$r2" "$r1" \
         --soloType CB_UMI_Complex \
         --soloAdapterSequence ACTGGCCTGCGANNNNNNNNNGGTAGCGGTGACA \
         --soloCBposition 2_-9_2_-1 2_12_2_20 2_34_2_42 \
         --soloUMIposition 3_10_3_17 \
         --soloUMIlen 8 \
         --soloCBwhitelist "$WHITES"/CLS1.txt "$WHITES"/CLS2.txt "$WHITES"/CLS3.txt \
	 --soloCBmatchWLtype 1MM \
	 --soloCellFilter EmptyDrops_CR \
         --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY \
         --soloCellReadStats Standard \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
	 --sjdbOverhang 61 \
         --sjdbGTFfile "$GTF"

	samtools index -@ "$NTHREADS" "$item"_Aligned.sortedByCoord.out.bam
	
	samtools view -H "$item"_Aligned.sortedByCoord.out.bam > header.sam	
	samtools view "$item"_Aligned.sortedByCoord.out.bam | grep "egfp_WPRE" | cat header.sam - | samtools view -Sb > "$item"_subset.bam
	
	rm -rf *.sam
done

rm -rf $WD/star_24/*STARgenome
