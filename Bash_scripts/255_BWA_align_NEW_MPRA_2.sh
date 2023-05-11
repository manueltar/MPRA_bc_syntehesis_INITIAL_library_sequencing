#!/bin/bash


#### PRIOR

# ######################### FASTQ FILE DOWNLOAD

# /software/team151/bin/bs download project -i 143809676 -o /lustre/scratch115/teams/soranzo/projects/NEW_MPRA/236/
# /software/team151/bin/bs download project -i 144266130 -o /lustre/scratch115/teams/soranzo/projects/NEW_MPRA/240/


# ######################### Get the quality


# NEW_MPRA$ sh ~/Scripts/Wraper_scripts/254_Illumina_MiSeq_FASTQ_qualities.sh



path=$1
mem=$2
pc=$3

tag=$7

REFERENCE=$4
fq_r1=$5
fq_r2=$6

cd $path

bsub -G team151 -o $path/Prep.out -M 4000 -J Quality_fq_r1_path1 -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"fastqc $fq_r1 --outdir $path"

bsub -G team151 -o $path/Prep.out -M 4000 -J Quality_fq_r2_path1 -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"fastqc $fq_r2 --outdir $path"


## first trimm 20 bp at the end of the R1 and the R2 reads


output_trimmomatic_r1_paired=$(echo "Trimmed_paired_""$fq_r1")
output_trimmomatic_r1_unpaired=$(echo "Trimmed_unpaired_""$fq_r1")

output_trimmomatic_r2_paired=$(echo "Trimmed_paired_""$fq_r2")
output_trimmomatic_r2_unpaired=$(echo "Trimmed_unpaired_""$fq_r2")


bsub -G team151 -o  $path/Prep.out -M 4000 -J Trimmomatic -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n 1 -q normal -- \
"java -jar /nfs/users/nfs_m/mt19/sOFTWARE/Trimmomatic-0.39/trimmomatic-0.39.jar PE $fq_r1 $fq_r2 $output_trimmomatic_r1_paired $output_trimmomatic_r1_unpaired \
 $output_trimmomatic_r2_paired $output_trimmomatic_r2_unpaired CROP:115"


 
## Align in stringent conditions 



#r1_to_align=$(echo "$output_trimmomatic_r1_paired" | sed -r 's/\.gz//g')
#r2_to_align=$(echo "$output_trimmomatic_r2_paired" | sed -r 's/\.gz//g')

r1_to_align=$(echo "$output_trimmomatic_r1_paired")
r2_to_align=$(echo "$output_trimmomatic_r2_paired")



bsub -G team151 -o $path/Prep.out -M 4000 -w"done(Trimmomatic)" -J Quality2_fq_r1_path1 -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"fastqc $r1_to_align --outdir $path"

bsub -G team151 -o $path/Prep.out -M 4000 -w"done(Trimmomatic)"  -J Quality2_fq_r2_path1 -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"fastqc $r2_to_align --outdir $path"


#### FLASH to merger reads



bsub -G team151 -o $path/FLASH.out -M 4000 -w"done(Trimmomatic)" -J FLASH -R"select[mem>=4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"flash -z --output-prefix=$tag -r 115 -f 185 -s 20 -m 10 -x 0.122 $r1_to_align $r2_to_align"

#exit

#####

flash_output=$(echo "$tag"".extendedFrags.fastq.gz")
OUTPUT_BWA=$path/$tag.sam
read_group_id=$(echo "'@RG"'\t'"ID:1"'\t'"SM:$tag"'\t'"PL:Illumina"'\t'"LB:1"'\t'"PU:1'")

bsub -G team151 -o  $path/BWA_alingment.out -M $mem -w"done(FLASH)" -J BWA -R"select[mem >$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q normal -- \
"bwa mem -M -R $read_group_id  -t 2 $REFERENCE $flash_output > $OUTPUT_BWA"

######

OUTPUT_BWA_StoB=$path/$tag.bam

#bsub -G team151 -o  $path/BWA_alingment.out -M 4000 -J StoB -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
bsub -G team151 -o  $path/BWA_alingment.out -M $mem -w"done(BWA)" -J StoB -R"select[mem >$mem] rusage[mem=$mem] span[hosts=1]" -n$pc -q normal -- \
"/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools view -b $OUTPUT_BWA > $OUTPUT_BWA_StoB"

bsub -G team151 -o  $path/BWA_alingment.out -M 4000 -w"done(StoB)" -J clear1 -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm $OUTPUT_BWA"



OUTPUT_BWA_StoB_sorted=$(echo "$path""/""$tag""_sorted.bam")

bsub -G team151 -o  $path/BWA_alingment.out -M 4000 -w"done(clear1)" -J SORT -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools sort -o $OUTPUT_BWA_StoB_sorted $OUTPUT_BWA_StoB"

bsub -G team151 -o  $path/BWA_alingment.out -M 4000 -w"done(SORT)" -J INDEX -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"/nfs/users/nfs_m/mt19/sOFTWARE/samtools-1.6/bin/samtools index $OUTPUT_BWA_StoB_sorted"

bsub -G team151 -o  $path/BWA_alingment.out -M 4000 -w"done(INDEX)" -J clear2 -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal -- \
"rm $OUTPUT_BWA_StoB"
