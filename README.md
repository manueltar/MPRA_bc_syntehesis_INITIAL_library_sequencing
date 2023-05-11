
########################### INITIAL_library

# The PE150 reads are in the folder: 

/nfs/team151_data03/MPRA_Programmed/INITIAL_library_sequencing/ former 240

# The manifest file is: Manifest_19_11_2019.csv

# The global file reads are the result of cats:

cat INITIAL_library_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a/INITIAL-library-PrepShendure_S2_L001_R1_001.fastq.gz INITIAL_library_PrepShendure_L001_ds.72543c18e4f44ba3819030ac539d5d78/INITIAL-library-PrepShendure_S2_L001_R1_001.fastq.gz INITIAL_library_2STEP_Prep_L001_ds.8584a561088340828a4c86\
2eabcbb725/INITIAL-library-2STEP-Prep_S1_L001_R1_001.fastq.gz INITIAL_library_2STEP_Prep_L001_ds.9a3e73d96cc8460f8bec4d281833fed7/INITIAL-library-2STEP-Prep_S1_L001_R1_001.fastq.gz > global_R1.fastq.gz


cat INITIAL_library_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a/INITIAL-library-PrepShendure_S2_L001_R2_001.fastq.gz INITIAL_library_PrepShendure_L001_ds.72543c18e4f44ba3819030ac539d5d78/INITIAL-library-PrepShendure_S2_L001_R2_001.fastq.gz INITIAL_library_2STEP_Prep_L001_ds.8584a561088340828a4c86\
 2eabcbb725/INITIAL-library-2STEP-Prep_S1_L001_R2_001.fastq.gz INITIAL_library_2STEP_Prep_L001_ds.9a3e73d96cc8460f8bec4d281833fed7/INITIAL-library-2STEP-Prep_S1_L001_R2_001.fastq.gz > global_R2.fastq.gz



# The reference for the PE150 alignment is: /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/reference_files/Library_TRIMMED_15bp_RMV_M.fasta


ALIGN with trimming and BWA

25$ sh ~/Scripts/Wraper_scripts/255_BWA_align_NEW_MPRA.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/240/INITIAL_library_2STEP_Prep_L001_ds.8584a561088340828a4c862eabcbb725/ 16000 4 INITIAL-library-2STEP-Prep_S1_L001_R1_001.fastq.gz  INITIAL-library-2STEP-Prep_S1_L001_R2_001.fastq.gz INITIAL_library_2ST\
EP_Prep_L001_ds.8584a561088340828a4c862eabcbb725


ed7$ sh ~/Scripts/Wraper_scripts/255_BWA_align_NEW_MPRA.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/240/INITIAL_library_2STEP_Prep_L001_ds.9a3e73d96cc8460f8bec4d281833fed7/ 16000 4 INITIAL-library-2STEP-Prep_S1_L001_R1_001.fastq.gz INITIAL-library-2STEP-Prep_S1_L001_R2_001.fastq.gz INITIAL_library_2ST\
EP_Prep_L001_ds.9a3e73d96cc8460f8bec4d281833fed7



a$ sh ~/Scripts/Wraper_scripts/255_BWA_align_NEW_MPRA.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/240/INITIAL_library_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a/  16000 4 INITIAL-library-PrepShendure_S2_L001_R1_001.fastq.gz INITIAL-library-PrepShendure_S2_L001_R2_001.fastq.gz INITIAL_librar\
y_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a


78$ sh ~/Scripts/Wraper_scripts/255_BWA_align_NEW_MPRA.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/240/INITIAL_library_PrepShendure_L001_ds.72543c18e4f44ba3819030ac539d5d78/ 16000 4 INITIAL-library-PrepShendure_S2_L001_R1_001.fastq.gz INITIAL-library-PrepShendure_S2_L001_R2_001.fastq.gz INITIAL_librar\
y_PrepShendure_L001_ds.72543c18e4f44ba3819030ac539d5d78




#### Pysam to calculate summaries of the PE150 alignments

240$ ./process_mpra_MTS_edited.py --b INITIAL_library_2STEP_Prep_L001_ds.8584a561088340828a4c862eabcbb725/INITIAL_library_2STEP_Prep_L001_ds.8584a561088340828a4c862eabcbb725_sorted.bam --r ../reference_files/Library_TRIMMED_15bp_RMV_M.fasta --o Eugene_2STEP_Prep_L001_ds_8.txt

240$ ./process_mpra_MTS_edited.py --b INITIAL_library_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a/INITIAL_library_PrepShendure_L001_ds.ef7bc87c2fd84cd0922f09f950a2996a_sorted.bam --r ../reference_files/Library_TRIMMED_15bp_RMV_M.fasta --o Eugene_PrepShendure_L001_ds_e.txt



#### R scripts to characterise the capture in two different methods 2_STEP and PrepShendure

240$ bsub -G team151 -o 2_STEP.out -M 4000 -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q long "/software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/124_MPRA_PE150_Rscript_edition.R --INPUT_TABLE ~/RareVar_2019/DEF_design/Resynthesis/Characterization_PE150_INITIAL_TABLE.txt --Eugene_summar\
y_1 Eugene_2STEP_Prep_L001_ds_8.txt  --Eugene_summary_2 Eugene_2STEP_Prep_L001_ds_9.txt --type 2_STEP --out 2_STEP"

240$ bsub -G team151 -o 2_STEP_PartII.out -M 4000 -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal "/software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/125_MPRA_PE150_PartII_Rscript_edition.R --INPUT_TABLE ~/RareVar_2019/DEF_design/Resynthesis/Characterization_PE150_INITIAL_TABLE.txt\
 --CAPTURED Captured_2_STEP.txt --DROPPED_OUT Dropped_out_2_STEP.txt --AFFECTED_BC MM_in_BC_2_STEP.txt --AFFECTED_ELSEWHERE MM_elsewhere_2_STEP.txt --type 2_STEP --out 2_STEP"

240$/software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/128_MPRA_PE150_2STEP_Part_III_Rscript_edition.R --INPUT_TABLE ~/RareVar_2019/DEF_design/Resynthesis/Characterization_PE150_INITIAL_TABLE.txt --Shendure DEF_MUT_BC_2_STEP.txt --type 2_STEP --out 2_STEP




240$ bsub -G team151 -o Shendure.out -M 12000 -R"select[mem >12000] rusage[mem=12000] span[hosts=1]" -n3 -q normal "/software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/124_MPRA_PE150_Rscript_edition.R --INPUT_TABLE ~/RareVar_2019/DEF_design/Resynthesis/Characterization_PE150_INITIAL_TABLE.txt --Eugene\
_summary_1 Eugene_PrepShendure_L001_ds_e.txt  --Eugene_summary_2 Eugene_PrepShendure_L001_ds_7.txt --type PrepShendure --out PrepShendure"

240$ bsub -G team151 -o PrepShendure_PartII.out -M 4000 -R"select[mem >4000] rusage[mem=4000] span[hosts=1]" -n1 -q normal "/software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/125_MPRA_PE150_PartII_Rscript_edition.R --INPUT_TABLE /nfs/users/nfs_m/mt19/RareVar_2019/DEF_design/Resynthesis/Characterizat\
ion_PE150_INITIAL_TABLE.txt --CAPTURED Captured_PrepShendure.txt --DROPPED_OUT Dropped_out_PrepShendure.txt --AFFECTED_BC MM_in_BC_PrepShendure.txt --AFFECTED_ELSEWHERE MM_elsewhere_PrepShendure.txt --type PrepShendure --out PrepShendure"

240$ /software/R-3.5.1/bin/Rscript /nfs/users/nfs_m/mt19/Scripts/R/128_MPRA_PE150_2STEP_Part_III_Rscript_edition.R --INPUT_TABLE ~/RareVar_2019/DEF_design/Resynthesis/Characterization_PE150_INITIAL_TABLE.txt --Shendure DEF_MUT_BC_PrepShendure.txt --type PrepShendure --out PrepShendure



# This is the alignment line for the global version after fusing the reads

bash ~/Scripts/Wraper_scripts/255_BWA_align_NEW_MPRA_2.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/240/ 16000 4 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/reference_files/Library_TRIMMED_15bp_RMV_M.fasta  global_R1.fastq.gz global_R2.fastq.gz global


########################### First assay of the library FINAL in the plasmids

# The PE 25 reads are in /nfs/team151_data03/MPRA_Programmed/FINAL_library_sequencing/ former 236

# The manifest file is Manifest_5_11_2019.csv

NEW_MPRA/246$ scp -r /lustre/scratch117/sciops/team117/npg/srl/mt9/191210_MS2_MiSeq_walk-up_246_A_MS8539685-050V2/* ./

236$ sh ~/Scripts/Wraper_scripts/260_extractor_of_counts.sh /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/236/ 32000 8 /lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/reference_files/Library_bc_RMV_M.fasta 1_S1_L001_R1_001.fastq.gz 1_S1_L001_R2_001.fastq.gz FINAL plasmid

## Count the abundancies in the final library to see if they differ from the INITIAL library
