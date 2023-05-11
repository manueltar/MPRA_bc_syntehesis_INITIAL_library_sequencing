

suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("Biostrings", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("splitstackshape", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("data.table", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))

opt = NULL

create_intermediate_files_1 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  #### READ and transform INPUT_TABLE LIABILITY!!!!! ----
  
  INPUT_table = read.table(opt$INPUT_TABLE, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_0_\n")
  cat(sprintf(as.character(INPUT_table$seq_name_bc[1])))
  cat("\n")
  
  #### INPUT check
  
  check<-paste(INPUT_table$KEY,
                           INPUT_table$Tile,
                           INPUT_table$Type,
                           INPUT_table$barcode,
                           sep=";")
  cat("Manuel_0.5_\n")
  cat(sprintf(as.character(check[1])))
  cat("\n")
  
  #### READ CAPTURED ----
  
  CAPTURED = read.table(opt$CAPTURED, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.5_\n")
  cat(sprintf(as.character(CAPTURED[1,])))
  cat("\n")
  
  # cat("Manuel_1.5_\n")
  # cat(sprintf(as.character(CAPTURED$summary$element[1])))
  # cat("\n")
  
  
  #### READ DROPPED ----
  
  DROPPED_OUT = read.table(opt$DROPPED_OUT, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.75_\n")
  cat(sprintf(as.character(DROPPED_OUT[1,])))
  cat("\n")
  
  # cat("Manuel_1.75_\n")
  # cat(sprintf(as.character(Eugene_summary_2$summary$element[1])))
  # cat("\n")
  
  #### READ AFFECTED_BC ----
  
  Affected_BC = read.table(opt$AFFECTED_BC, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.5_\n")
  cat(sprintf(as.character(Affected_BC[1,])))
  cat("\n")
  
  # cat("Manuel_1.5_\n")
  # cat(sprintf(as.character(CAPTURED$summary$element[1])))
  # cat("\n")
  
  #### READ AFFECTED_ELSEWHERE ----
  
  Affected_elsewhere = read.table(opt$AFFECTED_ELSEWHERE, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.5_\n")
  cat(sprintf(as.character(Affected_elsewhere[1,])))
  cat("\n")
  
  # cat("Manuel_1.5_\n")
  # cat(sprintf(as.character(AFFECTED_ELSEWHERE$summary$element[1])))
  # cat("\n")
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("Manuel_3_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("Manuel_4_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### collapse elsewhere mutations to mutations per seq-bc ----
    
  Affected_elsewhere.dt<-data.table(Affected_elsewhere, key="seq_name_bc")
  
  Affected_elsewhere_COLLAPSED<-as.data.frame(Affected_elsewhere.dt[, 
                                                                    .(Mut_alleles=.N, 
                                                                      Mut_string=paste(mutations_STRING, collapse="|"),
                                                                      Mut_string_reads=paste(reads_Element_mutated, collapse="|"),
                                                                      Mut_total_reads=sum(reads_Element_mutated),
                                                                      WT_reads=max(reads_Element_perfect)), 
                                                                    by = seq_name_bc])
  
  Affected_elsewhere_COLLAPSED$TOTAL_READS<-Affected_elsewhere_COLLAPSED$Mut_total_reads + Affected_elsewhere_COLLAPSED$WT_reads
  
  Affected_elsewhere_COLLAPSED$Perc_Mutated<-100*(Affected_elsewhere_COLLAPSED$Mut_total_reads/Affected_elsewhere_COLLAPSED$TOTAL_READS)
  
  
  
  SUMMARY_table_1<-setDT(Affected_elsewhere_COLLAPSED)[, .N, by=Perc_Mutated]
  
  
  SUMMARY_table_1$Factor3<-cut(SUMMARY_table_1$Perc_Mutated, breaks = c(-Inf,1,20,50,75,99.99,Inf))
  
  SUMMARY_table_2<-setDT(SUMMARY_table_1)[, sum(N), by=Factor3]
  
  
  colnames(SUMMARY_table_2)<-c("Perc_mutated_reads","seq-bc")
  
  #### loop mutations in the SNP ----
  
  Gather_MUT_AFFECTS_REL_POS<-data.frame(matrix(vector(), 0, 
                                                4,
                                                dimnames=list(c(),c("seq_name_bc",
                                                                    "Mut_string","Mutation_affecting",
                                                                    "Mut_string_reads")
                                                )),
                                         stringsAsFactors=F)
  
  for(i in 1:length(INPUT_table$seq_name_bc))
  {
    seq_name_bc_plus_index_SEL<-INPUT_table$seq_name_bc[i]
    
    INPUT_table_SEL<-INPUT_table[i,]
    
    seq_name_bc_SEL<-paste(INPUT_table_SEL$KEY,
                           INPUT_table_SEL$Tile,
                           INPUT_table_SEL$Type,
                           INPUT_table_SEL$barcode,
                           sep=";")
    
    REL_POS_SEL<-as.numeric(INPUT_table_SEL$RELPOS_TRIMMED)
    
    ### Mutations affect key SNP POS?
    
    Affected_elsewhere_COLLAPSED_SEL<-Affected_elsewhere_COLLAPSED[which(Affected_elsewhere_COLLAPSED$seq_name_bc == seq_name_bc_SEL),]
    
    if(length(Affected_elsewhere_COLLAPSED_SEL$seq_name_bc) > 0 & !is.na(REL_POS_SEL))
    {
      TOTAL_READS_SEL<-as.numeric(Affected_elsewhere_COLLAPSED_SEL$TOTAL_READS)
      
      mut_string_COMPOSITE<-as.character(Affected_elsewhere_COLLAPSED_SEL$Mut_string)
      
      read_string_COMPOSITE<-as.character(Affected_elsewhere_COLLAPSED_SEL$Mut_string_reads)
      
      
      mut_string_COMPOSITE_array<-unlist(strsplit(mut_string_COMPOSITE, "\\|"))
      
      read_string_COMPOSITE<-unlist(strsplit(read_string_COMPOSITE, "\\|"))
      
      
      for(j in 1:length(mut_string_COMPOSITE_array))
      {
        mut_string_SEL<-mut_string_COMPOSITE_array[j]
        read_string_SEL<-read_string_COMPOSITE[j]
        
        mutations_array<-unlist(strsplit(mut_string_SEL, "\\;"))
        #read_array<-unlist(strsplit(read_string_SEL, "\\;"))
        read_SEL<-as.numeric(read_string_SEL)
        
        for(k in 1:length(mutations_array))
        {
          mutation_SEL<-mutations_array[k]
          
          
          if(sum(grep("-",mutation_SEL))>0)
          {
            INTERVAL_MUT<-as.character(gsub("_\\w+\\:\\w+","",mutation_SEL, perl=T))
            
            INTERVAL_MUT_START<-as.numeric(gsub("-\\w+","",INTERVAL_MUT, perl=T))
            INTERVAL_MUT_END<-as.numeric(gsub("\\w+-","",INTERVAL_MUT, perl=T))
            
            INTERVAL_POS<-seq(INTERVAL_MUT_START,INTERVAL_MUT_END,by=1)
            
            for(l in 1:length(INTERVAL_POS))
            {
              POS_MUT<-INTERVAL_POS[l]
              if(POS_MUT == REL_POS_SEL)
              {
                
                A_muts<-as.data.frame(cbind(seq_name_bc_SEL,mut_string_SEL,
                                            mutation_SEL,read_SEL))
                
                colnames(A_muts)<-colnames(Gather_MUT_AFFECTS_REL_POS)
                
                Gather_MUT_AFFECTS_REL_POS<-rbind(Gather_MUT_AFFECTS_REL_POS,A_muts)
                
              }else{
                
                # Do nothing
              }
              
            }
            
            
          } #INTERVAL DEL
          else{
            if(sum(grep(">",mutation_SEL))>0)
            {
              POS_MUT<-as.numeric(gsub("_\\w+\\>\\w+","",mutation_SEL, perl=T))
              
              if(is.na(POS_MUT))
              {
                
                POS_MUT<-as.numeric(gsub("_\\w+\\>","",mutation_SEL, perl=T))
              }
              
              if(POS_MUT == REL_POS_SEL)
              {
                A_muts<-as.data.frame(cbind(seq_name_bc_SEL,mut_string_SEL,
                                            mutation_SEL,read_SEL))
                
                colnames(A_muts)<-colnames(Gather_MUT_AFFECTS_REL_POS)
                
                Gather_MUT_AFFECTS_REL_POS<-rbind(Gather_MUT_AFFECTS_REL_POS,A_muts)
                
              }else{
                
                # Do nothing
              }
              
              
            }# Mismatch
            else{
              
              if(sum(grep(":",mutation_SEL))>0)
              {
                POS_MUT<-as.numeric(gsub("_\\w+\\:\\w+","",mutation_SEL, perl=T))
                if(POS_MUT == REL_POS_SEL)
                {
                  A_muts<-as.data.frame(cbind(seq_name_bc_SEL,mut_string_SEL,
                                              mutation_SEL,read_SEL))
                  
                  colnames(A_muts)<-colnames(Gather_MUT_AFFECTS_REL_POS)
                  
                  Gather_MUT_AFFECTS_REL_POS<-rbind(Gather_MUT_AFFECTS_REL_POS,A_muts)
                  
                  
                }else{
                  
                  # Do nothing
                }
              } 
              else{
                
                if(sum(grep("NONE",mutation_SEL))>0)
                {
                  # Do nothing
                  
                } # NONE mutation
                else{
                  
                  # ANY INSTANCE NOT CONTEMPLATED?
                }
              }# DEL or insertion no INTERVAL
            } # Mismatch
          }# #INTERVAL DEL
          
        }# k for individual mutations
        
      }# j mutation_strings in COMPOSITE
    }# if the variant is in the Affected_elsewhere_COLLAPSED set
  }#i
  
  #### COLLAPSE mutations in the SNP ----
  
  Gather_MUT_AFFECTS_REL_POS$Mut_string_reads<-as.numeric(as.character(Gather_MUT_AFFECTS_REL_POS$Mut_string_reads))
  
  Gather_MUT_AFFECTS_REL_POS.dt<-data.table(Gather_MUT_AFFECTS_REL_POS, key="seq_name_bc")
  
  
  COLLAPSE.SNP<-Gather_MUT_AFFECTS_REL_POS.dt[, .(mut_string_COMPOSITE=paste(Mut_string, collapse="|"),
                                                  mut_affecting_COMPOSITE=paste(Mutation_affecting, collapse="|"),
                                                  SNP_affected_reads=sum(Mut_string_reads)), by=seq_name_bc]
  
  #### loop add mutations in the SNP and mutations elsewhere to the INPUT_TABLE ----
  
  
  DEF_Mut<-data.frame(matrix(vector(), 0, 
                             23,
                             dimnames=list(c(),c(colnames(INPUT_table),
                                                 "seq_name_bc_NO_INDEX","Mut_alleles","Mut_string","Mut_string_reads",
                                                 "Mut_total_reads","WT_reads","TOTAL_READS","Perc_Mutated",
                                                 "SNP.affected_string","Perc.SNP.affected")
                             )),
                      stringsAsFactors=F)
  
  
  for(i in 1:length(INPUT_table$seq_name_bc))
  {
    seq_name_bc_plus_index_SEL<-INPUT_table$seq_name_bc[i]
    
    INPUT_table_SEL<-INPUT_table[i,]
    
    seq_name_bc_SEL<-paste(INPUT_table_SEL$KEY,
                           INPUT_table_SEL$Tile,
                           INPUT_table_SEL$Type,
                           INPUT_table_SEL$barcode,
                           sep=";")
    
    REL_POS_SEL<-as.numeric(INPUT_table_SEL$RELPOS_TRIMMED)
    
    ### Mutations affect key SNP POS?
    
    Affected_elsewhere_COLLAPSED_SEL<-Affected_elsewhere_COLLAPSED[which(Affected_elsewhere_COLLAPSED$seq_name_bc == seq_name_bc_SEL),]
    
    if(length(Affected_elsewhere_COLLAPSED_SEL$seq_name_bc) > 0)
    {
      COLLAPSE.SNP.SEL<-COLLAPSE.SNP[which(COLLAPSE.SNP$seq_name_bc == seq_name_bc_SEL),]
      TOTAL_READS_SEL<-Affected_elsewhere_COLLAPSED_SEL$TOTAL_READS
      Perc.SNP.affected<-0
      SNP.affected_string<-"NONE"
      
      if(length(COLLAPSE.SNP.SEL$seq_name_bc) > 0)
      {
        Perc.SNP.affected<-100*(COLLAPSE.SNP.SEL$SNP_affected_reads/TOTAL_READS_SEL)
        SNP.affected_string<-COLLAPSE.SNP.SEL$mut_affecting_COMPOSITE
        
      }else{
        
        # Do nothing
      }
      
      Affected_elsewhere_COLLAPSED_SEL$SNP.affected_string<-SNP.affected_string
      Affected_elsewhere_COLLAPSED_SEL$Perc.SNP.affected<-Perc.SNP.affected
      
      B_muts<-as.data.frame(cbind(INPUT_table_SEL,Affected_elsewhere_COLLAPSED_SEL))
      
      colnames(B_muts)<-colnames(DEF_Mut)
      DEF_Mut<-rbind(DEF_Mut,B_muts)
      
    }# if exist Affected_elsewhere_COLLAPSED_SEL
    else{
      
      # Dropout
      # NA.vector<-as.data.frame(rep("NA",10))
      # NA.vector.t<-t(NA.vector)
      # 
      # B_muts<-as.data.frame(cbind(INPUT_table_SEL,NA.vector.t))
      # 
      # colnames(B_muts)<-colnames(DEF_Mut)
      # DEF_Mut<-rbind(DEF_Mut,B_muts)
    }
  } #i
  
  
  #### Add BC affected reads ----
  
  colnames(Affected_BC)[which(colnames(Affected_BC) == "seq_name_bc_sel")]<-"seq_name_bc_NO_INDEX"
  colnames(Affected_BC)[which(colnames(Affected_BC) == "affected_reads")]<-"BC_affected_reads"
  
  indexes_subset<-c(which(colnames(Affected_BC) == "seq_name_bc_NO_INDEX"),which(colnames(Affected_BC) == "BC_affected_reads"))
  
  Affected_BC_subset<-Affected_BC[,indexes_subset]
  
  DEF_Mut_BC<-merge(DEF_Mut,
                    Affected_BC_subset,
                    by="seq_name_bc_NO_INDEX",
                    ALL=T)
  
  DEF_Mut_BC$Perc.BC.affected<-100*(DEF_Mut_BC$BC_affected_reads/DEF_Mut_BC$TOTAL_READS)
  
  #### SAVE captured and dropout ----
  
  filename_1<-paste("DEF_MUT_BC_",type,".txt", sep='')

  write.table(DEF_Mut_BC,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}

#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--INPUT_TABLE"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CAPTURED"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DROPPED_OUT"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AFFECTED_BC"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AFFECTED_ELSEWHERE"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "125_MPRA_PE150_PartII_Rscript_edition.R
                        --INPUT_TABLE FILE.txt
                        --CAPTURED FILE.txt 
                        --DROPPED_OUT FILE.txt 
                        --AFFECTED_BC FILE.txt
                        --AFFECTED_ELSEWHERE FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  create_intermediate_files_1(opt)
  
}


###########################################################################

system.time( main() )