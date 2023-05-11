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
  
  INPUT_TABLE = read.table(opt$INPUT_TABLE, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_0_\n")
  cat(sprintf(as.character(INPUT_TABLE[1,])))
  cat(sprintf(as.character(INPUT_TABLE$seq_name_bc[1])))
  cat("\n")
  
  #### INPUT check
  
  INPUT_TABLE$check<-paste(INPUT_TABLE$KEY,
                           INPUT_TABLE$Tile,
                           INPUT_TABLE$Type,
                           INPUT_TABLE$barcode,
                           sep=";")
  cat("Manuel_0.5_\n")
  cat(sprintf(as.character(INPUT_TABLE$check[1])))
  cat("\n")
  
  #### READ and transform Eugene_summary_1 ----
  
  Eugene_summary_1 = read.table(opt$Eugene_summary_1, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.5_\n")
  cat(sprintf(as.character(Eugene_summary_1[1,])))
  cat("\n")
  
  # cat("Manuel_1.5_\n")
  # cat(sprintf(as.character(Eugene_summary_1$summary$element[1])))
  # cat("\n")
  
  
  #### READ and transform Eugene_summary_2 ----
  
  Eugene_summary_2 = read.table(opt$Eugene_summary_2, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.75_\n")
  cat(sprintf(as.character(Eugene_summary_2[1,])))
  cat("\n")
  
  # cat("Manuel_1.75_\n")
  # cat(sprintf(as.character(Eugene_summary_2$summary$element[1])))
  # cat("\n")
  
  #### Merge both summary files ----
  
  summary<-rbind(Eugene_summary_1,Eugene_summary_2)
  
  cat("Manuel_1.95_\n")
  cat(sprintf(as.character(summary$element[1])))
  cat("\n")
  
  #### Summary file obtain seq_name_bc
  
  summary$seq_name_bc<-paste(summary$element, summary$tag, sep=":")
  
  summary$seq_name_bc<-gsub(":",";",summary$seq_name_bc)
  
  cat("Manuel_2_\n")
  cat(sprintf(as.character(summary$seq_name_bc[1])))
  cat("\n")
  
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
  
  #### Table 1 TOTAL reads found
  
  cat("Manuel_5.0_\n")
  cat(sprintf(as.character(summary$seq_name_bc[1])))
  str(summary)
  cat("\n")
  
  TABLE_1<-as.data.frame(table(summary$seq_name_bc))
  
  colnames(TABLE_1)<-c("seq_name_bc","reads")
  
  TABLE_1$seq_name_bc<-as.character(TABLE_1$seq_name_bc)
  
  
  cat("Manuel_5.1_\n")
  str(TABLE_1)
  cat(sprintf(as.character(TABLE_1[1,])))
  cat("\n")
  
  #### captured seq-bc
  
  Captured<-INPUT_TABLE[which(INPUT_TABLE$check%in%TABLE_1$seq_name_bc),]
  
  #### Dropout seq-bc
  
  DO<-INPUT_TABLE[-which(INPUT_TABLE$check%in%TABLE_1$seq_name_bc),]
  
  #### SAVE captured and dropout ----

  filename_1<-paste("Captured_",type,".txt", sep='')

  write.table(Captured, 
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  filename_2<-paste("Dropped_out_",type,".txt", sep='')
  
  write.table(DO, 
              file=filename_2, sep="\t", quote=F, row.names = F)
  
  #### Table 2 reads with MM in the bc ----
  
  summary(as.factor(summary$tag_matches))
  
  summary.dt<-data.table(summary, key="seq_name_bc")
  
  Interaction_table_2<-as.data.frame(summary.dt[, .N ,by = .(tag_matches,seq_name_bc)])
  
  colnames(Interaction_table_2)[which(colnames(Interaction_table_2) == "tag_matches")]<-"Perfect_bc"
  colnames(Interaction_table_2)[which(colnames(Interaction_table_2) == "N")]<-"Freq"
  
  #### loop
  
  # cat("Manuel_5.75_\n")
  # str(Interaction_table_2)
  # cat("\n")
  
  Gather_perfect<-data.frame(matrix(vector(), 0, 2,
                                    dimnames=list(c(),c("seq_name_bc",
                                                        "reads_BC_perfect")
                                    )),
                             stringsAsFactors=F)
  
  Gather_affected<-data.frame(matrix(vector(), 0, 2,
                                     dimnames=list(c(),c("seq_name_bc",
                                                         "reads_BC_affected")
                                     )),
                              stringsAsFactors=F)
  
  for(i in 1:length(Interaction_table_2$seq_name_bc))
  {
    seq_name_bc_sel<-as.character(Interaction_table_2$seq_name_bc[i])
    Tag_sel<-as.character(Interaction_table_2$Perfect_bc[i])
    
    #str(Tag_sel)
    
    if(Tag_sel == "True")
    {
      perfect_reads<-as.numeric(Interaction_table_2$Freq[i])
      
      A<-as.data.frame(cbind(seq_name_bc_sel,perfect_reads))
      
      Gather_perfect<-rbind(Gather_perfect,A)
      
    }else{
      
      affected_reads<-as.numeric(Interaction_table_2$Freq[i])
      
      B<-as.data.frame(cbind(seq_name_bc_sel,affected_reads))
      
      Gather_affected<-rbind(Gather_affected,B)
      
    }
    
  }
  
  # cat("Manuel_6_\n")
  # str(Gather_affected)
  # cat("\n")
  # 
  # cat("Manuel_7_\n")
  # str(Gather_perfect)
  # cat("\n")
  
  DEF.Table2<-merge(Gather_affected,
                    Gather_perfect,
                    by="seq_name_bc_sel",
                    all=T)
  
  DEF.Table2$affected_reads<-as.numeric(as.character(DEF.Table2$affected_reads))
  DEF.Table2$affected_reads[which(is.na(DEF.Table2$affected_reads))]<-0
  
  DEF.Table2$perfect_reads<-as.numeric(as.character(DEF.Table2$perfect_reads))
  DEF.Table2$perfect_reads[which(is.na(DEF.Table2$perfect_reads))]<-0
  
  
  DEF.Table2$TOTAL<-DEF.Table2$affected_reads+DEF.Table2$perfect_reads
  
  DEF.Table2$Perc_affected<-100*(DEF.Table2$affected_reads/DEF.Table2$TOTAL)
  
  DEF.Table2$Perc_perfect<-100*(DEF.Table2$perfect_reads/DEF.Table2$TOTAL)
  
  #### SAVE MM in barcode ----
  
  filename_3<-paste("MM_in_BC_",type,".txt", sep='')
  
  write.table(DEF.Table2, 
              file=filename_3, sep="\t", quote=F, row.names = F)
  
  
  #### Table 3 reads with MM ELSEWHERE ----
  
  summary(as.factor(summary$tag_matches))
  
  summary.dt<-data.table(summary, key="seq_name_bc")
  
  NEW<-as.data.frame(summary.dt[, .N ,by = .(missmatched_read_bases,seq_name_bc)])
  
  Interaction_table_3<-as.data.frame(cbind(paste(NEW$missmatched_read_bases,
                                                 NEW$seq_name_bc, sep='.'),
                                           NEW$N))  
  
  colnames(Interaction_table_3)<-c("Var1","Freq")
  
  #### cSplit
  
  Interaction_table_3<-cSplit(Interaction_table_3,"Var1", 
                              sep = '.', fixed =T, drop=T)
  
  colnames(Interaction_table_3)[which(colnames(Interaction_table_3) == "Var1_1")]<-"mutations_element"
  
  colnames(Interaction_table_3)[which(colnames(Interaction_table_3) == "Var1_2")]<-"seq_name_bc"
  
  
  #### Loop #1 Processing mutation string to a more readable format
  
  Gather_perfect<-data.frame(matrix(vector(), 0, 2,
                                    dimnames=list(c(),c("seq_name_bc","reads_Element_perfect")
                                    )),
                             stringsAsFactors=F)
  
  Gather_mutated<-data.frame(matrix(vector(), 0, 3,
                                    dimnames=list(c(),c("seq_name_bc",
                                                        "reads_Element_mutated",
                                                        "mutations")
                                    )),
                             stringsAsFactors=F)
  
  for(i in 1:length(Interaction_table_3$seq_name_bc))
  {
    seq_name_bc_sel<-as.character(Interaction_table_3$seq_name_bc[i])
    Tag_sel<-as.character(Interaction_table_3$mutations_element[i])
    
    if(Tag_sel == "")
    {
      perfect_reads<-as.numeric(Interaction_table_3$Freq[i])
      
      A<-as.data.frame(cbind(seq_name_bc_sel,perfect_reads))
      
      colnames(A)<-colnames(Gather_perfect)
      
      Gather_perfect<-rbind(Gather_perfect,A)
      
    }else{
      
      mutated_reads<-as.numeric(Interaction_table_3$Freq[i])
      mutation<-Tag_sel
      
      B<-as.data.frame(cbind(seq_name_bc_sel,mutated_reads,Tag_sel))
      
      colnames(B)<-colnames(Gather_mutated)
      
      Gather_mutated<-rbind(Gather_mutated,B)
      
    }
    
  }
  
  
  DEF.Table3<-merge(Gather_mutated,
                    Gather_perfect,
                    by="seq_name_bc",
                    all=T)
  
  DEF.Table3$reads_Element_mutated<-as.numeric(as.character(DEF.Table3$reads_Element_mutated))
  DEF.Table3$reads_Element_mutated[which(is.na(DEF.Table3$reads_Element_mutated))]<-0
  
  DEF.Table3$reads_Element_perfect<-as.numeric(as.character(DEF.Table3$reads_Element_perfect))
  DEF.Table3$reads_Element_perfect[which(is.na(DEF.Table3$reads_Element_perfect))]<-0
  
  DEF.Table3$mutations<-as.character(DEF.Table3$mutations)
  DEF.Table3$mutations[which(is.na(DEF.Table3$mutations))]<-"NONE"
  
  
  DEF.Table3$TOTAL<-DEF.Table3$reads_Element_mutated+DEF.Table3$reads_Element_perfect
  
  DEF.Table3$Perc_mutated<-100*(DEF.Table3$reads_Element_mutated/DEF.Table3$TOTAL)
  
  DEF.Table3$Perc_perfect<-100*(DEF.Table3$reads_Element_perfect/DEF.Table3$TOTAL)
  
  #### Processing mutation string to a summarized table
  
  DEF.Table3_REP<-cSplit(DEF.Table3,"mutations", 
                         sep = ';', fixed =T, drop=F)
  
  
  mini.indexes<-c(grep("seq_name_bc", colnames(DEF.Table3_REP)),
                  which(colnames(DEF.Table3_REP) == "mutations"),
                  grep("mutations_", colnames(DEF.Table3_REP)))
  
  
  Gather_MUT_STRING<-data.frame(matrix(vector(), 0, 3,
                                       dimnames=list(c(),c("seq_name_bc","mutations",
                                                           "mutations_STRING")
                                       )),
                                stringsAsFactors=F)
  
  for(i in 1:length(DEF.Table3_REP$seq_name_bc))
  {
    mini_table_sel<-as.data.frame(DEF.Table3_REP)[i,mini.indexes]
    
    seq_name_bc_sel<-as.character(mini_table_sel$seq_name_bc)
    mutation_sel<-as.character(mini_table_sel$mutations)
    
    
    mini_table_sel.t<-as.data.frame(t(mini_table_sel[,-c(1:2)]))
    colnames(mini_table_sel.t)<-"individual_mutations"
    
    mini_table_sel.t$mutations<-row.names(mini_table_sel.t)
    
    mini_table_sel.t.NO.NA<-mini_table_sel.t[!is.na(mini_table_sel.t$individual_mutations),]
    
    
    ## iterate within the mutations
    
    mutations_vector<-as.character(mini_table_sel.t.NO.NA$individual_mutations)
    
    SUPER_STRING<-NULL
    
    for(k in 1:length(mutations_vector))
    {
      mutations_vector_sel<-mutations_vector[k]
      
      if(sum(grep("del",mutations_vector_sel)) >0)
      {
        a0<-gsub("del","",mutations_vector_sel, perl=T)
        
        if(sum(grep("-",a0)) >0)
        {
          a1<-gsub("\\d+","",a0, perl=T)
          
          a1<-gsub("-","DEL:",a1, perl=T)
          
        }else{
          
          a1<-gsub("\\d+","DEL:",a0, perl=T)
          
        }
        
        a2<-gsub("[atgc]","",a0, perl=T)
        
        partial_string<-paste(a2,a1, sep="_")
        
        SUPER_STRING[k]<-partial_string
        
      }else{
        
        if(sum(grep("ins",mutations_vector_sel)) >0)
        {
          a0<-gsub("ins","",mutations_vector_sel, perl=T)
          
          a1<-gsub("\\d+","INS:",a0, perl=T)
          
          a2<-gsub("[atgc]","",a0, perl=T)
          
          partial_string<-paste(a2,a1, sep="_")
          
          SUPER_STRING[k]<-partial_string
          
          
        }else{
          
          if(sum(grep("NONE",mutations_vector_sel)) >0)
          {
            SUPER_STRING[k]<-"NONE"
            
          }else{
            
            a1<-gsub("\\d+",">",mutations_vector_sel, perl=T)
            
            a2<-gsub("[atgc]","",mutations_vector_sel, perl=T)
            
            partial_string<-paste(a2,a1, sep="_")
            
            SUPER_STRING[k]<-partial_string  
          }
        }
      }# first else
    }# k
    
    SUPER_STRING<-paste(SUPER_STRING, collapse=";")
    A_STRING<-as.data.frame(cbind(seq_name_bc_sel,mutation_sel,SUPER_STRING))
    colnames(A_STRING)<-colnames(Gather_MUT_STRING)
    
    Gather_MUT_STRING<-rbind(Gather_MUT_STRING,A_STRING)
    
    
  }# i
  
  DEF.Table4<-merge(DEF.Table3,
                    Gather_MUT_STRING,
                    by=c("seq_name_bc","mutations"),
                    all=T)
  
  ###### SAVE  -----
  
  filename_4<-paste("MM_elsewhere_",type,".txt", sep='')
  
  write.table(DEF.Table4, 
              file=filename_4, sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--Eugene_summary_1"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Eugene_summary_2"), type="character", default=NULL, 
                metavar="FILE.txt", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
    
  )
  parser = OptionParser(usage = "124_MPRA_PE150_Rscript_edition.R
                        --INPUT_TABLE FILE.txt
                        --Eugene_summary_1 FILE.txt 
                        --Eugene_summary_2 FILE.txt 
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  create_intermediate_files_1(opt)
  
}


###########################################################################

system.time( main() )