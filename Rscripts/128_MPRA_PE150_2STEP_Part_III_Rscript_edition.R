udunits_dir <- file.path(Sys.getenv("HOME"), "udunits")
dyn.load(paste0(udunits_dir, "/local/lib/libudunits2.so.0"))

suppressMessages(library("plyr", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("data.table", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("Biostrings", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("optparse", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("splitstackshape", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs"))
suppressMessages(library("data.table", lib.loc = "/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("udunits2", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))
suppressMessages(library("ggforce", lib.loc="/nfs/users/nfs_m/mt19/sOFTWARE/R_libs/"))


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
  
  #### READ Shendure ----
  
  Shendure = read.table(opt$Shendure, sep="\t", stringsAsFactors = F, header = T)
  
  cat("Manuel_1.5_\n")
  cat(sprintf(as.character(Shendure[1,])))
  cat("\n")
  
  # cat("Manuel_1.5_\n")
  # cat(sprintf(as.character(Shendure$summary$element[1])))
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
  
  #### Eliminate duplicate seq-bcs from the file ----
  
  DUP.examples<-Shendure[duplicated(Shendure$seq_name_bc),]
  
  
  DUP<-Shendure[which(Shendure$seq_name_bc%in%DUP.examples$seq_name_bc),]
  
  Shendure.dt<-data.table(Shendure, key="seq_name_bc")
  colnames(Shendure.dt)
  
  Shendure.dt.max<-Shendure.dt[,.SD[which.max(Perc.SNP.affected)], by=seq_name_bc]
  
  sum(duplicated(Shendure.dt.max$seq_name_bc)) # now it is 0
  
  
  Shendure.NO.DUP<-as.data.frame(Shendure.dt.max)
  
  #### Real Tile including INDEX ----
  
  Shendure.NO.DUP$Real_Tile<-paste(Shendure.NO.DUP$KEY,
                                   Shendure.NO.DUP$Tile,
                                   Shendure.NO.DUP$Index,
                                   Shendure.NO.DUP$Type,
                                   sep=";")
  
  #### First discretization of data thresholds for mutated reads "(11,50]"="Q50" ----
  
  Shendure.NO.DUP$Q90<-cut(Shendure.NO.DUP$Perc_Mutated, breaks = c(-Inf,11,50,75,99.99,Inf))
  
  Shendure.NO.DUP$Q90<-factor(Shendure.NO.DUP$Q90,
                              levels=c("(-Inf,11]","(11,50]","(50,75]","(75,100]","(100, Inf]"),
                              ordered = T)
  
  
  Shendure.NO.DUP$Q90_R<-revalue(Shendure.NO.DUP$Q90,
                                 c("(-Inf,11]"="Q90",
                                   "(11,50]"="Q50",
                                   "(50,75]"="Q25",
                                   "(75,100]"="Q1",
                                   "(100, Inf]"="Q0"))
  
  
  Shendure.NO.DUP$Q90_R<-factor(Shendure.NO.DUP$Q90_R,
                                levels=c("Q90","Q50","Q25","Q1","Q0"),
                                ordered = T)
  
  #### Second discretization of the data mutations affecting SNP of interest "(-Inf,11]"="SNP90" ----
  
  Shendure.NO.DUP$SNP_90<-cut(Shendure.NO.DUP$Perc.SNP.affected, breaks = c(-Inf,11,Inf))
  
  
  Shendure.NO.DUP$SNP_90<-factor(Shendure.NO.DUP$SNP_90,
                                 levels=c("(-Inf,11]","(11, Inf]"),
                                 ordered = T)
  
  Shendure.NO.DUP$SNP_90_R<-revalue(Shendure.NO.DUP$SNP_90,
                                    c("(-Inf,11]"="SNP90",
                                      "(11, Inf]"="SNP10"))
  
  #### Third discretization reads affecting BC  "(-Inf,11]"="BC90" ----
  
  
  Shendure.NO.DUP$BC_90<-cut(Shendure.NO.DUP$Perc.BC.affected, breaks = c(-Inf,11,Inf))
  
  
  Shendure.NO.DUP$BC_90<-factor(Shendure.NO.DUP$BC_90,
                                levels=c("(-Inf,11]","(11, Inf]"),
                                ordered = T)
  
  Shendure.NO.DUP$BC_90_R<-revalue(Shendure.NO.DUP$BC_90,
                                   c("(-Inf,11]"="BC90",
                                     "(11, Inf]"="BC10"))
  
  #### Interactions between Q,SNP and BC ----
  
  Shendure.NO.DUP$interaction_1<-interaction(Shendure.NO.DUP$Q90_R,
                                             Shendure.NO.DUP$SNP_90_R,
                                             drop=F)
  
  Shendure.NO.DUP$interaction_2<-interaction(Shendure.NO.DUP$interaction_1,
                                             Shendure.NO.DUP$BC_90_R,
                                             drop=F)
  
  #### SAVE captured and dropout ----
  
  filename_1<-paste("NO_DUP_BC",type,".txt", sep='')

  write.table(Shendure.NO.DUP,
              file=filename_1, sep="\t", quote=F, row.names = F)
  
  
}

opt = NULL

opt = NULL

create_intermediate_files_2 = function(option_list)
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
  
  #### READ Shendure ----
  
  filename_1<-paste("NO_DUP_BC",type,".txt", sep='')
  
  Shendure.NO.DUP<-read.table(file=filename_1,
                              sep="\t",
                              header=T,
                              stringsAsFactors = F)
  
  Shendure.NO.DUP$interaction_2_R<-factor(Shendure.NO.DUP$interaction_2,
                                          levels=c("Q90.SNP90.BC90","Q50.SNP90.BC90","Q25.SNP90.BC90","Q1.SNP90.BC90","Q0.SNP90.BC90",
                                                   "Q90.SNP90.BC10","Q50.SNP90.BC10","Q25.SNP90.BC10","Q1.SNP90.BC10","Q0.SNP90.BC10",
                                                   "Q90.SNP10.BC90","Q50.SNP10.BC90","Q25.SNP10.BC90","Q1.SNP10.BC90","Q0.SNP10.BC90",
                                                   "Q90.SNP10.BC10","Q50.SNP10.BC10","Q25.SNP10.BC10","Q1.SNP10.BC10","Q0.SNP10.BC10"),
                                          ordered = T)
  
  ####  Real TILE condensation  ----
  
  cat("Manuel_CALM_1_\n")
  str(Shendure.NO.DUP)
  cat("\n")
  
  TILE.condensation.dt<-data.table(Shendure.NO.DUP, key="Real_Tile")
  
  TILE.condensation<-as.data.frame(TILE.condensation.dt[,.(Tiles_recovered=.N,
                                                           Tiles_Q_string=paste(interaction_2_R, collapse="|"),
                                                           Tiles_mut_reads=paste(Mut_total_reads, collapse="|"),
                                                           Tiles_wt_reads=paste(WT_reads, collapse="|")), by="Real_Tile"])
  
  
  
  #### Loop to get the summarizing factors for a TILE ----
  
  Gather_TILE<-data.frame(matrix(vector(), 0,5,
                                 dimnames=list(c(),c("Real_Tile",
                                                     "Q90_Freq","Q90_Reads","Q25_Freq","Q25_Reads")
                                 )),
                          stringsAsFactors=F)
  
  for(i in 1:length(TILE.condensation$Real_Tile))
  {
    
    Real_Tile_SEL<-as.character(TILE.condensation$Real_Tile[i])
    
    midi_dt<-as.data.frame(cbind(as.character(unlist(strsplit(TILE.condensation$Tiles_Q_string[i], "\\|"))),
                                 as.numeric(unlist(strsplit(TILE.condensation$Tiles_mut_reads[i], "\\|"))),
                                 as.numeric(unlist(strsplit(TILE.condensation$Tiles_wt_reads[i], "\\|")))),
                           stringsAsFactors = F)
    
    colnames(midi_dt)<-c("Q_TILE","Mut_reads","WT_reads")
    
    
    midi_dt$Factor3<-factor(midi_dt$Q_TILE,
                            levels=levels(Shendure.NO.DUP$interaction_2_R),
                            ordered=T)
    
    midi_dt$Mut_reads<-as.numeric(midi_dt$Mut_reads)
    midi_dt$WT_reads<-as.numeric(midi_dt$WT_reads)
    
    
    # Summarize Q90 and Q25 Freq and reads
    
    midi_dt.dt<-data.table(midi_dt, key="Factor3")
    
    
    midi_dt.dt_condensed<-as.data.frame(midi_dt.dt[,
                                                   .(Freq=.N,
                                                     WT_reads_c=sum(WT_reads),
                                                     Mut_reads_c=sum(Mut_reads)),
                                                   by=Factor3])
    
    Q90<-levels(Shendure.NO.DUP$interaction_2_R)[1]
    
    Q90_Freq<-0
    Q90_Reads<-0
    
    if(length(midi_dt.dt_condensed$Freq[which(midi_dt.dt_condensed$Factor3 == Q90)]) !=0)
    {
      Q90_Freq<-midi_dt.dt_condensed$Freq[which(midi_dt.dt_condensed$Factor3 == Q90)]
      Q90_Reads<-midi_dt.dt_condensed$WT_reads_c[which(midi_dt.dt_condensed$Factor3 == Q90)]
    }
    
    Q25<-c(levels(Shendure.NO.DUP$interaction_2_R)[3:5])
    
    Q25_Freq<-0
    Q25_Reads<-0
    
    if(length(midi_dt.dt_condensed$Freq[which(midi_dt.dt_condensed$Factor3%in%Q25)]) !=0)
    {
      Q25_Freq<-sum(midi_dt.dt_condensed$Freq[which(midi_dt.dt_condensed$Factor3%in%Q25)])
      Q25_Reads<-sum(midi_dt.dt_condensed$Mut_reads_c[which(midi_dt.dt_condensed$Factor3%in%Q25)])
    }
    
    
    A<-as.data.frame(cbind(Real_Tile_SEL,Q90_Freq,Q90_Reads,Q25_Freq,Q25_Reads))
    
    colnames(A)<-colnames(Gather_TILE)
    
    Gather_TILE<-rbind(Gather_TILE,A)
    
  } #i
  
  #### Categorize Gather_TILE ----
  
  cat("Manuel_CALM_\n")
  str(Gather_TILE)
  cat("\n")
  
  Gather_TILE$Real_Tile<-as.character(Gather_TILE$Real_Tile)
  Gather_TILE$Q90_Freq<-as.numeric(as.character(Gather_TILE$Q90_Freq))
  Gather_TILE$Q90_Reads<-as.numeric(as.character(Gather_TILE$Q90_Reads))
  Gather_TILE$Q25_Freq<-as.numeric(as.character(Gather_TILE$Q25_Freq))
  Gather_TILE$Q25_Reads<-as.numeric(as.character(Gather_TILE$Q25_Reads))
  
  
  
  Gather_TILE$Q90_Freq_R<-cut(Gather_TILE$Q90_Freq, 
                              breaks = c(-Inf,4,10,Inf))
  Gather_TILE$Q25_Freq_R<-cut(Gather_TILE$Q25_Freq, 
                              breaks = c(-Inf,2,4,10,Inf))
  
  Gather_TILE$Interaction_1<-interaction(Gather_TILE$Q90_Freq_R,
                                         Gather_TILE$Q25_Freq_R, drop=F)
  
  Gather_TILE$Interaction_1_R<-factor(Gather_TILE$Interaction_1,
                                      levels=c("(10, Inf].(-Inf,2]","(10, Inf].(2,4]","(10, Inf].(4,10]","(10, Inf].(10, Inf]",
                                               "(4,10].(-Inf,2]","(4,10].(2,4]","(4,10].(4,10]","(4,10].(10, Inf]",
                                               "(-Inf,4].(-Inf,2]","(-Inf,4].(2,4]","(-Inf,4].(4,10]","(-Inf,4].(10, Inf]"),
                                      ordered=T)
  
  Gather_TILE$Interaction_1_R_Transformed<-revalue(Gather_TILE$Interaction_1_R, 
                                                   c("(4,10].(-Inf,2]" = "10 > Q90 >= 5 |Q25 <= 2",
                                                     "(4,10].(2,4]" = "10 > Q90 >= 5 | 2 < Q25 <= 4",
                                                     "(4,10].(4,10]" = "10 > Q90 >= 5 | 5 <= Q25 < 10",
                                                     "(10, Inf].(-Inf,2]" = "Q90 > 10 | Q25 <= 4",
                                                     "(10, Inf].(2,4]" = "Q90 > 10 | Q25 <= 4",
                                                     "(-Inf,4].(4,10]" = "4 >= Q90 | 5 <= Q25 < 10",
                                                     "(-Inf,4].(2,4]" = "FEW Tiles detected",
                                                     "(-Inf,4].(-Inf,2]" = "FEW Tiles detected",
                                                     "(10, Inf].(4,10]" = "Problem too many tiles",
                                                     "(10, Inf].(10, Inf]" = "Problem too many tiles",
                                                     "(4,10].(10, Inf]" = "Problem too many mutated",
                                                     "(-Inf,4].(10, Inf]" = "Problem too many mutated"))
  
  Gather_TILE$Interaction_1_R_Transformed_R<-factor(Gather_TILE$Interaction_1_R_Transformed,
                                                    levels=c("Q90 > 10 | Q25 <= 4","10 > Q90 >= 5 |Q25 <= 2","10 > Q90 >= 5 | 2 < Q25 <= 4",
                                                             "10 > Q90 >= 5 | 5 <= Q25 < 10","4 >= Q90 | 5 <= Q25 < 10",
                                                             "Problem too many mutated","FEW Tiles detected","Problem too many tiles"),
                                                    ordered=T)
  
  #### Merge Gather_TILE with  Real Tile INPUT Table ----
  
  indexes<-c(
    which(colnames(INPUT_table) == "seq_name_bc"),
    which(colnames(INPUT_table) == "barcode"),
    which(colnames(INPUT_table) == "START_TRIMMED"),
    which(colnames(INPUT_table) == "RELPOS_TRIMMED"),
    which(colnames(INPUT_table) == "END_TRIMMED"))
  
  Real_Tile_INPUT_table<-unique(INPUT_table[,-c(indexes)])
  
  Real_Tile_INPUT_table$Real_Tile<-paste(Real_Tile_INPUT_table$KEY,
                                         Real_Tile_INPUT_table$Tile,
                                         Real_Tile_INPUT_table$Index,
                                         Real_Tile_INPUT_table$Type,sep=";")
  
  
  DEF_Tile<-merge(Real_Tile_INPUT_table,
                  Gather_TILE,
                  by="Real_Tile",
                  all=T)
  
  #### SAVE captured and dropout ----
  
  filename_2<-paste("DEF_Tile_",type,".txt", sep='')
  
  write.table(DEF_Tile,
              file=filename_2, sep="\t", quote=F, row.names = F)
  
}

create_intermediate_files_3 = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
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
  
  
  #### READ DEF_TILE ----
  
  filename_2<-paste("DEF_Tile_",type,".txt", sep='')
  
  cat("Manuel_SUPER_CALM_\n")
  str(filename_2)
  cat("\n")
  
  DEF_Tile<-read.table(file=filename_2,
                              sep="\t",
                              header=T,
                              stringsAsFactors = F)
  
  cat("Manuel_SUPER_CALM_2_\n")
  str(DEF_Tile)
  cat("\n")
  
  #### dropout tiles ----
  
  drop_out_reads<-DEF_Tile[which(DEF_Tile$Q25_Reads <=10 &
                                   DEF_Tile$Q90_Reads <=10),]
  
  #### DEF Tile condense Allele Ratios ----
  
  DEF_Tile_AR<-DEF_Tile[-which(DEF_Tile$Real_Tile%in%drop_out_reads$Real_Tile),]
  
  DEF_Tile_AR$Tile_AR<-paste(DEF_Tile_AR$KEY,
                             DEF_Tile_AR$Tile,
                             sep=";")
  
  DEF_Tile_AR.dt<-data.table(DEF_Tile_AR, key="Tile_AR")
  
  DEF_Tile_AR_CONDENSED<-as.data.frame(DEF_Tile_AR.dt[,.(Alleles=.N,
                                                         Type_c=paste(Type, collapse="__"),
                                                         Q90_Reads_c=paste(Q90_Reads, collapse="__"),
                                                         Q25_Reads_c=paste(Q25_Reads, collapse="__"),
                                                         CLASS=paste(Interaction_1_R_Transformed_R, collapse="__")), 
                                                      by=Tile_AR])
  ##### LOOP for Allele Ratios ----
  
  Gather_Allele_Ratios<-data.frame(matrix(vector(), 0,6,
                                          dimnames=list(c(),c("Tile_AR","Alleles","CLASS",
                                                              "Q90_REF_ALT_Log_Ratio",
                                                              "Q25_REF_ALT_Log_Ratio",
                                                              "Global_Log_Ratio")
                                          )),
                                   stringsAsFactors=F)
  
  for(i in 1:length(DEF_Tile_AR_CONDENSED$Tile_AR))
  {
    DEF_Tile_SEL<-as.data.frame(DEF_Tile_AR_CONDENSED[i,])
    
    Tile_AR_SEL<-as.character(DEF_Tile_SEL$Tile_AR)
    Alleles_SEL<-as.numeric(as.character(DEF_Tile_SEL$Alleles))
    CLASS_SEL<-unlist(strsplit(DEF_Tile_SEL$CLASS, "__"))
    
    CLASS_DEF<-NULL
    
    if(length(CLASS_SEL) > 2)
    {
      CLASS_DEF<-"Too_many_tiles"
      
    }else{
      
      if(length(CLASS_SEL) == 1)
      {
        CLASS_DEF<-"Too_few_tiles"
        
      }else{
        
        if(CLASS_SEL[1] == CLASS_SEL[2])
        {
          CLASS_DEF<-CLASS_SEL[1]
        }else{
          
          if(CLASS_SEL[1] == "FEW Tiles detected")
          {
            CLASS_DEF<-"FEW Tiles detected_ALT"
            
          }else{
            
            if(CLASS_SEL[2] == "FEW Tiles detected")
            {
              CLASS_DEF<-"FEW Tiles detected_REF"
              
            }else{
              
              
              if(CLASS_SEL[2] == "Problem too many mutated")
              {
                CLASS_DEF<-"Problem too many mutated_REF"
                
              }else{
                
                if(CLASS_SEL[1] == "Problem too many mutated")
                {
                  CLASS_DEF<-"Problem too many mutated_ALT"
                  
                }else{
                  
                  if(sum(grep("Q90 > 10",CLASS_SEL[1])) >0 & sum(grep("10 > Q90 >= 5",CLASS_SEL[2])))
                  {
                    CLASS_DEF<-"Excellent_REF_Good_ALT" 
                  }else{
                    
                    if(sum(grep("10 > Q90 >= 5",CLASS_SEL[1])) >0 & sum(grep("10 > Q90 >= 5",CLASS_SEL[2])))
                    {
                      CLASS_DEF<-"Good_REF_Good_ALT" 
                    }else{
                      
                      if(sum(grep("10 > Q90 >= 5",CLASS_SEL[1])) >0 & sum(grep("Q90 > 10",CLASS_SEL[2])))
                      {
                        CLASS_DEF<-"Good_REF_Excellent_ALT" 
                      }else{
                        
                        if(sum(grep("10 > Q90 >= 5",CLASS_SEL[1])) >0 & sum(grep("4 >= Q90",CLASS_SEL[2])))
                        {
                          CLASS_DEF<-"Good_REF_Bad_ALT" 
                        }else{
                          if(sum(grep("4 >= Q90",CLASS_SEL[1])) >0 & sum(grep("10 > Q90 >= 5",CLASS_SEL[2])))
                          {
                            CLASS_DEF<-"Bad_REF_Good_ALT" 
                          }else{
                            CLASS_DEF<-DEF_Tile_SEL$CLASS
                          }
                        }
                        
                        
                      }
                    }
                    
                  }
                }
              }
            }
          }
        }
        
      }
      
    }
    
    
    
    Q90_allele_reads<-unlist(strsplit(DEF_Tile_SEL$Q90_Reads_c, "__"))
    Q90_allele_reads_T<-NULL
    
    for(h in 1:2)
    {
      if(Q90_allele_reads[h] == "NA" |is.na(Q90_allele_reads[h]))
      {
        Q90_allele_reads_T[h]<-0
      }else{
        
        Q90_allele_reads_T[h]<-Q90_allele_reads[h]
      }
      
    }
    
    Q90_allele_reads_T<-as.numeric(Q90_allele_reads_T)
    
    Q90_allele_reads_Log_Ratio<-log((Q90_allele_reads_T[2]+1)/(Q90_allele_reads_T[1]+1))
    
    Q25_allele_reads<-unlist(strsplit(DEF_Tile_SEL$Q25_Reads_c, "__"))
    Q25_allele_reads_T<-NULL
    
    for(h in 1:2)
    {
      if(Q25_allele_reads[h] == "NA" |is.na(Q25_allele_reads[h]))
      {
        Q25_allele_reads_T[h]<-0
      }else{
        
        Q25_allele_reads_T[h]<-Q25_allele_reads[h]
      }
      
    }
    
    Q25_allele_reads_T<-as.numeric(Q25_allele_reads_T)
    
    Q25_allele_reads_Log_Ratio<-log((Q25_allele_reads_T[2]+1)/(Q25_allele_reads_T[1]+1))
    
    TOT_REF<-sum(Q90_allele_reads_T[2],Q25_allele_reads_T[2])
    TOT_ALT<-sum(Q90_allele_reads_T[1],Q25_allele_reads_T[1])
    
    Global_Log_Ratio<-log((TOT_REF+1)/(TOT_ALT+1))
    
    A_SEL<-as.data.frame(cbind(Tile_AR_SEL,Alleles_SEL,CLASS_DEF,
                               Q90_allele_reads_Log_Ratio,
                               Q25_allele_reads_Log_Ratio,
                               Global_Log_Ratio))
    colnames(A_SEL)<-colnames(Gather_Allele_Ratios)
    
    Gather_Allele_Ratios<-rbind(Gather_Allele_Ratios,A_SEL)
  }
  
  
  #### transform Gather_Allele_Ratios ----
  
  Gather_Allele_Ratios$Tile_AR<-as.character(Gather_Allele_Ratios$Tile_AR)
  Gather_Allele_Ratios$Alleles<-as.numeric(as.character(Gather_Allele_Ratios$Alleles))
  Gather_Allele_Ratios$Q90_REF_ALT_Log_Ratio<-as.numeric(as.character(Gather_Allele_Ratios$Q90_REF_ALT_Log_Ratio))
  Gather_Allele_Ratios$Q25_REF_ALT_Log_Ratio<-as.numeric(as.character(Gather_Allele_Ratios$Q25_REF_ALT_Log_Ratio))
  Gather_Allele_Ratios$Global_Log_Ratio<-as.numeric(as.character(Gather_Allele_Ratios$Global_Log_Ratio))
  
  Gather_Allele_Ratios$Tile_AR[is.na(Gather_Allele_Ratios$Q90_REF_ALT_Ratio)]
  
  cat("Manuel_SUPER_CALM_\n")
  str(Gather_Allele_Ratios)
  cat("\n")
  
  #### SAVE captured and dropout ----
  
  filename_3<-paste("Gather_Allele_Ratios_",type,".txt", sep='')
  
  write.table(Gather_Allele_Ratios,
              file=filename_3, sep="\t", quote=F, row.names = F)
  
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
    make_option(c("--Shendure"), type="character", default=NULL, 
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
                        --Shendure FILE.txt 
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  create_intermediate_files_1(opt)
  #graph_1(opt) Abandonned could not make it work
  create_intermediate_files_2(opt)
  create_intermediate_files_3(opt)
  
}


###########################################################################

system.time( main() )