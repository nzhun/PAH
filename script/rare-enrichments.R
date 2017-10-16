#setwd("~/server")

setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")
source("Pipeline/Enrichment/Denovo_enrichment.R")
#### start select individual ###
get_lgd <- function(dat){
  index1 <- grep("stop",dat$VarClass)
  index2 <- grep("^frame",dat$VarClass)
  index3 <- grep("splicing",dat$VarFunc)
  return(unique(c(index1,index2,index3)))
}
get_indel <- function(dat){
  len_ref <- unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt <- unlist(lapply(dat$ALT,FUN = function(x)nchar(x)))
  if(length(len_ref) != length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref != len_alt))
}
get_snv <- function(dat){
  len_ref <- unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt <- unlist(lapply(dat$ALT,FUN = function(x)nchar(x)))
  if(length(len_ref) != length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref==len_alt))
}

cohrts_gene_merge<- function(cohorts,solved.samples,ofile){
  cohorts$solved<-"-"
  for(id in 1:length(cohorts$proband)){ 
    index<-which(solved.samples$ID==cohorts$proband[id]); 
    if(length(index)>0){
      rs="";for(i in index){rs<-paste(solved.samples$GeneName[i],rs,sep=",")}
      cohorts$solved[id]<- rs;
    } 
  }
 # ofile <- "PAH/Result/Data/PAH_ped_adult.list.csv"
  write.csv(cohorts,ofile,row.names = F,quote = T)
}

VQSR_filter <- function(data){
  
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00"
                       , "VQSRTrancheSNP99.80to99.90" ,"VQSRTrancheSNP99.70to99.80",
                       "VQSRTrancheSNP99.60to99.70" ,  "VQSRTrancheSNP99.50to99.60",
                    #   "VQSRTrancheINDEL99.50to99.60",
                    #   "VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00"
                    )
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  return(data)
}


filter_allfreq_CHD <- function(data,freq_avg,freq_max){

  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                     &as.numeric(data$gnomAD_genome_EAS)<freq_max
                     &as.numeric(data$gnomAD_genome_AFR)<freq_max
                     &as.numeric(data$gnomAD_genome_NFE)<freq_max
                    # &as.numeric(data$gnomAD_genome_FIN)<freq_max
                     &as.numeric(data$gnomAD_genome_OTH)<freq_max
                     &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     &as.numeric(data$gnomAD_genome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_ALL)<freq_avg
                      &as.numeric(data$gnomAD_exome_EAS)<freq_max
                      &as.numeric(data$gnomAD_exome_NFE)<freq_max
                  #    &as.numeric(data$gnomAD_exome_FIN)<freq_max
                     &as.numeric(data$gnomAD_exome_OTH)<freq_max
                      &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                      &as.numeric(data$gnomAD_exome_AMR)<freq_max
                  &as.numeric(data$gnomAD_exome_AFR)<freq_max
  
 
  ),]
  
  return(data)
}

filter_allfreq_PAH <- function(data,freq_avg,freq_max){
  
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$`1000g2015aug_all`)< freq_max
                     & as.numeric(dat$esp6500siv2_all)<freq_max
                   #  &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                 #    &as.numeric(data$ExAC_OTH)< freq_max
                  #   &as.numeric(data$gnomAD_genome_EAS)<freq_max
                   #  &as.numeric(data$gnomAD_genome_AFR)<freq_max
                  #   &as.numeric(data$gnomAD_genome_NFE)<freq_max
                     # &as.numeric(data$gnomAD_genome_FIN)<freq_max
                  #   &as.numeric(data$gnomAD_genome_OTH)<freq_max
                  #   &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                #     &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     &as.numeric(data$gnomAD_genome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_ALL)<freq_avg
                #     &as.numeric(data$gnomAD_exome_EAS)<freq_max
                #     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     #    &as.numeric(data$gnomAD_exome_FIN)<freq_max
                #     &as.numeric(data$gnomAD_exome_OTH)<freq_max
                #     &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                #     &as.numeric(data$gnomAD_exome_AMR)<freq_max
                #     &as.numeric(data$gnomAD_exome_AFR)<freq_max
                     
                     
  ),]
  
  return(data)
}



filter_allfreq2 <- function(data){
  freq <- 0.0001
  freq2 <- 0.0001
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq2)
                     &as.numeric(data$ExAC_AFR)< freq2
                     &as.numeric(data$ExAC_NFE)< freq2
                     &as.numeric(data$ExAC_SAS)< freq2
                     &as.numeric(data$ExAC_EAS)< freq2
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  data <- data[which( as.numeric(data$AC)< 25
                      &as.numeric(data$AB)>0.2
  ),]
 
  return(data)
}

# map_filter<-function(data){
#   index<-grep("mappability",names(data),ignore.case = T)
#   if(length(index)>0){ 
#     data <- data[which(as.numeric(data[,index])==1),]
#   }
#   
#   if(length(grep("genomicSuperDups",names(data)))>0){
#     index<-grep("Score",data$genomicSuperDups)
#   
#     dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;"))[2])))>0.95)]
#     if(length(dup_indexs)>0){
#       data <- data[-dup_indexs,]
#     }
#   }
#   return (data)
# }

# VQSR_filter <- function(data){
# 
#   if(length(grep("VQSLOD",names(data)))>0){
#        data <- data[which(data$VQSLOD > -2.75),]
#   }
#   return(data)
# }

# COV_filter <- function(data){
#   if(length(grep("GnomAD_Genome_cov",names(data)))>0){
#     data<-data[which(as.numeric(data$GnomAD_Genome_cov10)>0.85),]
#   }
#   if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
#     data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.8 & 2*data$AC_xGen/data$AN_xGen >0.8 & 2*(data$AC_VCR/data$AN_VCR) >0.8  &data$GnomAD_Genome_cov10>0.8 ),]
#   }
#   if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
#     data <- data[which(2*data$AC_Control/data$AN_Control >0.8 & 2*data$AC_VU/data$AN_VU >0.8 ),]
#   }
#   return(data)
# }


Merge_test_single <- function(dat_vcr,dat_xgen, TestName){
  merge_test <- c()
  
  c_capture1  <-  dat_vcr
  c_capture2  <-  dat_xgen
  for (class in unique(c(c_capture1$class,c_capture2$class))){
    index1<-which(c_capture1$class==class)
    index2<-which(c_capture2$class==class)
    obv_case<-0
    obv_control<-0
    len_case<-0
    len_control<-0;
    ncase<-0
    ncontrol<-0
    
    if(length(index1)>0){
        obv_case  <-  obv_case + c_capture1$observed[index1]
        obv_control  <-  c_capture1$expected[index1]
        len_control  <-  c_capture1$Ncontrol[index1]
        len_case <- len_case + c_capture1$Ncase[which(c_capture1$class==class)]
        ncase  <-  ncase + c_capture1$Carrier[which(c_capture1$class==class)]
        ncontrol  <-  c_capture1$CtrCarrier[which(c_capture1$class==class)]
    }
    
    if(length(index2)>0){
      obv_case <- obv_case + c_capture2$observed[index2] 
      obv_control  <-  c_capture2$expected[index2]
      len_control  <-  c_capture2$Ncontrol[index2]
      len_case <- len_case + c_capture2$Ncase[which(c_capture2$class==class)] 
      ncase  <- ncase +c_capture2$Carrier[which(c_capture2$class==class)] 
      ncontrol  <-  c_capture2$CtrCarrier[which(c_capture2$class==class)]
    }
    
     
   # len_case  <-  c_capture1$Ncase[which(c_capture1$class==class)]+c_capture2$Ncase[which(c_capture2$class==class)] 
  #  len_control  <-  c_capture1$Ncontrol[which(c_capture1$class==class)]
    
 #   ncase  <-  c_capture1$Carrier[which(c_capture1$class==class)]+c_capture2$Carrier[which(c_capture2$class==class)] 
#    ncontrol  <-  c_capture1$CtrCarrier[which(c_capture1$class==class)]
    bt<-Btest(case = obv_case,case_total = len_case,control = obv_control,control_total = len_control,correct = 1)
    enrich_mis  <-bt$estimate  #call_enrichment_n(obv_case,obv_control,len_case,len_control)
    pvalue  <-bt$p.value #  call_pvalue(ncase,len_case,ncontrol,len_control)
   # pvalue=pvalue$p.value
    Effect=(1-1/enrich_mis)*ncase
    test  <-  data.frame(class=class,observed=obv_case,expected=obv_control,
                         enrichment=enrich_mis,pvalue=pvalue,
                         Ncase=len_case,Ncontrol=len_control,
                         Effect=Effect,Carrier=ncase,CtrCarrier=ncontrol,
                         type=TestName,Test=TestName,stringsAsFactors = F)
    
    merge_test  <-  rbind(merge_test,test)
  }
  
  return(merge_test)
  
}


Merge_test_group  <-  function(dat_vcr,dat_xgen, TestName){
  merge_test  <-  c()
  for(type in unique(dat_vcr$type)){
    c_capture1  <-  dat_vcr[which(dat_vcr$type==type),]
    c_capture2  <-  dat_xgen[which(dat_xgen$type==type),]
    for (class in unique(c_capture1$class[which(c_capture1$type==type)])){
      obv_case  <-  c_capture1$observed[which(c_capture1$class==class)]+c_capture2$observed[which(c_capture2$class==class)] 
      obv_control  <-  c_capture1$expected[which(c_capture1$class==class)]
      
      len_case  <-  c_capture1$Ncase[which(c_capture1$class==class)]+c_capture2$Ncase[which(c_capture2$class==class)] 
      len_control  <-  c_capture1$Ncontrol[which(c_capture1$class==class)]
      
      ncase  <-  c_capture1$Carrier[which(c_capture1$class==class)]+c_capture2$Carrier[which(c_capture2$class==class)] 
      ncontrol  <-  c_capture1$CtrCarrier[which(c_capture1$class==class)]
      enrich_mis  <-  call_enrichment_n(obv_case,obv_control,len_case,len_control)
      pvalue  <-  call_pvalue(ncase,len_case,ncontrol,len_control)
      pvalue=pvalue$p.value
      Effect=(1-1/enrich_mis)*ncase
      test  <-  data.frame(class=class,observed=obv_case,expected=obv_control,
                           enrichment=enrich_mis,pvalue=pvalue,
                           Ncase=len_case,Ncontrol=len_control,
                           Effect=Effect,Carrier=ncase,CtrCarrier=ncontrol,
                           type=type,Test=TestName,stringsAsFactors = F)
      
      merge_test  <-  rbind(merge_test,test)
    }
  }
  return(merge_test)
}



#### associated gene ### 

#pah.list<-read.csv("PAH/Result/Data/PAH_pediatric_adult_list.csv",header = 1,stringsAsFactors = F,strip.white = T,check.names = F)
#pah_solved_vars<-read.csv("PAH/Result/Data/PAH_ped_adult_All.known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")






geneset_enrich<-function(control_data,xgen_data,vcr_data,gene_asso,xgen_size,vcr_size,control_size){
  correction<-1
  correction2<-1
  sets <- load_dataset()
  exac <- load_exac()
  lib <- load_genesets(sets,exac)
  init_all<-NA
  init<-NA
  enrich_asso<-NA
  enrich_minor_asso<-NA
  enrich_asso_VCR<-NA
  enrich_minor_asso_VCR<-NA
  df<-NA
  merged_known_test <-NA
  merged_minor_known_test<-NA
  df_unsolved_xgen<-NA
  df_unsolved_VCR <-NA
  df_vcr<-NA
  merge_test_unsolved_all<-NA

  
  
  if(dim(xgen_data)[1]>0){
    init_all <- enrich_func(xgen_data,control_data,"All",F,c(),1,1,xgen_size,control_size)
    init_all$type="General_beforeCorrecting"
    init_all$Test="ALL"
    correction <- init_all$enrichment[which(init_all$class=="SYN")]
    correction2 <- length(get_indel(xgen_data))/length(unique(xgen_data$ProbandName))/(length(get_indel(control_data))/length(unique(control_data$ProbandName)))
    
  #  init <- enrich_func(xgen_data,control_data,"All",F,c(),correction,correction2)
  #  init$type <- "General"
  #  init$Test="ALL"
    enrich_asso <- enrich_func(xgen_data,control_data,gene_asso,F,c(),correction,correction2,xgen_size,control_size)
    if(dim(enrich_asso)[1]>0){
      enrich_asso$type = "XGEN_knownGene"
      enrich_asso$Test="XGEN_KnownGene"
    }

  }
  if(dim(vcr_data)[1]>0){
    enrich_asso_VCR <- enrich_func(vcr_data,control_data,gene_asso,F,c(),correction,correction2,vcr_size,control_size)
    if(dim(enrich_asso_VCR)[1]>0){
      enrich_asso_VCR$type="VCR_knownGene"
      enrich_asso_VCR$Test="VCR_KnownGene"
    }
  }
  
  merged_known_test<-c()
 
  if(dim(vcr_data)[1]>0  && dim(xgen_data)[1]>0 ){
    if( dim(enrich_asso_VCR)[1]>0 && dim(enrich_asso)[1]>0){
      merged_known_test <- Merge_test_single(enrich_asso_VCR,enrich_asso,"merged_known_gene")
    }else if( dim(enrich_asso_VCR)[1]>0){
      merged_known_test<-enrich_asso_VCR
    }else if(dim(enrich_asso)[1]>0){
      merged_known_test<-enrich_asso
    }
  }
  return(rbind(enrich_asso, enrich_asso_VCR,
                  merged_known_test)
         )
}

enrich_test<-function(xgen_size,vcr_size,control_size,control_data,xgen_data,vcr_data,gene_asso,minor_gene_asso,solved.samples,outfile){
      ### start  do the enrichment ###
      #tdata,t_control,genesets,invert,gexcludes,correction,correction2
    correction<-1
    correction2<-1
    sets <- load_dataset()
    exac <- load_exac()
    lib <- load_genesets(sets,exac)
    init_all<-NA
    init<-NA
    enrich_asso<-NA
    enrich_minor_asso<-NA
    enrich_asso_VCR<-NA
    enrich_minor_asso_VCR<-NA
    df<-NA
    merged_known_test <-NA
    merged_minor_known_test<-NA
    df_unsolved_xgen<-NA
    df_unsolved_VCR <-NA
    df_vcr<-NA
    merge_test_unsolved_all<-NA
    
    merge_test_all<-NA
    
    if(dim(xgen_data)[1]>0){
        init_all <- enrich_func(xgen_data,control_data,"All",F,c(),1,1,xgen_size,control_size)
        init_all$type="General_beforeCorrecting"
        init_all$Test="ALL"
        correction <- init_all$enrichment[which(init_all$class=="SYN")]
        correction2 <- length(get_indel(xgen_data))/length(unique(xgen_data$ProbandName))/(length(get_indel(control_data))/length(unique(control_data$ProbandName)))
      
        init <- enrich_func(xgen_data,control_data,"All",F,c(),correction,correction2,xgen_size,control_size)
        init$type <- "General"
        init$Test="ALL"
        enrich_asso <- enrich_func(xgen_data,control_data,gene_asso,F,c(),correction,correction2,xgen_size,control_size)
        if(dim(enrich_asso)[1]>0){
          enrich_asso$type = "XGEN_knownGene"
          enrich_asso$Test="XGEN_KnownGene"
        }
        enrich_minor_asso <- enrich_func(xgen_data,control_data,minor_gene_asso,F,c(),correction,correction2,xgen_size,control_size)
        if(dim(enrich_minor_asso)[1]>0){
          enrich_minor_asso$type="XGEN_minor_knownGene"
          enrich_minor_asso$Test="XGEN_minor_KnownGene"
        }
        
        df <- enrichSet(xgen_data,control_data,gene_asso,correction,correction2,lib,xgen_size,control_size)  ## gene set 
        if(dim(df)[1]>0){
          df$Test="all-xGen"
        }
        
       # print(xgen_data[which(!xgen_data$ProbandName %in% solved.samples))
        unsolved_xgen_size<- xgen_size-length(unique(xgen_data$ProbandName[which(xgen_data$ProbandName %in%  solved.samples$ID)]))
        
        df_unsolved_xgen <- enrichSet(xgen_data[which(!xgen_data$ProbandName %in% solved.samples$ID),],control_data,gene_asso,correction,correction2,
                                      lib,unsolved_xgen_size,control_size)
        if(dim(df_unsolved_xgen)[1]>0){
          df_unsolved_xgen$Test="all-xGen-unsolved"
        }
      }
      if(dim(vcr_data)[1]>0){
          enrich_asso_VCR <- enrich_func(vcr_data,control_data,gene_asso,F,c(),correction,correction2,vcr_size,control_size)
          if(dim(enrich_asso_VCR)[1]>0){
            enrich_asso_VCR$type="VCR_knownGene"
            enrich_asso_VCR$Test="VCR_KnownGene"
          }
          
          enrich_minor_asso_VCR <- enrich_func(vcr_data[which(!vcr_data$ProbandName %in% solved.samples$ID[solved.samples$GeneName %in% c("BMPR2","TBX4")]),],
                                               control_data,minor_gene_asso,F,c(),correction,correction2,
                                               vcr_size,control_size)
          if(dim( enrich_minor_asso_VCR )[1]>0){
            enrich_minor_asso_VCR$type="VCR_minor_knownGene"
            enrich_minor_asso_VCR$Test="VCR_minor_knownGene"
          }
          df_vcr <- enrichSet(vcr_data,control_data,gene_asso,correction,correction2,lib,
                              vcr_size,control_size)  ## gene set 
          if(dim(df_vcr)>0){df_vcr$Test="all-VCR"}
          
          unsolved_vcr_size<- vcr_size-length(unique(vcr_data$ProbandName[which(vcr_data$ProbandName %in%  solved.samples$ID)]))
          df_unsolved_VCR <- enrichSet(vcr_data[which(!vcr_data$ProbandName %in%  solved.samples$ID),],control_data,gene_asso,1,1,lib,
                                  unsolved_vcr_size,control_size)
          if(dim(df_unsolved_VCR)>0){df_unsolved_VCR$Test <- "all-VCR-unsolved"}
      }
    
      merged_known_test<-c()
      merged_minor_known_test<-c()
      merge_test_unsolved_all<-c()
      if(dim(vcr_data)[1]>0  && dim(xgen_data)[1]>0 ){
         if( dim(enrich_asso_VCR)[1]>0 && dim(enrich_asso)[1]>0){
            merged_known_test <- Merge_test_single(enrich_asso_VCR,enrich_asso,"merged_known_gene")
         }else if( dim(enrich_asso_VCR)[1]>0){
           merged_known_test<-enrich_asso_VCR
         }else if(dim(enrich_asso)[1]>0){
           merged_known_test<-enrich_asso
         }
        if( dim(enrich_minor_asso_VCR)[1]>0 && dim(enrich_minor_asso)[1]>0){
          merged_minor_known_test <- Merge_test_single(enrich_minor_asso_VCR,enrich_minor_asso,"merged_minor_known_gene")
        }else if(dim(enrich_minor_asso_VCR)[1]>0){
          merged_minor_known_test<-enrich_minor_asso_VCR
        }else if(dim(enrich_minor_asso)[1]>0){
          merged_minor_known_test <- enrich_minor_asso
        }
        
        
        ### correction_merged
        if( dim(df_unsolved_VCR)[1]>0 && dim(df_unsolved_xgen)[1]>0){
          merge_test_unsolved_all <- Merge_test_group(df_unsolved_VCR,df_unsolved_xgen,"unsolved_merged_all")
        }else if(dim(df_unsolved_VCR)[1]>0){
          merge_test_unsolved_all<-df_unsolved_VCR
        }else if(dim(df_unsolved_xgen)[1]>0){
          merge_test_unsolved_all<-df_unsolved_xgen
        }
        
        ### correction_merged
        if( dim(df_vcr)[1]>0 && dim(df)[1]>0){
          merge_test_all <- Merge_test_group(df_vcr,df,"merged_all")
        }else if(dim(df_vcr)[1]>0){
          merge_test_all<-df_vcr
        }else if(dim(df)[1]>0){
          merge_test_all<-df
        }
        
      }
      write.csv(rbind(init_all,init,
                      enrich_asso,
                      enrich_minor_asso,
                      df,
                      enrich_asso_VCR,
                      enrich_minor_asso_VCR,
                      df_vcr,
                      merged_known_test,
                      merged_minor_known_test,
                      df_unsolved_xgen,df_unsolved_VCR,merge_test_unsolved_all,merge_test_all
      ),outfile,row.names = F)
      
      print(paste("output to ",outfile,sep=""))
      
}



rare_enrich<-function(cohorts,mutlist,fasso,out_prefix,data,doset){

  known_vars<-read.csv(mutlist,header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
  
  gene_asso_set  <-  read.csv(fasso,header=1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T)
  gene_asso  <-  gene_asso_set[,1]
  minor_gene_asso  <-  setdiff(gene_asso,c("BMPR2","TBX4"))
  
  #unsolved_denovo<-denovo[!which(denovo$ProbandName%in%solved_denovos),]
  solved.samples<-known_vars[,c("proband","Gene.refGene")]
  names(solved.samples)<-c("ID","GeneName")
  
  conf_solved_samples<-solved.samples[which(solved.samples$GeneName%in%gene_asso),]
  
  ### VCR(CU+VU)
  pop_vcr  <-  read.table("PAH/Result/Data/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
  vcr_european<-cohorts$proband[which(cohorts$proband%in%pop_vcr$V1 )]
  
  ###XGEN
  pop <- read.csv("PAH/Result/Data/source/xgen_pop.csv",header=1,stringsAsFactors = F,check.names = F)
  xgen_european <- intersect(cohorts$proband,pop$sampleID[which(pop$pop=="EUR")])
  
  #XGEN_solved<-xgen_european[which(xgen_european%in%conf_solved_samples)]
  
  ### control
  control_IDs  <-  read.csv("PAH/Result/Data/source/european-control.txt",header = F,check.names = F,stringsAsFactors = F)
  control_IDs   <-   control_IDs$V1
  
  ### indel IGV check xGen ###  the indels failed IGV check

  #ped<-read.table()
  
  ### start import variant #####
  
  
  control_dat <- data[which(data$ProbandName %in%  control_IDs),]
  xgen_dat<-data[which(data$ProbandName%in% xgen_european),]
  vcr_dat<-data[which(data$ProbandName%in% vcr_european),]
  
  
  control_size<-length(unique(control_dat$ProbandName))
  vcr_size<-length(unique(vcr_dat$ProbandName))
  xgen_size<-length(unique(xgen_dat$ProbandName))

  ### do geneset test
  if(doset ){
    set0<-gene_asso_set[,1]
    rs <- geneset_enrich(control_dat,xgen_dat,vcr_dat,gene_asso = set0)
    rs$Geneset<-"ALL"
    rs$NGene<-length(set0)
    arr<-rs
    if(length(grep("HHE",names(gene_asso_set)))>0 && length(grep("pLI",names(gene_asso_set)))>0 ){
      set1<-gene_asso_set$Gene[which(gene_asso_set$HHE_Rank>75)]
      set2<-gene_asso_set$Gene[which(gene_asso_set$pLI_Score>=0.9)]
      set3<-intersect(set1,set2)
      he_rs1 <- geneset_enrich(control_dat,xgen_dat,vcr_dat,gene_asso = set1)
      he_rs1$Geneset<-"HHE"
      he_rs1$NGene<-length(set1)
      hpin <- geneset_enrich(control_dat,xgen_dat,vcr_dat,gene_asso = set2)
      hpin$Geneset <-"pLi"
      hpin$NGene <-length(set3)
      hpin_he<-geneset_enrich(control_dat,xgen_dat,vcr_dat,gene_asso = set3)
      hpin_he$Geneset <-"PLI&HHE"
      hpin_he$NGene <-length(set3)
      arr<-rbind(arr,he_rs1,hpin,hpin_he)
    }
    ### genelist2 
    write.csv(arr,file = paste("PAH/Result/Data/",out_prefix,"_knowngeneset.csv",sep=""))
  }
  
  
  outfile=paste("PAH/Result/Data/",out_prefix,"_ALL.csv",sep="")
  enrich_test(xgen_size,vcr_size,control_size,control_dat,xgen_dat,vcr_dat,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  
  # outfile=paste("PAH/Result/Data/",out_prefix,"_child.csv",sep="")
  # child_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE,ignore.case = T)]),]
  # child_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE,ignore.case = T)]),]
  # enrich_test(control_dat,child_xgen_dat,child_vcr_dat,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  # 
  # outfile=paste("PAH/Result/Data/",out_prefix,"_child.csv",sep="")
  # child_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
  # child_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
  # enrich_test(control_dat,child_xgen_dat,child_vcr_dat,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  
  # for( dd in unique(cohorts$Disease)){
  #   outfile=paste("PAH/Result/Data/",dd,"_",out_prefix,"_child.csv",sep="")
  #   child_xgen_dat2<-child_xgen_dat[which(child_xgen_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
  #   child_vcr_dat2<-child_vcr_dat[which(child_vcr_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
  #   enrich_test(control_dat,child_xgen_dat2,child_vcr_dat2,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  # }
  # 
  # outfile=paste("PAH/Result/Data/",out_prefix,"_adult.csv",sep="")
  # adult_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("adult",cohorts$TYPE,ignore.case = T)]),]
  # adult_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("adult",cohorts$TYPE,ignore.case = T)]),]
  # enrich_test(control_dat,adult_xgen_dat,adult_vcr_dat,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  # 
  # for( dd in unique(cohorts$Disease)){
  #   outfile=paste("PAH/Result/Data/",dd,"_",out_prefix,"_child.csv",sep="")
  #   adult_xgen_dat2<-adult_xgen_dat[which(adult_xgen_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
  #   adult_vcr_dat2<-adult_vcr_dat[which(adult_vcr_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
  #   enrich_test(control_dat,adult_xgen_dat2,adult_vcr_dat2,gene_asso,minor_gene_asso,conf_solved_samples,outfile)
  # }
  
}


# 
# IGV_match<-function(){
#   rs<-read.table("PAH/PAH-CHD/PAH-CHD.253KNG.TP.list",header = F,comment.char = "",check.names = F,stringsAsFactors = F)
#   all<-read.csv("PAH/PAH-CHD/Data/PAH-CHD.253KNG.csv",header = 1,comment.char = "",check.names = F,stringsAsFactors = F)
#   ped<-read.table("PAH/Joint_calls_20170715/VCX_Control.ped",header=1,comment.char = "",check.names = F,stringsAsFactors = F,fill=T)
#   rem<-c();
#   
#   for(j in 1:dim(rs)[1]){
#     index<-which(all$`#CHROM`==rs$V2[j] & all$FROM==rs$V3[j])
#     if(length(index)==0){
#       print(paste(rs[j,]))
#     }else{
#       index2<-which(all$`#CHROM`==rs$V2[j] & all$FROM==rs$V3[j] & all$proband==rs$V1[j])  
#       if(length(index2)>0){index=index2}
#       temp<-cbind(all[index,],rs[j,])
#       temp$INHERITANCE="-"
#       for(m in 1:dim(temp)[1]){
#         #    print(temp[m,])
#         if(temp$proband[m]==temp$V1[m]){temp$INHERITANCE[m]='UN'}
#         if(length(which(ped$ID==temp$proband[m]))>0){
#           if(temp$V1[m]==ped$Father[which(ped$ID==temp$proband[m])]){temp$INHERITANCE[m]='P-inherited'}
#           if(temp$V1[m]==ped$Mother[which(ped$ID==temp$proband[m])]){temp$INHERITANCE[m]='M-inherited'}
#         }
#         ref=temp$REF[m];
#         alt=temp$ALT[m];
#         
#         if(nchar(ref)>1 & nchar(alt)>1){
#           print(paste("1 ",ref,alt,sep="  "))
#           s=substr(ref,2,nchar(ref))
#           if(nchar(ref)>nchar(alt)){
#             s=substr(alt,2,nchar(alt))
#           }
#           ref=gsub(paste(s,"$",sep=""),"",ref)
#           alt=gsub(paste(s,"$",sep=""),"",alt)
#           print(paste("2  ",ref,alt,sep="  "))
#           temp$REF[m]=ref
#           temp$ALT[m]=alt
#         }
#       }
#       rem<-rbind(rem,temp)
#     }
#   }
#   write.table(unique(rem),file = "PAH/Result//Data/output/PAH-CHD.253KNG.TrueIGV.txt",row.names = F,sep="\t",quote=F)
#   
# }

cbioportal<-function(fwes,fsanger){
  wes<-read.csv(fwes,header = 1,check.names = F,comment.char = "",stringsAsFactors = F)
  sanger<-read.csv(fsanger,header=1,check.names = F,comment.char = "",stringsAsFactors = F)
  sanger$NucleotideChange<-unlist(lapply(sanger$Mutation,FUN = function(x) gsub(" ","",unlist(strsplit(x,split = ":"))[1])))
  sanger$ProteinChange<-unlist(lapply(sanger$Mutation,FUN = function(x) gsub(" ","",unlist(strsplit(x,split = ":"))[2])))
  nm<-names(wes)
  nm[which(nm=="Gene.refGene")]="Gene"
  nm[which(nm=="proband")]="ID"
  names(wes)<-nm
  type="child"
  child_cb<-rbind(wes[intersect(grep(type,wes$age,ignore.case = T),grep("PAH",wes$type,ignore.case = T)),c("Gene","ID","NucleotideChange","ProteinChange","Mut_Type")],
                  sanger[intersect(grep(type,sanger$Type,ignore.case = T),which(sanger$Group=="PAH")),c("Gene","ID","NucleotideChange","ProteinChange","Mut_Type")])
  
  
  type="adult"
  adult_cb<-rbind(wes[intersect(grep(type,wes$age,ignore.case = T),grep("PAH",wes$type,ignore.case = T)),c("Gene","ID","NucleotideChange","ProteinChange","Mut_Type")],
                  sanger[intersect(grep(type,sanger$Type,ignore.case = T),which(sanger$Group=="PAH")),c("Gene","ID","NucleotideChange","ProteinChange","Mut_Type")])
  
  adult_cb$ProteinChange[which(is.na(adult_cb$ProteinChange))]<-adult_cb$NucleotideChange[which(is.na(adult_cb$ProteinChange))]
  child_cb$ProteinChange[which(is.na(child_cb$ProteinChange))]<-child_cb$NucleotideChange[which(is.na(child_cb$ProteinChange))]
  adult_cb<-adult_cb[order(adult_cb$Gene),]
  child_cb<-child_cb[order(child_cb$Gene),]
  write.table(adult_cb[,c(1,2,4,5)],file = "PAH/Result/Data/output/Pedvsadult.cbioportal.adult.txt",quote = F,row.names = F,sep="\t")
  write.table(child_cb[,c(1,2,4,5)],file = "PAH/Result/Data/output/Pedvsadult.cbioportal.child.txt",quote = F,row.names = F,sep="\t")
  
}
fsanger="PAH/Result/Data/source/sanger_HPAH.csv"
fwes="PAH/Result/Data/source/PAH_ped_adult_All.known.csv"
cbioportal(fwes,fsanger )


start<-function(){
      indel_check  <-  read.table("PAH/JointCalls/IGV.check.indel.txt",header = 1,check.names = F,stringsAsFactors = F)
      indel_check  <-  indel_check[which(indel_check$Flag==0),]
      
      case_false_indels  <-  read.table("PAH/Image/Remote.false.list",header = F,stringsAsFactors = F,check.names = F)
      ctr_false_indels  <-  read.table("PAH/Image/control.false.indel.txt",header=F,stringsAsFactors = F,check.names = F)
      exclude_set_cs  <-  unlist(lapply(1:dim(case_false_indels)[1],FUN = function(x) paste(case_false_indels$V1[x],case_false_indels$V2[x],case_false_indels$V3[x]) ))
      exclude_set_ctr  <-  unlist(lapply(1:dim(ctr_false_indels)[1],FUN = function(x) paste(ctr_false_indels$V1[x],ctr_false_indels$V2[x],ctr_false_indels$V3[x]) ))
      
      
      out_prefixs<-c("PAH_CHD_rare_enrichment_253CHD","PAH_ped_adult_rare_enrichment")
      fdat<-"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed"
      fcohorts<-c("PAH/Result/Data/source/PAH-CHD-list.csv","PAH/Result/Data/source/PAH_pediatric_adult_list.csv")
      mutlists<-c("PAH/PAH-CHD/Data/PAH-CHD.solved.0608.csv" ,"PAH/Result/Data/source/PAH_ped_adult_All.known.csv")
      fassos  <- c("PAH/Result/Data/source/CHD.253GeneList.csv","PAH/Result/Data/source/PAH_associated11-13.txt") #"PAH/Result/Data/source/CHD.253kngenes.list"
      
      
      dat <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")
      #dat<-seperate_proband_info(dat)
      dat<-formatFreq_new(dat)
      ######## PAH_settting ########
      filter_dat <- filter_allfreq_PAH(dat,freq_avg = 0.001,0.001)  ### PAH_setting
      filter_dat<- filter_dat[which(as.numeric(filter_dat$AC)<25),] ## PAH_setting
      filter_dat <-  VQSR_filter(filter_dat) ## PAH_setting
      filter_dat<- COV_filter(filter_dat) ## PAH_setting
      filter_dat <- filter_dat[which(!paste(filter_dat$ProbandName,filter_dat$CHROM,filter_dat$POS)  %in% c( exclude_set_ctr,exclude_set_cs)),] ## PAH_setting
      filter_dat <- filter_dat[which(filter_dat$ALT!="*"),] ## PAH_setting
      filter_dat<- map_filter(filter_dat) ## PAH_setting
      filter_dat<- filter_dat[which(as.numeric(filter_dat$RD_ind)/as.numeric(filter_dat$DP_ind)<0.8),] ### PAH_setting
      
      
      fcohort<-fcohorts[2]
      mutlist<-mutlists[2]
      fasso<-fassos[2]
      out_prefix<-out_prefixs[2]
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      
      
      cohorts<-cohorts[intersect(which(cohorts$Disease=="FPAH"),grep("adult",cohorts$TYPE,ignore.case = T)),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"_FPAH_adult",sep=""),filter_dat,doset = 1)
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[intersect(which(cohorts$Disease=="FPAH" ),grep("child",cohorts$TYPE,ignore.case = T)),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"_FPAH_child",sep=""),filter_dat,doset = 1)
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[grep("adult",cohorts$TYPE,ignore.case = T),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"_adult",sep=""),filter_dat,doset = 1)
      
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[grep("child", cohorts$TYPE,ignore.case = T ),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"_child",sep=""),filter_dat,doset = 1)
      
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[intersect(which(cohorts$Disease=="IPAH"),grep("adult",cohorts$TYPE,ignore.case = T)),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"IPAH_adult",sep=""),filter_dat,doset = 1)
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[intersect(which(cohorts$Disease=="IPAH" ),grep("child", cohorts$TYPE,ignore.case = T )),]
      rare_enrich(cohorts,mutlist,fasso,paste(out_prefix,"IPAH_child",sep=""),filter_dat,doset = 1)
      
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      cohorts<-cohorts[intersect(which(cohorts$Disease=="IPAH" ),grep("child", cohorts$TYPE,ignore.case = T )),]
      ycohorts<-cohorts[which(cohorts$Age_FT<=5 & cohorts$Age_FT>-0.1 ),]
      rare_enrich(ycohorts,mutlist,fasso,paste(out_prefix,"LET5",sep=""),filter_dat,doset = 1)
      
      
      ocohorts<-cohorts[intersect(which(cohorts$Age_FT>5 ),grep("child",cohorts$TYPE,ignore.case = T)),]
      rare_enrich(ocohorts,mutlist,fasso,paste(out_prefix,"GT5",sep=""),filter_dat,doset = 1)
      
      ###########################
      
      # ############ PAH_CHD setting ############
      # 
       filter_dat <- filter_allfreq_CHD(dat,freq_avg = 0.0001,0.0001)  ### PAH-CHD _setting
       filter_dat<- filter_dat[which(as.numeric(filter_dat$AF)<0.01),] ## PAH-CHD _setting
       filter_dat <-  VQSR_filter(filter_dat) ## PAH-CHD_setting
       filter_dat<- COV_filter(filter_dat) ## PAH-CHD_setting
       filter_dat <- filter_dat[which(!paste(filter_dat$ProbandName,filter_dat$CHROM,filter_dat$POS)  %in% c( exclude_set_ctr,exclude_set_cs)),] ## PAH-CHD_setting
       filter_dat <- filter_dat[which(filter_dat$ALT!="*"),] ## PAH-CHD _setting
       filter_dat<- map_filter(filter_dat) ## PAH-CHD_setting
       filter_dat<- filter_dat[which(as.numeric(filter_dat$RD_ind)/as.numeric(filter_dat$DP_ind)<0.9),] ### PAH-CHD_setting
      # 
      # 
      # ########## PAH-CHD setting ###############
      # cohorts<-read.csv(fcohorts[1],header = 1,stringsAsFactors = F,strip.white = T)
      fcohort<-fcohorts[1]
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      rare_enrich(cohorts,mutlists[1],fassos[2],paste(out_prefixs[1],"_11PAH",sep=""),filter_dat,doset = 1)
      rare_enrich(cohorts,mutlists[1],fassos[1],paste(out_prefixs[1],"_253PAH",sep=""),filter_dat,doset = 1)
       

    
      
     

}
