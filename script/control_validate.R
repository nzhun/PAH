## check control

#setwd("~/server/")
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

filter_allfreq <- function(data){
  freq <- 0.001
  freq2 <- 0.001
  data <- data[which(na.pass(as.numeric(data$ExACfreq)< freq)
                     &na.pass(as.numeric(data$ExAC.amr.freq)< freq2)
                     &as.numeric(data$ExAC.afr.freq)< freq2
                     &as.numeric(data$ExAC.nfe.freq)< freq2
                     &as.numeric(data$ExAC.sas.freq)< freq2
                     &as.numeric(data$ExAC.eas.freq)< freq2
                     & as.numeric(data$`1KGfreq`) < freq2
                     &as.numeric(data$ESPfreq)< freq2
                     &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  data <- data[which( as.numeric(data$AC)< 25
                      &as.numeric(data$AB)>0.2
  ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  return(data)
}

VQSR_filter <- function(data){
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.80to99.90","VQSRTrancheSNP99.90to100.00","VQSRTrancheSNP99.70to99.80",
                       #"VQSRTrancheSNP99.60to99.70","VQSRTrancheSNP99.50to99.60","VQSRTrancheINDEL99.50to99.60","VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  return(data)
}

COV_filter <- function(data){
  if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
    data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.9 & 2*data$AC_xGen/data$AN_xGen >0.9 & 2*(data$AC_VCR/data$AN_VCR) >0.9  &data$GnomAD_Genome_cov10>0.9 ),]
  }
  if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
    data <- data[which(2*data$AC_Control/data$AN_Control >0.9 & 2*data$AC_VU/data$AN_VU >0.9 ),]
  }
  return(data)
}


Merge_test_single <- function(dat_vcr,dat_xgen, TestName){
  merge_test <- c()
  
  c_capture1  <-  dat_vcr
  c_capture2  <-  dat_xgen
  for (class in unique(c_capture1$class)){
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
fasso  <-  "PAH/DOC/HongJian/source/PAH_associated11-13.txt"
gene_asso  <-  read.table(fasso,header = F,stringsAsFactors = F,check.names = F)
gene_asso  <-  gene_asso$V1
minor_gene_asso  <-  setdiff(gene_asso,c("BMPR2","TBX4"))
####

control_IDs  <-  read.csv("PAH/DOC/HongJian/source/european-control.txt",header = F,check.names = F,stringsAsFactors = F)
control_IDs   <-   control_IDs$V1


### indel IGV check xGen ###  the indels failed IGV check
indel_check  <-  read.table("PAH/JointCalls/IGV.check.indel.txt",header = 1,check.names = F,stringsAsFactors = F)
indel_check  <-  indel_check[which(indel_check$Flag==0),]

case_false_indels  <-  read.table("PAH/Image/Remote.false.list",header = F,stringsAsFactors = F,check.names = F)
ctr_false_indels  <-  read.table("PAH/Image/control.false.indel.txt",header=F,stringsAsFactors = F,check.names = F)
exclude_set_cs  <-  unlist(lapply(1:dim(case_false_indels)[1],FUN = function(x) paste(case_false_indels$V1[x],case_false_indels$V2[x],case_false_indels$V3[x]) ))
exclude_set_ctr  <-  unlist(lapply(1:dim(ctr_false_indels)[1],FUN = function(x) paste(ctr_false_indels$V1[x],ctr_false_indels$V2[x],ctr_false_indels$V3[x]) ))



PAH_case_control <- read.table("PAH/JointCalls/PAHCASE_CONTROl_ALL.inherited.csv.european.slim.csv",header = 1,sep = "\t",stringsAsFactors = F,check.names = F)

control_dat <- PAH_case_control[which(PAH_case_control$ProbandName %in%  control_IDs),]


control_dat <-  formatFreq(control_dat)
### start filter the variants ###

control_dat <- filter_allfreq(control_dat)

control_dat <-  VQSR_filter(control_dat)

control_dat <- COV_filter(control_dat)


control_dat <- control_dat[which(!paste(control_dat$ProbandName,control_dat$CHROM,control_dat$POS)  %in%  exclude_set_ctr),]




