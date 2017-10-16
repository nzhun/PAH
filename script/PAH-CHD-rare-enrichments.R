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

filter_allfreq <- function(data){
  freq <- 0.0001
  freq2 <- 0.0001
  data <- data[which(na.pass(as.numeric(data$ExACfreq)< freq)
                     &na.pass(as.numeric(data$ExAC.amr.freq)< freq2)
                     &as.numeric(data$ExAC.afr.freq)< freq2
                     &as.numeric(data$ExAC.nfe.freq)< freq2
                     &as.numeric(data$ExAC.sas.freq)< freq2
                     &as.numeric(data$ExAC.eas.freq)< freq2
                   #  & as.numeric(data$`1KGfreq`) < freq2
                   #  &as.numeric(data$ESPfreq)< freq2
                   #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  data <- data[which( as.numeric(data$AC)< 25
                      &as.numeric(data$AB)>0.2
  ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
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

map_filter<-function(data){
  
  if(length(grep("mappability",names(data)))>0){ 
    data <- data[which(data$mappability==1),]
  }
  if(length(grep("genomicSuperDups",names(data)))>0){
    index<-grep("Score",data$genomicSuperDups)
  
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return (data)
}

VQSR_filter <- function(data){
  # if(length(grep("FILTER",names(data)))>0){
  #   filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
  #                      #"VQSRTrancheSNP99.70to99.80","VQSRTrancheSNP99.80to99.90",
  #                      #"VQSRTrancheSNP99.60to99.70","VQSRTrancheSNP99.50to99.60","VQSRTrancheINDEL99.50to99.60","VQSRTrancheINDEL99.60to99.70",                 
  #                      #"VQSRTrancheINDEL99.70to99.80",  
  #                      "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
  #   
  #   data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  # }
  if(length(grep("VQSLOD",names(data)))>0){
       data <- data[which(data$VQSLOD > -2.75),]
  }
  return(data)
}

COV_filter <- function(data){
  if(length(grep("GnomAD_Genome_cov",names(data)))>0){
    data<-data[which(as.numeric(data$GnomAD_Genome_cov10)>0.85),]
  }
  if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
    data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.8 & 2*data$AC_xGen/data$AN_xGen >0.8 & 2*(data$AC_VCR/data$AN_VCR) >0.8  &data$GnomAD_Genome_cov10>0.8 ),]
  }
  if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
    data <- data[which(2*data$AC_Control/data$AN_Control >0.8 & 2*data$AC_VU/data$AN_VU >0.8 ),]
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

#pah.list<-read.csv("PAH/Result/Data/PAH_pediatric_adult_list.csv",header = 1,stringsAsFactors = F,strip.white = T,check.names = F)
#pah_solved_vars<-read.csv("PAH/Result/Data/PAH_ped_adult_All.known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")



rare_enrich<-function(fdat,fcohort,mutlist,fasso,out_prefix,data){
      cohorts<-read.csv(fcohort,header = 1,stringsAsFactors = F,strip.white = T)
      known_vars<-read.csv(mutlist,header = 1,stringsAsFactors = F,check.names = F,comment.char = "")

      gene_asso  <-  read.table(fasso,header = F,stringsAsFactors = F,check.names = F)
      gene_asso  <-  gene_asso$V1
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
      indel_check  <-  read.table("PAH/Result/Data/source/IGV.check.indel.txt",header = 1,check.names = F,stringsAsFactors = F)
      indel_check  <-  indel_check[which(indel_check$Flag==0),]
      
      case_false_indels  <-  read.table("PAH/Result/Data/source/Remote.false.list",header = F,stringsAsFactors = F,check.names = F)
      ctr_false_indels  <-  read.table("PAH/Result/Data/source/control.false.indel.txt",header=F,stringsAsFactors = F,check.names = F)
      exclude_set_cs  <-  unlist(lapply(1:dim(case_false_indels)[1],FUN = function(x) paste(case_false_indels$V1[x],case_false_indels$V2[x],case_false_indels$V3[x]) ))
      exclude_set_ctr  <-  unlist(lapply(1:dim(ctr_false_indels)[1],FUN = function(x) paste(ctr_false_indels$V1[x],ctr_false_indels$V2[x],ctr_false_indels$V3[x]) ))
      
      
      
      ### start import variant #####
     
      
      control_dat <- data[which(data$ProbandName %in%  control_IDs),]
      xgen_dat<-data[which(data$ProbandName%in% xgen_european),]
      vcr_dat<-data[which(data$ProbandName%in% vcr_european),]
      outfile=paste("PAH/Result/Data/output/",out_prefix,".csv",sep="")
      enrich_test(control_dat,xgen_dat,vcr_dat,gene_asso,minor_gene_asso,solved.samples,outfile)
      
      outfile=paste("PAH/Result/Data/output/",out_prefix,"_child.csv",sep="")
      child_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
      child_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
      enrich_test(control_dat,child_xgen_dat,child_vcr_dat,gene_asso,minor_gene_asso,solved.samples,outfile)
      
      outfile=paste("PAH/Result/Data/output/",out_prefix,"_child.csv",sep="")
      child_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
      child_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("child",cohorts$TYPE)]),]
      enrich_test(control_dat,child_xgen_dat,child_vcr_dat,gene_asso,minor_gene_asso,solved.samples,outfile)
      
      for( dd in unique(cohorts$Disease)){
        outfile=paste("PAH/Result/Data/output/",dd,"_",out_prefix,"_child.csv",sep="")
        child_xgen_dat2<-child_xgen_dat[which(child_xgen_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
        child_vcr_dat2<-child_vcr_dat[which(child_vcr_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
        enrich_test(control_dat,child_xgen_dat2,child_vcr_dat2,gene_asso,minor_gene_asso,solved.samples,outfile)
      }
      
      outfile=paste("PAH/Result/Data/output/",out_prefix,"_adult.csv",sep="")
      adult_xgen_dat<-xgen_dat[which(xgen_dat$ProbandName%in%cohorts$proband[grep("adult",cohorts$TYPE)]),]
      adult_vcr_dat<-vcr_dat[which(vcr_dat$ProbandName%in%cohorts$proband[grep("adult",cohorts$TYPE)]),]
      enrich_test(control_dat,adult_xgen_dat,adult_vcr_dat,gene_asso,minor_gene_asso,solved.samples,outfile)
      
      for( dd in unique(cohorts$Disease)){
        outfile=paste("PAH/Result/Data/output/",dd,"_",out_prefix,"_child.csv",sep="")
        adult_xgen_dat2<-adult_xgen_dat[which(adult_xgen_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
        adult_vcr_dat2<-adult_vcr_dat[which(adult_vcr_dat$ProbandName%in%cohorts$proband[which(cohorts$Disease==dd)]),]
        enrich_test(control_dat,adult_xgen_dat2,adult_vcr_dat2,gene_asso,minor_gene_asso,solved.samples,outfile)
      }
      
}


enrich_test<-function(control_data,xgen_data,vcr_data,gene_asso,minor_gene_asso,solved.samples,outfile){
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
    df<-NA
    merged_known_test <-NA
    merged_minor_known_test<-NA
    df_unsolved_xgen<-NA
    df_unsolved_VCR <-NA
    merge_test_unsolved_all<-NA
    
    
    if(dim(xgen_data)[1]>0){
        init_all <- enrich_func(xgen_data,control_data,"All",F,c(),1,1)
        init_all$type="General_beforeCorrecting"
        init_all$Test="ALL"
        correction <- init_all$enrichment[which(init_all$class=="SYN")]
        correction2 <- length(get_indel(xgen_data))/length(unique(xgen_data$ProbandName))/(length(get_indel(control_data))/length(unique(control_data$ProbandName)))
      
        init <- enrich_func(xgen_data,control_data,"All",F,c(),correction,correction2)
        init$type <- "General"
        init$Test="ALL"
        enrich_asso <- enrich_func(xgen_data,control_data,gene_asso,F,c(),correction,correction2)
        enrich_asso$type = "XGEN_knownGene"
        enrich_asso$Test="XGEN_KnownGene"
        
        enrich_minor_asso <- enrich_func(xgen_data,control_data,minor_gene_asso,F,c(),correction,correction2)
        enrich_minor_asso$type="XGEN_minor_knownGene"
        enrich_minor_asso$Test="XGEN_minor_KnownGene"
        df <- enrichSet(xgen_data,control_data,gene_asso,correction,correction2,lib)  ## gene set 
        df$Test="all-xGen"
        
        df_unsolved_xgen <- enrichSet(xgen_data[which(!xgen_data$ProbandName %in% conf_solved_samples),],control_data,gene_asso,correction,correction2,lib)
        df_unsolved_xgen$Test="all-xGen-unsolved"
        
      }
      if(dim(vcr_data)[1]>0){
          enrich_asso_VCR <- enrich_func(vcr_data,control_data,gene_asso,F,c(),correction,correction2)
          enrich_asso_VCR$type="VCRknownGene"
          enrich_asso_VCR$Test="VCRKnownGene"
          
          enrich_minor_asso_VCR <- enrich_func(vcr_data[which(!vcr_data$ProbandName %in% solved.samples$ID[solved.samples$GeneName %in% c("BMPR2","TBX4")]),],control_data,minor_gene_asso,F,c(),correction,correction2)
          enrich_minor_asso_VCR$type="VCR_minor_knownGene"
          enrich_minor_asso_VCR$Test="VCR_minor_knownGene"
          
          df_vcr <- enrichSet(vcr_data,control_data,gene_asso,correction,correction2,lib)  ## gene set 
          df_vcr$Test="all-VCR"
          
          df_unsolved_VCR <- enrichSet(vcr_data[which(!vcr_data$ProbandName %in%  conf_solved_samples),],control_data,gene_asso,1,1,lib)
          df_unsolved_VCR$Test <- "all-VCR-unsolved"
      }
    
    
      if(dim(vcr_data)[1]>0 && dim(xgen_data)[1]>0){
    
        merged_known_test <- Merge_test_single(enrich_asso_VCR,enrich_asso,"merged_known_gene")
      
        merged_minor_known_test <- Merge_test_single(enrich_minor_asso_VCR,enrich_minor_asso,"merged_minor_known_gene")
        
        ### correction_merged
        merge_test_unsolved_all <- Merge_test_group(df_unsolved_VCR,df_unsolved_xgen,"unsolved_merged_all")
        
      }  
     
      write.csv(rbind(init_all,init,
                      enrich_asso,
                      enrich_minor_asso,
                      df,
                      merged_known_test,
                      merged_minor_known_test,
                      df_unsolved_xgen,df_unsolved_VCR,merge_test_unsolved_all
      ),outfile,row.names = F)
      
      print(paste("output to ",outfile,sep=""))
      
}


out_prefixs<-c("PAH_CHD_rare_enrichment","PAH_ped_adult_rare_enrichment")
fdat<-"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed"
fcohorts<-c("PAH/Result/Data/source/PAH-CHD-list.csv","PAH/Result/Data/source/PAH_pediatric_adult_list.csv")
mutlists<-c("PAH/PAH-CHD/Data/PAH-CHD.solved.0608.csv" ,"PAH/Result/Data/source/PAH_ped_adult_All.known.csv")
fasso  <-  "PAH/Result/Data/source/PAH_associated11-13.txt"


dat <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")
#dat<-seperate_proband_info(dat)
dat<-formatFreq_new(dat)

filter_dat <- filter_allfreq_new(dat,freq_avg = 0.0001,0.0001)
filter_dat<-filter_dat[which(filter_dat$AF<0.01),]
filter_dat <-  VQSR_filter(filter_dat)
filter_dat<-COV_filter(filter_dat)
filter_dat <- filter_dat[which(!paste(filter_dat$ProbandName,filter_dat$CHROM,filter_dat$POS)  %in% c( exclude_set_ctr,exclude_set_cs)),]
filter_dat <- filter_dat[which(filter_dat$ALT!="*"),]
filter_dat<-map_filter(filter_dat)
filter_dat<-filter_dat[which(filter_dat$AD_ind/filter_dat$DP_ind>0.2),]



for(i in 1:2){
  out_prefix<-out_prefixs[i]
  fcohort<-fcohorts[i]
  mutlist<-mutlists[i]
  rare_enrich(fdat,fcohort,mutlist,fasso,out_prefix,filter_dat)
}
