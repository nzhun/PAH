#setwd("~/server/")
setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")
source("Pipeline/Enrichment/Denovo_enrichment.R")

filter_allfreq_EUR <- function(data,freq_avg,freq_max){
  
  data <- data[which(as.numeric(data$ExAC_NFE)< freq_max
                     
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     
  ),]
  
  return(data)
}

filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(
                     na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                     #    &as.numeric(data$gnomAD_genome_EAS)<freq_max
                      #    &as.numeric(data$gnomAD_genome_NFE)<freq_max
                     #    &as.numeric(data$gnomAD_genome_FIN)<freq_max
                     #   &as.numeric(data$gnomAD_genome_OTH)<freq_max
                     #  &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                     # &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     #&as.numeric(data$gnomAD_genome_ALL)<freq_max
                     &
                    as.numeric(data$gnomAD_exome_ALL)<freq_max
                     &as.numeric(data$gnomAD_exome_EAS)<freq_max
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     &as.numeric(data$gnomAD_exome_FIN)<freq_max
                     &as.numeric(data$gnomAD_exome_OTH)<freq_max
                     &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_exome_AMR)<freq_max
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  # data <- data[which( as.numeric(data$AC)< 25
  #                      &as.numeric(data$AB)>0.2
  #  ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  if(length(grep("genomicSuperDups",names(data)))>0){
    index<-grep("Score",data$genomicSuperDups)
    as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))>0.9)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
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

VQSR_filter <- function(data){
  
  # if(length(grep("VQSLOD",names(data)))>0){
  #   #data <- data[which(data$VQSLOD > -2.75),]
  #   data <- data[which(data$VQSLOD > -0.5),]
  # }
  
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
                       "VQSRTrancheSNP99.80to99.90",
                    #   "VQSRTrancheSNP99.70to99.80",
                       #"VQSRTrancheSNP99.60to99.70",#"VQSRTrancheSNP99.50to99.60",
                       #"VQSRTrancheINDEL99.50to99.60",
                       #"VQSRTrancheINDEL99.60to99.70",                 
                    #   "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", 
                       "VQSRTrancheINDEL99.90to100.00")
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  return(data)
}

map_filter<-function(data){
  index<-grep("mappability",names(data),ignore.case = T)
  if(length(index)>0){ 
    data <- data[which(as.numeric(data[,index])==1),]
  }
  if(length(grep("genomicSuperDups",names(data),ignore.case = T))>0){
    index<-grep("Score",data$genomicSuperDups)
    
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return (data)
}

qcase<-function(data){
  return(dim(unique(data[,c("ProbandName","Gene.refGene")]))[1])
}

test_geneset<-function(total_case,total_control,xGen_PAH,chr_vars,genesets,correct){
  
  types<-unique(c(xGen_PAH$ExonicFunc.refGene,chr_vars$ExonicFunc.refGene))
  
  lof_class<-types[grep("^stop|^frame",types)]  #c("stopgain","frameshiftinsertion","frameshiftdeletion","stoploss","frameshift substitution","frameshift_deletion","frameshift_insertion")
  lof_fun<-"splicing"
  binom<-c()
  rs<-c()
  subcase<-c()
  subcontrol<-c()
  #  print (gene)
  if(length(genesets)>1){
    index_case<-which(xGen_PAH$Gene.refGene%in%genesets)
    index_control<-which(chr_vars$Gene.refGene%in%genesets)
    
    if(length(index_case)<1){len_case=0}else{subcase<-xGen_PAH[index_case,] }
    if(length(index_control)<1){len_control=0;}else{  subcontrol<-chr_vars[index_control,]}
  }else{
    subcase<-xGen_PAH
    subcontrol<-chr_vars
  }
  ## synonymous
  #if(length(index_case)<1){next}
  
  
  keys<-c("ExonicFunc.refGene","VarClass")
  if(length(subcase)>0){
    len_case=qcase(subcase[grep("^synonymous",subcase$ExonicFunc.refGene),c("ProbandName","Gene.refGene")]);
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol$ExonicFunc.refGene)])
  }
  test_syn<-Btest(len_case,total_case,len_control,total_control,correct) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-c("SYN",len_case,len_control,test_syn$p.value,test_syn$estimate)

  ### LGD
  if(length(subcase)>0){
    len_case=qcase(subcase[which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun),c("ProbandName","Gene.refGene")]);
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  }
  test_lof<-Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("LGD",len_case,len_control,test_lof$p.value,test_lof$estimate))
  
  
  # mis
  if(length(subcase)>0){
    len_case=qcase(subcase[grep("^nonsynonymous",subcase$ExonicFunc.refGene),c("ProbandName","Gene.refGene")]);
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[grep("^nonsynonymous",subcontrol$ExonicFunc.refGene)])
  }
  test_mis<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
  rs<-rbind(rs,c("MIS",len_case,len_control,test_mis$p.value,test_mis$estimate))
  
  
  ### missense D (CADD>=25)
  if(length(subcase)>0){
    len_case=qcase(subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,  which(as.numeric(subcase$CADD_phred)>=25)),c("ProbandName","Gene.refGene")]);
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene), which(as.numeric(subcontrol$CADD13_PHRED)>=25))])
  }
  test_cadd<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
  rs<-rbind(rs,c("CADD25",len_case,len_control,test_cadd$p.value,test_cadd$estimate))
  
  ## missense D +LOF
  if(length(subcase)>0){
    len_case=qcase(subcase[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,
                                                           which(as.numeric(subcase$CADD13_PHRED)>=25)), which(
                                                             subcase$ExonicFunc.refGene %in% lof_class |
                                                               subcase$Func.refGene%in%lof_fun)
    ),c("ProbandName","Gene.refGene")]) ;
  }
  if(length(subcontrol)>0){
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) ,
                                                  which(as.numeric(subcontrol$CADD13_PHRED)>=25)),which( 
                                                    subcontrol$ExonicFunc.refGene%in%lof_class |
                                                      subcontrol$Func.refGene%in%lof_fun))])
  }
  test_cadd_lof<-Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+LGD",len_case,len_control,test_cadd_lof$p.value,test_cadd_lof$estimate))
  
  
  
  ### missense D (CADD>=25 and metasvm ==d D)
  if(length(subcase)>0){
    len_case=qcase(subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),which(as.numeric(subcase$CADD_phred)>=25  & subcase$MetaSVM_pred=="D" )),c("ProbandName","Gene.refGene")]);
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) , which(as.numeric(subcontrol$CADD13_PHRED)>=25 & subcontrol$MetaSVM_pred=="D"))])
  }
  test_cadd_meta<-Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+METASVM",len_case,len_control,test_cadd_meta$p.value,test_cadd_meta$estimate))
  
  ## missense D +LOF
  if(length(subcase)>0){
    len_case=qcase(subcase[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),
                                                           which((as.numeric(subcase$CADD13_PHRED)>=25  & subcase$MetaSVM_pred=="D"))) ,
                                                 which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)),c("ProbandName","Gene.refGene")])[1];
    
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$CADD13_PHRED)>=25   & subcontrol$MetaSVM_pred=="D"))
                                        ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
  }
  
  test_cadd_meta_lof<- Btest(len_case,total_case,len_control,total_control,correct)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+METASVM+LGD",len_case,len_control,test_cadd_meta_lof$p.value,test_cadd_meta_lof$estimate))
  
  
  
  
  ## missense D( mcap>=0.05)
  if(length(subcase)>0){
    index<-intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05))
    len_case<-qcase(subcase[index,]);
    #len_case=dim(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05)),c("ProbandName","Gene.refGene")]))[1];
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))])
  }
  test_mcap<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("mcap>=0.5",len_case,len_control,test_mcap$p.value,test_mcap$estimate))
  
  ## missense D +LOF
  #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$MCAP)>=0.05)
  #                                                  |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
  #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$MCAP)>=0.05)
  #                                       | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  
  if(length(subcase)>0){
    len_case=qcase((subcase[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                           which((as.numeric(subcase$MCAP)>=0.05  ))) ,
                                                 which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)),]) );
  }
  if(length(subcontrol)>0){
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))
                                        ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
  }
  test_mcap_lof<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("MCAP>=0.05+LGD",len_case,len_control,test_mcap_lof$p.value,test_mcap_lof$estimate))
  
  
  ## missense D( revel>=0.75)
  ## missense D +LOF
  if(length(subcase)>0){
    
    len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.5)),]));
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))])
  }
  #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
  #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
  test_revel<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL>=0.5",len_case,len_control,test_revel$p.value,test_revel$estimate))
  
  ## missense D +LOF
  #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
  #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  if(length(subcase)>0){
    len_case=qcase((subcase[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                           which((as.numeric(subcase$REVEL)>=0.5  ))) ,
                                                 which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)),]) );
  }
  if(length(subcontrol)>0){
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))
                                        ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
  }
  test_revel_lof<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL0.5+LGD",len_case,len_control,test_revel_lof$p.value,test_revel_lof$estimate))
  
  ## missense D(mcap>=0.75)
  #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75)]));
  #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
  if(length(subcase)>0){
    len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.75)),]));
  #  len_case=length(unique(subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.75)),c("ProbandName","Gene.refGene")]));
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))])
  }
  
  test_revel2<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL0.75",len_case,len_control,test_revel2$p.value,test_revel2$estimate))
  
  ## missense D +LOF
  # len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
  #  len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.75) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  if(length(subcase)>0){
    len_case=qcase((subcase[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                           which((as.numeric(subcase$REVEL)>=0.75  ))) ,
                                                 which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)),]) );
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))
                                        ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
  }
  test_revel2_lof<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL0.75+LGD",len_case,len_control,test_revel2_lof$p.value,test_revel2_lof$estimate))
  
  
  binom<-rs
  ## missense D +LOF  
  #}
  #print("done")
  binom<-data.frame(as.matrix(binom),stringsAsFactors = F,check.names = F,row.names = NULL)
  
  names(binom)<-c("Type",paste("N_case/",total_case,sep="" ),paste("N_control/",total_control,sep=""),"Pvalue","OR" ) #,
                #  "LOF_case","LOF_control","LOF_Pvalue","LOF_OR",
                #  "MISN_case","MISN_control","MIS_Pvalue","MIS_OR",
                #  "CaddN_case","CaddN_control","Cadd_Pvalue","Cadd_OR",
                #  "Cadd_LOF_N_case","Cadd_LOF_N_control","Cadd_LOF_Pvalue","Cadd_LOF_OR",
               #   "CaddN_Meta_case","CaddN_Meta_control","Cadd_Meta_Pvalue","Cadd_Meta_OR",
               #   "Cadd_Meta_LOF_N_case","Cadd_Meta_LOF_N_control","Cadd_Meta_LOF_Pvalue","Cadd_Meta_LOF_OR",
              #    "MCAP_N_case","MCAP_N_control","MCAP_Pvalue","MCAP_OR",
               #   "MCAP_LOF_N_case","MCAP_LOF_N_control","MCAP_LOF_Pvalue","MCAP_LOF_OR",
              #    "REVEL_N_case","REVEL_N_control","REVEL_Pvalue","REVEL_OR",
             #     "REVEL_LOF_N_case","REVEL_LOF_N_control","REVEL_LOF_Pvalue","REVEL_LOF_OR",
            #      "REVEL2_N_case","REVEL2_N_control","REVEL2_Pvalue","REVEL2_OR",
            #      "REVE2L_LOF_N_case","REVEL2_LOF_N_control","REVEL2_LOF_Pvalue","REVEL2_LOF_OR"
 # )
  return(binom)
}

process<-function(total_case,total_control,xGen_PAH,outname,chr_vars){
  #chr=8;
  print(outname)
  lib_exp<-load_dataset();
  exac<-load_exac()
  sets<-load_genesets(lib_exp,exac)
  genesets="ALL";
  group="ALL--before"
  All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,1)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
 # correct<-as.numeric(sbinom[which(sbinom$Type=="SYN"),"OR"])
  correct<-1
  print(correct)
  genesets="ALL";
  group="ALL"
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  ####### 25
  hhe_genesets=sets$gene_hhe;
  group="Heart_top25"
 # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhe_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hhe_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  ### channel or pathway 
  
  #### chanel ### 
  chanelgenes<-read.csv("PAH/documents/Channelopathy_gene_list.csv",stringsAsFactors = F)
  genesets<-chanelgenes$GeneName;
  group=paste(dim(chanelgenes)[1],"channel genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  #### chanel ### 
  chanelgenes<-read.table("PAH/documents/Channel_highlyExps.txt",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"channel hig expression genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  chanelgenes<-read.table("PAH/documents/All_channelpluspotassium_genes.txt",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"All_channelpluspotassium_genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  chanelgenes<-read.table("PAH/documents/All_ion_plus_channel.txt",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"All_ion_lus_channel")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  chanelgenes<-read.csv("PAH/documents/ion_channel_paecs.txt",stringsAsFactors = F,header=F)
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"ION_Channel_PAECs")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  
  chanelgenes<-read.table("PAH/documents/Potassium.channle.txt",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"Potassium")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  
  chanelgenes<-read.csv("PAH/documents/ION_Chanel.anno.csv",stringsAsFactors = F)
  genesets<-chanelgenes$GeneName;
  group=paste(dim(chanelgenes)[1],"ION_Channel")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  chanelgenes<-read.csv("PAH/documents/human_SMC_genelists.csv",stringsAsFactors = F)
  genesets<-chanelgenes$GeneName;
  group=paste(dim(chanelgenes)[1],"Human SMC genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  chanelgenes<-read.csv("PAH/documents/human_endothelial_genelists.csv",stringsAsFactors = F)
  genesets<-chanelgenes$GeneName;
  group=paste(dim(chanelgenes)[1],"Human endothelial genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  chanelgenes<-read.csv("PAH/documents/human_epithelial_genelists.csv",stringsAsFactors = F)
  genesets<-chanelgenes$GeneName;
  group=paste(dim(chanelgenes)[1],"Human epithelial genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  chanelgenes<-read.table("PAH/documents/BMP_signaling-HSA-201451].gene.list",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"BMP signaling genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  print("BMP")
  chanelgenes<-read.table("PAH/documents/PAH-UK-DMIS-LOF.gene.list",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"UK_DmiS_lof genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  print("UK")
  
  chanelgenes<-read.table("PAH/documents/TGF-Beta.pathway.genes.list",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"TGF-beta pathway genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  print("TGFB") 
  
  chanelgenes<-read.table("PAH/documents/VEGF_gene_set.txt",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"VEGF")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  print("VEGF") 
  
  chanelgenes<-read.table("PAH/documents/Wnt-signal.pathway.genes.list",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes$V1;
  group=paste(dim(chanelgenes)[1],"WNT signaling")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  print("WNT")
  
  
  hle_genesets=sets$gene_hle;
  group="lung_top25"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hle_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hle_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  
  
  he_genesets=sets$gene_he;
  group="lung&Heart_top25"
 # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,he_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(he_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  hhle_genesets=sets$gene_hhle;
  group="heart or lung_top25"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhle_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hhle_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)

  
  ########50####################
  
  hhe_50_genesets=sets$gene_hhe_50;
  group="Heart_top50"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhe_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hhe_50_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  hle_50_genesets=sets$gene_hle_50;
  group="lung_top50"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hle_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hle_50_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  
  
  he_50_genesets=sets$gene_he_50;
  group="lung&Heart_top50"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,he_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(he_50_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  hhle_50_genesets=sets$gene_hhle_50;
  group="heart or lung_top50"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhle_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(hhle_50_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  ##############################
  print(group)
  
  
  misz_genesets=sets$misz_sets;
  group="misZ >3"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(misz_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  print(group)
  pli_genesets=sets$pli_sets;
  group="pLi >0.9"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,pli_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(pli_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  genesets=sets$lofz_sets;
  group="lofZ >3"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  chdgenes<-read.csv("PAH/Result/Data/source/CHD.253GeneList.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-chdgenes$Gene;
  group="253CHD"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  print(group)
  pahgenes<-read.csv("PAH/Result/Data/source/PAH_associated11-13.txt",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-pahgenes$Gene;
  group="11 PAH genes"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  print(group)
  
sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]
group="SOX17_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

## SOX17 and PLI >0.9
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_pli_genesets<- intersect(SOX17_genesets,pli_genesets)
group="SOX17_pLi_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


## SOX17 and misZ>3
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_misz_genesets<- intersect(SOX17_genesets,misz_genesets)
group="SOX17_misz_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_misz_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_misz_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)
print(group)
## SOX17 and lung top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
group="SOX17_lung25_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

print(group)

## SOX17 and heart top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_heart_genesets<- intersect(SOX17_genesets,hhe_genesets)
group="SOX17_heart25_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_heart_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


SOX17_lung_heart_genesets<- intersect(SOX17_genesets,he_genesets)
group="SOX17_lung_heart25_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_heart_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


############### top 50 ###################

## SOX17 and lung top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_lung_50_genesets<- intersect(SOX17_genesets,hle_50_genesets)
group="SOX17_lung50_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

print(group)

## SOX17 and heart top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_heart_50_genesets<- intersect(SOX17_genesets,hhe_50_genesets)
group="SOX17_heart50_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_heart_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


SOX17_lung_heart_50_genesets<- intersect(SOX17_genesets,he_50_genesets)
group="SOX17_lung_heart50_target"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_heart_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

print(group)

####### Exclude SOX17 #####


#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]
SOX17_genesets<-SOX17_genesets[which(SOX17_genesets!="SOX17")]
group="SOX17_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

## SOX17 and PLI >0.9
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_pli_genesets<- intersect(SOX17_genesets,pli_genesets)
group="SOX17_pLi_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

print(group)
## SOX17 and misZ>3
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_misz_genesets<- intersect(SOX17_genesets,misz_genesets)
group="SOX17_misz_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_misz_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_misz_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)
print(group)
## SOX17 and lung top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
group="SOX17_lung25_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)

print(group)

## SOX17 and heart top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_heart_genesets<- intersect(SOX17_genesets,hhe_genesets)
group="SOX17_heart25_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_heart_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


SOX17_lung_heart_genesets<- intersect(SOX17_genesets,he_genesets)
group="SOX17_lung_heart25_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_genesets,correct)
sbinom$Genesets<-paste(length(SOX17_lung_heart_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


############### top 50 ###################

## SOX17 and lung top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_lung_50_genesets<- intersect(SOX17_genesets,hle_50_genesets)
group="SOX17_lung50_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)



## SOX17 and heart top25
#sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_heart_50_genesets<- intersect(SOX17_genesets,hhe_50_genesets)
group="SOX17_heart50_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_heart_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)


SOX17_lung_heart_50_genesets<- intersect(SOX17_genesets,he_50_genesets)
group="SOX17_lung_heart50_target-1"
#All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_50_genesets,correct)
#    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
sbinom$Genesets<-paste(length(SOX17_lung_heart_50_genesets),group)
All_genes_binom<-rbind(All_genes_binom,sbinom)



## before 
# 
# sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,1)
# #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
# sbinom$Genesets<-paste(length(genesets),group,"_raw")
# All_genes_binom<-rbind(All_genes_binom,sbinom)
print(group)

  write.csv(All_genes_binom,file = paste("PAH/Result/Data/output/",outname,".gnomeAD.genesets.binom.csv",sep=""),row.names = F)
  return (All_genes_binom)
}

start<-function(){
    #fdat<-"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed"
    fdat<-"PAH/Result/Data/source/PAH_07102017.filteredvariants.anno.map.xgen_nbg.0508.bed.gz"
    flists<-c(  "PAH/Result/Data/source/PAH_EUR.list","PAH/Result/Data/source/PAH_CHD_IPAH_FPAH.list.csv",
              "PAH/Result/Data/source/PAH-CHD-list.csv",
              "PAH/Result/Data/source/PAH_FPAH.list.csv",
              "PAH/Result/Data/source/PAH_IPAH.list.csv",
            
              "PAH/Result/Data/source/PAH_pediatric_adult_list.csv")
    fctr="PAH/Result/Data/source/gnomad.genomes.r2.0.2.sites.NFE.0.01.anno.mapp.xgen_vcr.bed.gz"
    fexclude<-"PAH/Result/Data/source/ExcludingCUMCFromAsso.txt"
    excluding<-read.csv(fexclude,stringsAsFactors = F,header = F)
    excluding<-excluding$V1
    chr_vars_all<-read.table(fctr,stringsAsFactors = F,check.names = F,sep="\t",header = 1,comment.char = "",quote = "")
    
    
   
    total_control=7509
    
    #process(total_case,total_control,PAH_CHD_dat,"PAH-CHD",chr_vars)
    #process(1321,total_control,pcgc_parent,"PCGC-Parent",chr_vars)
    
    #chd_dat <- read.csv("PAH/PAH-CHD/CHD/CHD.WES.253known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
    chr_vars<-chr_vars_all[which(chr_vars_all$FILTER=="PASS" & chr_vars_all$ExonicFunc.refGene!="unknown"),]
    chr_vars <- formatFreq_new(chr_vars)
    chr_vars <- filter_allfreq_local(chr_vars,0.0001,0.0001)  ## Max gnomad exome <10^-4
    chr_vars <- chr_vars[which(chr_vars$AF<0.01),] ## gnomad WGS cohort frequency < 10^-3
    chr_vars<- map_filter(chr_vars)
    chr_vars<-chr_vars[which(chr_vars$AN_NFE > 0.9*max(chr_vars$AN_NFE) &chr_vars$AN > 0.9*max(chr_vars$AN) ),]
    chr_vars<-chr_vars[which(chr_vars$DP/chr_vars$AN > 15),]
  #  table(chr_vars$ExonicFunc.refGene)
    chr_vars<- chr_vars[which(chr_vars$MQ >=40),]
    chr_vars<- chr_vars[which(as.numeric(chr_vars$VQSLOD) >-2.5 & as.numeric(chr_vars$QD)  > 2 ),]
    remove<-intersect(grep('non',chr_vars$ExonicFunc.refGene),which(chr_vars$FS > 20 | chr_vars$SOR >2 |as.numeric(chr_vars$QD) <2 |as.numeric(chr_vars$VQSLOD) < -2.5 |chr_vars$MQ < 40 ))
    if(length(remove) >0){chr_vars<-chr_vars[-remove,]}
    
    print(sum(chr_vars$AC_NFE[grep("^syn",chr_vars$ExonicFunc.refGene)])/7509)
    print(sum(chr_vars$AC_NFE[grep("^nonsyn",chr_vars$ExonicFunc.refGene)])/7509)
    print(sum(chr_vars$AC_NFE[grep("non",chr_vars$ExonicFunc.refGene,invert = T)])/7509)
    
    PAH_case_control <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")
    if(length(grep("ProbandName",names(PAH_case_control)))==0){
      PAH_case_control$ProbandName<-PAH_case_control$proband
    }
    pcgc_european<-read.table("PAH/Result/Data/source/european-control.txt",stringsAsFactors = F,check.names = F,header =F )
    
    pcgc_european<-pcgc_european$V1
    
    pop <- read.csv("PAH/Result/Data/source/Phenotype_FREEZE6_pop.csv",header=1,stringsAsFactors = F,check.names = F)
    
    pop_vcr  <-  read.table("PAH/Result/Data/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
    chd.list<-read.csv(flists[1],header = 1,stringsAsFactors = F,strip.white = T)
    if(dim(chd.list)[2]<2){
      chd.list<-read.csv(flists[1],header = 1,stringsAsFactors = F,strip.white = T)
    }
    chd.list$FN<-chd.list$proband
    
    chd.list<-chd.list[which(!chd.list$FN%in%excluding),]
    
    vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]
    xgen_CHD_european <- intersect(chd.list$proband,pop$sampleID[which(pop$pop=="EUR")])
    case_IDs<-c(vcr_chd_european,xgen_CHD_european)
    total_case=length(xgen_CHD_european)+length(vcr_chd_european)
  
    
    filter_dat<- PAH_case_control[which(PAH_case_control$ProbandName%in%  case_IDs &PAH_case_control$ExonicFunc.refGene!="unknown"),]
 
    filter_dat <- formatFreq_new(filter_dat)
    filter_dat <- filter_allfreq_local(filter_dat,0.0001,0.0001) ## 
    #filter_dat<-COV_filter(filter_dat)
    filter_dat<-map_filter(filter_dat)  ## mappability ==1
    filter_dat<- filter_dat[which(filter_dat$AF<0.01),] ## cohort frequen < 10-3 # carefule when AN<10000
    filter_dat<-filter_dat[which(filter_dat$AN/max(filter_dat$AN) > 0.75),]
    filter_dat<- filter_dat[which(as.numeric(filter_dat$VQSLOD) > -5 ),]
    filter_dat<- filter_dat[which(filter_dat$MQ > 20 ),]
    filter_dat<- filter_dat[which(filter_dat$DP_ind >=7),]
    filter_dat<- filter_dat[which(filter_dat$AD_ind >=3),]
    filter_dat<- filter_dat[which(filter_dat$GQ >=30),]
   # filter_dat<- filter_dat[which(as.numeric(filter_dat$QD) > 2),]
    filter_dat<-filter_dat[which(filter_dat$GnomAD_Genome_cov10 > 0.9),]
    remove<-which(as.numeric(filter_dat$ReadPosRankSum) < -4 )
    if(length(remove)>0){filter_dat<-filter_dat[-remove,]}
    filter_dat<-filter_dat[grep('[3-9]',filter_dat$GT,invert = T),]
    filter_dat<- filter_dat[which((filter_dat$AB) >=0.2),]
   # remove<-intersect(grep("non",filter_dat$ExonicFunc.refGene,invert = T),which(filter_dat$AB <0.2 | filter_dat$AD_ind < 4| filter_dat$DP_ind < 7))
    #if(length(remove) >0){filter_dat<-filter_dat[-remove,]}
    dim(filter_dat[grep("^synony",filter_dat$ExonicFunc.refGene),])[1]/total_case
    dim(filter_dat[grep("^nonsynony",filter_dat$ExonicFunc.refGene),])[1]/total_case
    dim(filter_dat[grep("non",filter_dat$ExonicFunc.refGene,invert = T),])[1]/total_case
   # filter_dat<-VQSR_filter(filter_dat) ### VQSR >-2.75
    filter_dat$Is_homo<-0
    filter_dat$Is_homo[grep("1/1",filter_dat$GT)]<-1
    #PAH_case_control <- read.csv("PAH/Joint_calls_20170715/VCX_Control.inherited.csv_filtered_revel.csv",header = 1,sep = ",",stringsAsFactors = F,check.names = F)
    chr_vars$proband <- "."
    chr_vars$ProbandName <- "."
    filter_dat$AC_NFE <- 1
    filter_dat$Hom_NFE <- 1
    chr_vars$Is_homo <- F
    chr_vars$Genotype <- "."
    nms<-intersect(names(chr_vars),names(filter_dat))
    rnm<-c("AC_NFE","CADD_phred","ExonicFunc.refGene", "Gene.refGene","Hom_NFE", "proband",    "REVEL");
   # rnm[which(!rnm%in%nms)]
   # write.table(rbind(filter_dat[,nms],chr_vars[,nms]),row.names = F,quote =F,file = "PAH/Result/CUMC_gnomad/merged.CUMC.gnomad.filtered.bed")
    for(flist in flists){
      
        #"PAH/PAH-CHD/Data/PAH-CHD-list.csv"
        print(flist)
        chd.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
        if(dim(chd.list)[2]<2){
          chd.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
        }
        chd.list$FN<-chd.list$proband
        key=unlist(strsplit(basename(flist),split=".list"))[1]
        vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]
        
        
        xgen_CHD_european <- intersect(chd.list$proband,pop$sampleID[which(pop$pop=="EUR")])
        
        write.table(c(xgen_CHD_european,vcr_chd_european),file = paste("PAH/Result/CUMC_gnomad/VAT/REVEL/Bionomial/",key,"_EUR.list",sep=""),quote = F)
        total_case=length(xgen_CHD_european)+length(vcr_chd_european)
        PAH_CHD_dat<-filter_dat[which(filter_dat$ProbandName%in%c(xgen_CHD_european,vcr_chd_european)),]
        #pcgc_parent<-filter_dat[which(filter_dat$ProbandName%in%c(pcgc_european)),]
        #pah_chd<-read.table("PAH/PAH-CHD/PAH-CHD.0.0001.european.REVEL.mappablity.bed",comment.char = "",header = 1,stringsAsFactors = F,check.names = F)
    
        auto_gnomad<-chr_vars[which(chr_vars$`#CHROM`!="X"),]
        auto_pah<-PAH_CHD_dat[which(PAH_CHD_dat$`#CHROM`!="X" ),]
        #names(chd_dat)[grep("proband",names(chd_dat))]<-"ProbandName"
     
        result<-process(total_case,total_control,PAH_CHD_dat,key,chr_vars) #,"ALL","ALL")
     #   result<-process(total_case,total_control,auto_pah,paste(key,"_Auto",sep=""),auto_gnomad) #,"ALL","ALL")
  }
}
start()
# 
# 
# qqq_binom(len_case,len_control){
#   len_case=143
#   len_control=7509
#   data<-read.csv("~/server/PAH/PAH-CHD/PAH-CHD.gnomeAD.ALL.binom.csv",header = 1,check.names = F,stringsAsFactors = T,comment.char = "")
#   pvalues<-data$REVEL_LOF_Pvalue[which(data$REVEL_LOF_Pvalue>0)]
#   plot(-log10(pvalues),-log10(rbinom(n=length(pvalues),size=1,prob=len_case/(len_case+len_control))))
#   
# }