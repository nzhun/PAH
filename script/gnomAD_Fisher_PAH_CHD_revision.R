#setwd("~/server/")
setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")
source("Pipeline/Enrichment/Denovo_enrichment.R")

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
                                                  which(as.numeric(subcontrol$CADD13_PHRED)>=20)),which( 
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
  
  names(binom)<-c("Type","N_case","N_control","Pvalue","OR" ) #,
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

filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                  #   &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                  #   &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     & as.numeric(data$`1000g2015aug_all`) < 0.001
                  & as.numeric(data$`1000g2015aug_eur`) < 0.001
                  #   &as.numeric(data$ExAC_FIN)< freq_max
                  #   &as.numeric(data$ExAC_SAS)< freq_max
                  #   &as.numeric(data$ExAC_EAS)< freq_max
                  #   &as.numeric(data$ExAC_OTH)< freq_max
                     #    &as.numeric(data$gnomAD_genome_EAS)<freq_max
                     #     &as.numeric(data$gnomAD_genome_NFE)<freq_max
                     #    &as.numeric(data$gnomAD_genome_FIN)<freq_max
                     #   &as.numeric(data$gnomAD_genome_OTH)<freq_max
                     #  &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                     # &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     #&as.numeric(data$gnomAD_genome_ALL)<freq_max
                     &as.numeric(data$gnomAD_exome_ALL)<freq_max
                  #   &as.numeric(data$gnomAD_exome_EAS)<freq_max
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                  #   &as.numeric(data$gnomAD_exome_FIN)<freq_max
                  #   &as.numeric(data$gnomAD_exome_OTH)<freq_max
                  #   &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                  #   &as.numeric(data$gnomAD_exome_AMR)<freq_max
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
  # if(length(grep("genomicSuperDups",names(data)))>0){
  #   index<-grep("Score",data$genomicSuperDups)
  #   as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))
  #   dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))>0.9)]
  #   if(length(dup_indexs)>0){
  #     data <- data[-dup_indexs,]
  #   }
  # }
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


enrich_geneset <- function(total_case,total_control,xGen_PAH,outname,chr_vars){
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
  
  
  
  misz_genesets=sets$misz_sets;
  group="misZ >3"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(misz_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
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
  
  
  pahgenes<-read.csv("PAH/Result/Data/source/PAH_associated11-13.txt",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-pahgenes$Gene;
  group="11 PAH genes"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  
  sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]
#  sox17_target<-read.table("PAH/Result/Data/SOX17_enrichr.target.txt",header=F,check.names = F,stringsAsFactors = F)
#  SOX17_genesets<- sox17_target$V1
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
  
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
  group="SOX17_lung25_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  
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
  
  
  ## SOX17 and misZ>3
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_misz_genesets<- intersect(SOX17_genesets,misz_genesets)
  group="SOX17_misz_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(SOX17_misz_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
  group="SOX17_lung25_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  
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
  
  
  write.csv(All_genes_binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.genesets.binom.csv",sep=""),row.names = F)
  return (All_genes_binom)
}
VQSR_filter <- function(data){
  
  #  if(length(grep("VQSLOD",names(data)))>0){
  #    data <- data[which(data$VQSLOD > -0.5),]
  #  }
  
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
                       "VQSRTrancheSNP99.80to99.90",
                     #  "VQSRTrancheSNP99.70to99.80",
                       #"VQSRTrancheSNP99.60to99.70",#"VQSRTrancheSNP99.50to99.60",
                       #"VQSRTrancheINDEL99.50to99.60",
                     #  "VQSRTrancheINDEL99.60to99.70",                 
                   #    "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
    
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

process<-function(total_case,total_control,xGen_PAH,outname,chr_vars){
  #chr=8;
  print(outname)
  All_genes_binom<-c()
  #header<-read.table("PAH/PAH-CHD/gnomAD/gnomAD_filter/gnomad.header.txt",comment.char = "",header=1,stringsAsFactors = F,check.names = F)
  # for(chr in c(seq(1,22),"X")){
  #for(chr in c("X")){
  # print(chr)
  #  chr_vars<-read.table(paste("PAH/PAH-CHD/gnomAD/gnomAD_filter/gnomad.genomes.r2.0.1.sites.",chr,".corrected.bed.hg19_multianno.txt",sep=""),skip = 3,stringsAsFactors = F,check.names = F,sep="\t",header = F)
  # xGen_PAH<-pah_chd
  #names(chr_vars)<-names(header)
  
  types<-unique(c(xGen_PAH$ExonicFunc.refGene,chr_vars$ExonicFunc.refGene))
  
  lof_class<-types[grep("^stop|^frame",types)]  #c("stopgain","frameshiftinsertion","frameshiftdeletion","stoploss","frameshift substitution","frameshift_deletion","frameshift_insertion")
  lof_fun<-"splicing"
  binom<-c()
  for(gene in unique(c(chr_vars$Gene.refGene,xGen_PAH$Gene.refGene))){ #[which(xGen_PAH$`#CHROM`==chr)]
    rs<-c()
    subcase<-c()
    subcontrol<-c()
    nms<-c()
    #  print (gene)
    index_case<-which(xGen_PAH$Gene.refGene==gene)
    index_control<-which(chr_vars$Gene.refGene==gene)
    if(length(index_case)<1){len_case=0}else{subcase<-xGen_PAH[index_case,] }
    if(length(index_control)<1){len_control=0;}else{  subcontrol<-chr_vars[index_control,]}
    ## synonymous
    #if(length(index_case)<1){next}
    
    len_case=0
    len_control=0
    
    keys<-c("ExonicFunc.refGene","VarClass")
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[grep("^synonymous",subcase$ExonicFunc.refGene)]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol$ExonicFunc.refGene)])
    }
    test_syn<-Btest(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(gene,len_case,len_control,test_syn$p.value,test_syn$estimate)
    nms<- c("Gene","syN_case","syN_control","syn_Pvalue","syn_OR")
    
    ### LGD
    
    len_case=0
    len_control=0
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    }
    test_lof<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_lof$p.value,test_lof$estimate)
    nms<-c(nms,
           "LOF_case","LOF_control","LOF_Pvalue","LOF_OR")
    
    ### missense D (CADD>=25)
    for(cutoff in seq(15,35,5)){
          len_case=0
          len_control=0
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,  which(as.numeric(subcase$CADD_phred)>=cutoff))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene), which(as.numeric(subcontrol$CADD13_PHRED)>=cutoff))])
          }
          test_cadd<- Btest(len_case,total_case,len_control,total_control)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
          rs<-c(rs,len_case,len_control,test_cadd$p.value,test_cadd$estimate)
          nms<-c(nms,
                 paste("CADD",cutoff,"_N_case",sep=""),paste("CADD",cutoff,"_N_control",sep=""),paste("CADD",cutoff,"_Pvalue",sep=""),paste("CADD",cutoff,"_OR"))
          ## missense D +LOF
          len_case=0
          len_control=0
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,
                                                                   which(as.numeric(subcase$CADD13_PHRED)>=cutoff)), which(
                                                                     subcase$ExonicFunc.refGene %in% lof_class |
                                                                       subcase$Func.refGene%in%lof_fun)
            )])) ;
          }
          if(length(subcontrol)>0){
            
            len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) ,
                                                          which(as.numeric(subcontrol$CADD13_PHRED)>=cutoff)),which( 
                                                            subcontrol$ExonicFunc.refGene%in%lof_class |
                                                              subcontrol$Func.refGene%in%lof_fun))])
          }
          test_cadd_lof<-Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_cadd_lof$p.value,test_cadd_lof$estimate)
          
          nms<-c(nms,
                 "Cadd_LOF_N_case","Cadd_LOF_N_control","Cadd_LOF_Pvalue","Cadd_LOF_OR")
  }
  
    ### missense D (CADD>=25 and metasvm ==d D)
    
    len_case=0
    len_control=0
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),which(as.numeric(subcase$CADD_phred)>=20  & subcase$MetaSVM_pred=="D" ))]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) , which(as.numeric(subcontrol$CADD13_PHRED)>=20 & subcontrol$MetaSVM_pred=="D"))])
    }
    test_cadd_meta<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_cadd_meta$p.value,test_cadd_meta$estimate)
    nms<-c(nms,
           "CaddN_Meta_case","CaddN_Meta_control","Cadd_Meta_Pvalue","Cadd_Meta_OR")
    ## missense D +LOF
    
    len_case=0
    len_control=0
    
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),
                                                             which((as.numeric(subcase$CADD13_PHRED)>=20  & subcase$MetaSVM_pred=="D"))) ,
                                                   which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
      
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$CADD13_PHRED)>=20   & subcontrol$MetaSVM_pred=="D"))
                                          ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    }
    
    test_cadd_meta_lof<- Btest(len_case,total_case,len_control,total_control)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_cadd_meta_lof$p.value,test_cadd_meta_lof$estimate)
    nms<-c(nms,
           "Cadd_Meta_LOF_N_case","Cadd_Meta_LOF_N_control","Cadd_Meta_LOF_Pvalue","Cadd_Meta_LOF_OR")
    
    
    
    ## missense D( mcap>=0.05)
    for(cutoff in seq(0.025,0.2,0.025)){
        len_case=0
        len_control=0
        if(length(subcase)>0){
          len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=cutoff))]));
        }
        if(length(subcontrol)>0){
          len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=cutoff))])
        }
        test_mcap<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
        rs<-c(rs,len_case,len_control,test_mcap$p.value,test_mcap$estimate)
        nms<-c(nms,
               paste("MCAP",cutoff,"_N_case",sep=""),paste("MCAP",cutoff,"_N_control",sep=""),paste("MCAP",cutoff,"_Pvalue",sep=""),paste("MCAP",cutoff,"_OR",sep=""))
        ## missense D +LOF
        #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$MCAP)>=0.05)
        #                                                  |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
        #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$MCAP)>=0.05)
        #                                       | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
        len_case=0
        len_control=0
        if(length(subcase)>0){
          len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                                 which((as.numeric(subcase$MCAP)>=cutoff  ))) ,
                                                       which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
        }
        if(length(subcontrol)>0){
          
          len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=cutoff))
                                              ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
        }
        test_mcap_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
        rs<-c(rs,len_case,len_control,test_mcap_lof$p.value,test_mcap_lof$estimate)
        # nms<-c(nms,
        #        "MCAP_LOF_N_case","MCAP_LOF_N_control","MCAP_LOF_Pvalue","MCAP_LOF_OR")
        # }  
        # 
        nms<-c(nms,
               paste("MCAP",cutoff,"_LOF_N_case",sep=""),paste("MCAP",cutoff,"_LOF_N_control",sep=""),paste("MCAP",cutoff,"_LOF_Pvalue",sep=""),paste("MCAP",cutoff,"_LOF_OR",sep=""))
        
    }
    # ## missense D( revel>=0.75)
    # ## missense D +LOF
    # len_case=0
    # len_control=0
    # if(length(subcase)>0){
    #   
    #   len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.5))]));
    # }
    # if(length(subcontrol)>0){
    #   len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))])
    # }
    # #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
    # #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    # test_revel<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    # rs<-c(rs,len_case,len_control,test_revel$p.value,test_revel$estimate)
    # nms <- c(nms,
    #          "REVEL_N_case","REVEL_N_control","REVEL_Pvalue","REVEL_OR")
    # ## missense D +LOF
    # len_case=0
    # len_control=0
    # #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
    # #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    # if(length(subcase)>0){
    #   len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
    #                                                          which((as.numeric(subcase$REVEL)>=0.5  ))) ,
    #                                                which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
    # }
    # if(length(subcontrol)>0){
    #   
    #   len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))
    #                                       ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    # }
    # test_revel_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    # rs<-c(rs,len_case,len_control,test_revel_lof$p.value,test_revel_lof$estimate)
    # nms<-c(nms,
    #        "REVEL_LOF_N_case","REVEL_LOF_N_control","REVEL_LOF_Pvalue","REVEL_LOF_OR")
    
    ## missense D(mcap>=0.75)
    #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75)]));
    #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    
    for(cutoff in seq (0.5,1,0.05)){
    len_case=0
    len_control=0
    
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=cutoff))]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=cutoff))])
    }
    
    test_revel2<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel2$p.value,test_revel2$estimate)
    nms<-c(nms,paste("REVEL",cutoff,"_N_case",sep=""),paste("REVEL",cutoff,"_N_control",sep=""),paste("REVEL",cutoff,"_Pvalue",sep=""),paste("REVEL",cutoff,"_OR",sep=""))
    ## missense D +LOF
    # len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
    #  len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.75) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    len_case=0
    len_control=0
    if(length(subcase)>0){
      len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                             which((as.numeric(subcase$REVEL)>=cutoff  ))) ,
                                                   which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=cutoff))
                                          ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    }
    test_revel2_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel2_lof$p.value,test_revel2_lof$estimate)
    nms<-c(nms,paste("REVEL",cutoff,"_LOF_N_case",sep=""),paste("REVEL",cutoff,"_LOF_N_control",sep=""),paste("REVEL",cutoff,"_LOF_Pvalue",sep=""),paste("REVEL",cutoff,"_LOF_OR",sep=""))
  }
    
    binom<-rbind(binom,rs)
    ## missense D +LOF  
  }
  #print("done")
  binom<-data.frame(as.matrix(binom),stringsAsFactors = F,check.names = F,row.names = NULL)
  
  names(binom)<- nms
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  All_genes_binom<-rbind(All_genes_binom,binom)
  
  
  
  All_genes_binom$LOF_pvalue_adjust<-p.adjust(All_genes_binom$LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  All_genes_binom$REVEL_LOF_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  All_genes_binom$REVEL_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  outf=paste("PAH/PAH-CHD/",outname,".gnomeAD.ALL.binom.csv",sep="")
  write.csv(All_genes_binom,file = outf,row.names = F)
  All_genes_binom<-read.csv(outf,header = 1,check.names = T,stringsAsFactors = F)
  All_genes_binom<- All_genes_binom[grep("GTF2A1L",All_genes_binom$Gene,invert = T),]
  pdf(paste("PAH/PAH-CHD/",outname,".gnomeAD.QQ.pdf",sep=""),width = 5,height = 5)
  
  
  #pdf(paste("PAH/PAH-CHD/",outname,".QQ.pdf",sep=""),width = 5,height = 5)
  par(mai=c(0.5,0.5,0,0),mar=c(4,4,2,1))
  if(!is.null(All_genes_binom) && dim(All_genes_binom)[1]>0){
    
    p<-All_genes_binom
    total<-dim(All_genes_binom)[2]-6
    arr<-unlist(lapply(1:dim(All_genes_binom)[1],FUN = function(x){sum(All_genes_binom[x,seq(3,total,by=4)])}))
    sickg<-All_genes_binom$Gene[which(arr==0)]
    for(id in c(4,8,seq(28,total,by = 4))){
      key=names(p)[id]
      p<-p[order(p[,id]),]
      y=-log10(sort((p[,id]),decreasing = F))
      x=-log10((1:length(y))/length(y))
      plot(x,y,main="",xlab="",
           ylab="",pch=20,xlim=c(0,max(y,x)+0.1),
           ylim=c(0,max(y,x)+0.1),cex.axis=1.5,cex.lab=1.5)
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      title(main=key)
      title(xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",line = 2.5,cex.lab=1.5)
      if(y[1]>4){
        text(x[1],y[1],p$Gene[1],font=3,pos = 2,col="red",cex=1.5)
        points( max(x) ,max(y),pch=19,col="red",cex=2)
      }
    }
  }
  #dev.off()
  dev.off()
  
  
}


start<-function(){
  fdat<-"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.xgen_vcr_hg19.bed.gz"  #"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed"
  flist<-"PAH/Result/Data/source/PAH-CHD-list.csv"
  #flist<-"PAH/Result/Data/source/PAH_pediatric_adult_list.csv"
  PAH_case_control <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")
  filter_dat <- formatFreq_new(PAH_case_control)
  filter_dat <- filter_allfreq_local(filter_dat,0.0001,0.0001)
  #filter_dat<-COV_filter(filter_dat)
  filter_dat<-map_filter(filter_dat)
  #filter_dat<-VQSR_filter(filter_dat)
  filter_dat<- filter_dat[which(filter_dat$AF<0.001),] ## be carefule when AN<5000
  #PAH_case_control <- read.csv("PAH/Joint_calls_20170715/VCX_Control.inherited.csv_filtered_revel.csv",header = 1,sep = ",",stringsAsFactors = F,check.names = F)
  #"PAH/PAH-CHD/Data/PAH-CHD-list.csv"
  
  chd.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
  chd.list$FN<-chd.list$proband
  
  pcgc_european<-read.table("PAH/Result/Data/source/european-control.txt",stringsAsFactors = F,check.names = F,header =F )
  pcgc_european<-pcgc_european$V1
  
  pop_vcr  <-  read.table("PAH/DOC/HongJian/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
  
  vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]
  
  
  pop <- read.csv("PAH/PAH_2017/Relatedness/Phenotype_FREEZE6_pop.csv",header=1,stringsAsFactors = F,check.names = F)
  xgen_CHD_european <- intersect(chd.list$proband,pop$sampleID[which(pop$pop=="EUR")])
  
  PAH_CHD_dat<-filter_dat[which(filter_dat$ProbandName%in%c(xgen_CHD_european,vcr_chd_european)),]
  PAH_CHD_dat <-  PAH_CHD_dat[which(PAH_CHD_dat$ExonicFunc.refGene!="unknown" & PAH_CHD_dat$Func.refGene!="ncRNA_exonic"),]
  PAH_CHD_dat <-  PAH_CHD_dat[which(PAH_CHD_dat$GnomAD_Genome_cov10 > 0.9),]
  
  PAH_CHD_dat <- PAH_CHD_dat[which(PAH_CHD_dat$AN/max(PAH_CHD_dat$AN) >0.9),]
  
  #PAH_CHD_dat <- PAH_CHD_dat[which(PAH_CHD_dat$VQSLOD > -5),]
  
  PAH_CHD_dat <- PAH_CHD_dat[which(PAH_CHD_dat$MQ > 30),]
  
  PAH_CHD_dat <-  PAH_CHD_dat[which(PAH_CHD_dat$GQ>30 & PAH_CHD_dat$AD_ind >4 & PAH_CHD_dat$DP_ind > 9 & PAH_CHD_dat$AD_ind/PAH_CHD_dat$DP_ind >0.25),]
  PAH_CHD_dat<-VQSR_filter(PAH_CHD_dat)
  
  PAH_CHD_dat[which(PAH_CHD_dat$Gene.refGene=="CASC5"),"Gene.refGene"]<-"KNL1"
  #PAH_CHD_dat[PAH_CHD_dat$Gene.refGene=="CASC5"]<-
  
  
  pcgc_parent<-filter_dat[which(filter_dat$ProbandName%in%c(pcgc_european)),]
  #pah_chd<-read.table("PAH/PAH-CHD/PAH-CHD.0.0001.european.REVEL.mappablity.bed",comment.char = "",header = 1,stringsAsFactors = F,check.names = F)
  
  ### 
  
  
  ###
  
  
  gnomad_vars<-read.table("PAH/Result/Data/source/gnomad.genomes.r2.0.2.sites.NFE.0.01.anno.mapp.xgen_vcr.bed.gz",stringsAsFactors = F,check.names = F,sep="\t",quote = "\"",header = 1,comment.char = "")
  chr_vars<-gnomad_vars[which(gnomad_vars$FILTER=="PASS" & gnomad_vars$ExonicFunc.refGene!="unknown"),]
  
  chr_vars <- formatFreq_new(chr_vars)
  chr_vars <- filter_allfreq_local(chr_vars,0.0001,0.0001)
  chr_vars <- map_filter(chr_vars)
  remove <- which(chr_vars$AN< 0.9 * max(chr_vars$AN) & chr_vars$`#CHROM`!="X")
  remove <- c(remove, which((chr_vars$AN_Female < 0.9* max(chr_vars$AN_Female)  | chr_vars$AN_Male < 0.9 * max(chr_vars$AN_Male[which(chr_vars$`#CHROM`=="X")])) & chr_vars$`#CHROM`=="X"))
  if(length(remove) >0){
    chr_vars <- chr_vars[-remove,]
  }
  remove <- which(as.numeric(chr_vars$VQSLOD) < -3)
  if(length(remove) >0){
    chr_vars <- chr_vars[-remove,]
  }
  
  chr_vars<-chr_vars[which(as.numeric(chr_vars$FS) <20),]
  
  chr_vars<-chr_vars[which(as.numeric(chr_vars$QD) >3 & as.numeric(chr_vars$MQRankSum) > -2 ),]
  
  chr_vars<-chr_vars[which(as.numeric(chr_vars$InbreedingCoeff) > -0.1),]
  
  
  #remove <- intersect(which(as.numeric(chr_vars$VQSLOD) < -2),grep("non",chr_vars$ExonicFunc.refGene,invert = T))
  #if(length(remove) >0){
  #  chr_vars <- chr_vars[-remove,]
  #}
  
  
 # chr_vars<-chr_vars[which(as.numeric(chr_vars$InbreedingCoeff) > -0.1),]
  
  #chr_vars<-chr_vars[which(chr),]
  chr_vars<-chr_vars[which(chr_vars$AF < 0.001),]
  total_case=length(xgen_CHD_european)+length(vcr_chd_european)
  total_control=7509
  
  #process(total_case,total_control,PAH_CHD_dat,"PAH-pediatric_adult",chr_vars)
  print(length(grep("^syn",PAH_CHD_dat$ExonicFunc.refGene))/total_case)
  print(sum(chr_vars$AC_NFE[grep("^syn",chr_vars$ExonicFunc.refGene)])/total_control)
  
  print(length(grep("^nonsyn",PAH_CHD_dat$ExonicFunc.refGene))/total_case)
  print(sum(chr_vars$AC_NFE[grep("^nonsyn",chr_vars$ExonicFunc.refGene)])/ total_control)
  
  
  print(length(grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T))/total_case)
  print(sum(chr_vars$AC_NFE[grep("non",chr_vars$ExonicFunc.refGene,invert = T)])/total_control)
  
  process(total_case,total_control,PAH_CHD_dat,"PAH-CHD-gnomad_revision",chr_vars)
  
  
  
  result<-enrich_geneset(total_case,total_control,PAH_CHD_dat,"PAH_CHD_gnomad_revision",chr_vars) #,"ALL","ALL")
  #process(total_case,total_control,PAH_CHD_dat,"PAH-CHD",chr_vars)
  #process(1321,total_control,pcgc_parent,"PCGC-Parent",chr_vars)
  
  #chd_dat <- read.csv("PAH/PAH-CHD/CHD/CHD.WES.253known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
  #chd_dat <- filter_allfreq_local(chd_dat,0.0001,0.0001)
  #chd_dat<-COV_filter(chd_dat)
  #chd_dat<-map_filter(chd_dat)
  #chd_dat<-VQSR_filter(chd_dat)
  #names(chd_dat)[grep("proband",names(chd_dat))]<-"ProbandName"
  #process(1958,total_control,chd_dat,"CHD_cases",chr_vars)
}
start()
qqq_binom<- function(len_case,len_control){
  len_case=143
  len_control=7509
  data<-read.csv("~/server/PAH/PAH-CHD/PAH-CHD.gnomeAD.ALL.binom.csv",header = 1,check.names = F,stringsAsFactors = T,comment.char = "")
  
  pvalues<-data$REVEL_LOF_Pvalue
  pvalues[which(pvalues<0)]<-0
  plot(-log10(pvalues),-log10())
  
}