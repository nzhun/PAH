#setwd("~/server/")
setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")
source("Pipeline/Enrichment/Denovo_enrichment.R")

filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                 #    &as.numeric(data$gnomAD_genome_EAS)<freq_max
                #     &as.numeric(data$gnomAD_genome_NFE)<freq_max
                 #    &as.numeric(data$gnomAD_genome_FIN)<freq_max
                  #   &as.numeric(data$gnomAD_genome_OTH)<freq_max
                   #  &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                    # &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     #&as.numeric(data$gnomAD_genome_ALL)<freq_max
                     &as.numeric(data$gnomAD_exome_ALL)<freq_max
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
  
#  if(length(grep("VQSLOD",names(data)))>0){
#    data <- data[which(data$VQSLOD > -0.5),]
#  }
  
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
                       "VQSRTrancheSNP99.80to99.90","VQSRTrancheSNP99.70to99.80",
                       "VQSRTrancheSNP99.60to99.70",#"VQSRTrancheSNP99.50to99.60",
                       #"VQSRTrancheINDEL99.50to99.60",
                       "VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  
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
          #  print (gene)
          index_case<-which(xGen_PAH$Gene.refGene==gene)
          index_control<-which(chr_vars$Gene.refGene==gene)
          if(length(index_case)<1){len_case=0}else{subcase<-xGen_PAH[index_case,] }
          if(length(index_control)<1){len_control=0;}else{  subcontrol<-chr_vars[index_control,]}
          ## synonymous
          #if(length(index_case)<1){next}
          
        
          keys<-c("ExonicFunc.refGene","VarClass")
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[grep("^synonymous",subcase$ExonicFunc.refGene)]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol$ExonicFunc.refGene)])
          }
          test_syn<-Btest(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(gene,len_case,len_control,test_syn$p.value,test_syn$estimate)
          
          
          ### LGD
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
          }
          test_lof<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_lof$p.value,test_lof$estimate)
          
          
          ### missense D (CADD>=25)
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,  which(as.numeric(subcase$CADD_phred)>=20))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene), which(as.numeric(subcontrol$CADD13_PHRED)>=20))])
          }
          test_cadd<- Btest(len_case,total_case,len_control,total_control)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
          rs<-c(rs,len_case,len_control,test_cadd$p.value,test_cadd$estimate)
          
          ## missense D +LOF
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,
                                                               which(as.numeric(subcase$CADD13_PHRED)>=20)), which(
                                                                       subcase$ExonicFunc.refGene %in% lof_class |
                                                                       subcase$Func.refGene%in%lof_fun)
                                                       )])) ;
          }
          if(length(subcontrol)>0){
            
            len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) ,
                                                      which(as.numeric(subcontrol$CADD13_PHRED)>=20)),which( 
                                                        subcontrol$ExonicFunc.refGene%in%lof_class |
                                                        subcontrol$Func.refGene%in%lof_fun))])
          }
         test_cadd_lof<-Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_cadd_lof$p.value,test_cadd_lof$estimate)
          
          
          
          ### missense D (CADD>=25 and metasvm ==d D)
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),which(as.numeric(subcase$CADD_phred)>=20  & subcase$MetaSVM_pred=="D" ))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) , which(as.numeric(subcontrol$CADD13_PHRED)>=20 & subcontrol$MetaSVM_pred=="D"))])
          }
          test_cadd_meta<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_cadd_meta$p.value,test_cadd_meta$estimate)
          
          ## missense D +LOF
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
          
          
          
          
          ## missense D( mcap>=0.05)
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))])
          }
          test_mcap<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_mcap$p.value,test_mcap$estimate)
          
          ## missense D +LOF
          #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$MCAP)>=0.05)
         #                                                  |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
          #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$MCAP)>=0.05)
          #                                       | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
          
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                        which((as.numeric(subcase$MCAP)>=0.05  ))) ,
                                              which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
          }
          if(length(subcontrol)>0){
            
            len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))
                                              ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
          }
          test_mcap_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_mcap_lof$p.value,test_mcap_lof$estimate)
          
          
          ## missense D( revel>=0.75)
          ## missense D +LOF
          if(length(subcase)>0){
            
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.5))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))])
          }
          #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
          #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
          test_revel<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_revel$p.value,test_revel$estimate)
          
          ## missense D +LOF
          #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
          #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                                 which((as.numeric(subcase$REVEL)>=0.5  ))) ,
                                                       which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
          }
          if(length(subcontrol)>0){
            
            len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))
                                              ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
          }
          test_revel_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_revel_lof$p.value,test_revel_lof$estimate)
          
          ## missense D(mcap>=0.75)
          #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75)]));
          #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.75))]));
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))])
          }
          
          test_revel2<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_revel2$p.value,test_revel2$estimate)
          
          ## missense D +LOF
        # len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
        #  len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.75) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
          if(length(subcase)>0){
            len_case=length(unique(subcase$ProbandName[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                                 which((as.numeric(subcase$REVEL)>=0.75  ))) ,
                                                       which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
          }
          if(length(subcontrol)>0){
            len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))
                                              ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
          }
          test_revel2_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
          rs<-c(rs,len_case,len_control,test_revel2_lof$p.value,test_revel2_lof$estimate)
          
          
          binom<-rbind(binom,rs)
          ## missense D +LOF  
        }
       #print("done")
        binom<-data.frame(as.matrix(binom),stringsAsFactors = F,check.names = F,row.names = NULL)
        
        names(binom)<-c("Gene","syN_case","syN_control","syn_Pvalue","syn_OR",
                         "LOF_case","LOF_control","LOF_Pvalue","LOF_OR",
                         "CaddN_case","CaddN_control","Cadd_Pvalue","Cadd_OR",
                         "Cadd_LOF_N_case","Cadd_LOF_N_control","Cadd_LOF_Pvalue","Cadd_LOF_OR",
                         "CaddN_Meta_case","CaddN_Meta_control","Cadd_Meta_Pvalue","Cadd_Meta_OR",
                         "Cadd_Meta_LOF_N_case","Cadd_Meta_LOF_N_control","Cadd_Meta_LOF_Pvalue","Cadd_Meta_LOF_OR",
                         "MCAP_N_case","MCAP_N_control","MCAP_Pvalue","MCAP_OR",
                         "MCAP_LOF_N_case","MCAP_LOF_N_control","MCAP_LOF_Pvalue","MCAP_LOF_OR",
                         "REVEL_N_case","REVEL_N_control","REVEL_Pvalue","REVEL_OR",
                         "REVEL_LOF_N_case","REVEL_LOF_N_control","REVEL_LOF_Pvalue","REVEL_LOF_OR",
                         "REVEL2_N_case","REVEL2_N_control","REVEL2_Pvalue","REVEL2_OR",
                         "REVE2L_LOF_N_case","REVEL2_LOF_N_control","REVEL2_LOF_Pvalue","REVEL2_LOF_OR"
        )
    #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
        All_genes_binom<-rbind(All_genes_binom,binom)
    
      
      
      All_genes_binom$LOF_pvalue_adjust<-p.adjust(All_genes_binom$LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
      
      All_genes_binom$REVEL_LOF_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
      
      All_genes_binom$REVEL_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
      
      
      write.csv(All_genes_binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.ALL.binom.csv",sep=""),row.names = F)
      All_genes_binom<-read.csv("PAH/PAH-CHD/PAH-CHD.gnomeAD.ALL.binom.csv",header = 1,check.names = T,)
      pdf(paste("PAH/PAH-CHD/",outname,".gnomeAD.QQ.pdf",sep=""),width = 5,height = 5)
      par(mai=c(0.5,0.5,0,0),mar=c(4,4,2,1))
      x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
      y<--log10(sort(as.numeric(All_genes_binom$syn_Pvalue),decreasing = T))
      plot(x,y,main="SYN QQ plot",xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",xlim=c(0,max(y,x)),ylim=c(0,max(y,x)))
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      
      x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
      y<--log10(sort(as.numeric(All_genes_binom$LOF_Pvalue),decreasing = T))
      plot(x,y,main="LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",pch=20)
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
    #  points( max(x) ,max(y),pch=5,col="red")
      
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$LOF_Pvalue==min(as.numeric(All_genes_binom$LOF_Pvalue))],pos = 2,col="red")
      
      
      
      x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
      # y<--log10(sort(as.numeric(All_genes_binom$LOF_pvalue_adjust),decreasing = T))
      # plot(x,y,main="LGD QQ plot",xlab="-log10(expected)",ylab="-log10(observed_adjusted p-value)")
      # abline(0,1,lty=2,col="gray")
      # abline(h=-log10(0.05),lty=2,col="gray")
      # 
      # text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$LOF_pvalue_adjust==min(as.numeric(All_genes_binom$LOF_pvalue_adjust))],pos = 2,col="red")
      
      
      
      y<--log10(sort(as.numeric(All_genes_binom$Cadd_Pvalue),decreasing = T))
      plot(x,y,main="CADD>=20+ LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",xlim=c(0,max(y,x)),ylim=c(0,max(y,x)))
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$Cadd_Pvalue==min(as.numeric(All_genes_binom$Cadd_Pvalue))])
      
      y<--log10(sort(as.numeric(All_genes_binom$Cadd_LOF_Pvalue),decreasing = T))
      plot(x,y,main="CADD>=20 QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$Cadd_LOF_Pvalue==min(as.numeric(All_genes_binom$Cadd_LOF_Pvalue))],pos = 2,col="red")
      
      
      y<--log10(sort(as.numeric(All_genes_binom$REVEL_Pvalue),decreasing = T))
      plot(x,y,main="REVEL>=0.5 QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_Pvalue==min(as.numeric(All_genes_binom$REVEL_Pvalue))],pos = 2,col="red")
      
      
      y<--log10(sort(as.numeric(All_genes_binom$REVEL_LOF_Pvalue),decreasing = T))
      plot(x,y,main="",xlab="",
           ylab="",pch=20,xlim=c(0,max(y,x)),
           ylim=c(0,max(y,x)),cex.axis=1.5,cex.lab=1.5)
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      title(main="REVEL>0.5 and LGD QQ plot")
      title(xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",line = 2.5,cex.lab=1.5)
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_LOF_Pvalue==min(as.numeric(All_genes_binom$REVEL_LOF_Pvalue))],
           font=3,pos = 4,col="red",cex=1.5)
      points( max(x) ,max(y),pch=19,col="red",cex=2)
      
      # y<--log10(sort(as.numeric(All_genes_binom$REVEL_LOF_Pvalue_adjust),decreasing = T))
      # plot(x,y,main="REVEL>=0.5 +LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
      # abline(0,1,lty=2,col="gray")
      # abline(h=-log10(0.05),lty=2,col="gray")
      # text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_LOF_Pvalue_adjust==min(as.numeric(All_genes_binom$REVEL_LOF_Pvalue_adjust))],pos = 2,col="red")
      
      dev.off()
      
      fout=paste("PAH/Result/Data/source/",outname,".gnomeAD.11PAH.binom.csv",sep="")  #"~/PAH/PAH-CHD/PAH-CHD.gnomeAD.253CHD.binom.csv"
      fgene="PAH/Result/Data/source/PAH_associated11-13.txt" #"~/PAH/Result/Data/source/CHD.253GeneList.csv"
      
      genes<-read.csv(fgene,header=1,stringsAsFactors=F,check.names=F,comment.char="",quote="")
      sub<-All_genes_binom[which(All_genes_binom$Gene%in%genes$Gene),]
      allr<-c("Gene")
      if(dim(sub)[1]>0){
          for(i in seq(2,(dim(sub)[2]-3),by = 4)){
            len_case<-as.integer(sum(as.numeric(sub[,i])))
            len_control<-as.integer(sum(as.numeric(sub[,(i+1)])))
            allr<-c(allr,c(len_case,len_control))
            test<- Btest(len_case,total_case,len_control,total_control)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
            allr<-c(allr,c(as.numeric(test$p.value),as.numeric(test$estimate)))
          }
          allr<-c(allr,rep("-",dim(sub)[2]-length(allr)))
        #  allr<-data.frame(allr,stringsAsFactors = F)
          #names(allr)<-names(sub)
          write.csv(rbind(as.matrix(sub),allr),fout,row.names=F)
    }
      
      fout=paste("PAH/PAH-CHD/",outname,".gnomeAD.253CHD.binom.csv",sep="")  #"~/PAH/PAH-CHD/PAH-CHD.gnomeAD.253CHD.binom.csv"
      fgene="PAH/Result/Data/source/CHD.253GeneList.csv"
      
      genes<-read.csv(fgene,header=1,stringsAsFactors=F,check.names=F,comment.char="",quote="")
      sub<-All_genes_binom[which(All_genes_binom$Gene%in%genes$Gene),]
      allr<-c("Gene")
      if(dim(sub)[1]>0){
          for(i in seq(2,(dim(sub)[2]-3),by = 4)){
           # print(i)
            len_case<-as.integer(sum(as.numeric(sub[,i])))
            len_control<-as.integer(sum(as.numeric(sub[,(i+1)])))
            allr<-c(allr,c(len_case,len_control))
            test<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
            allr<-c(allr,c(as.numeric(test$p.value),as.numeric(test$estimate)))
          }
          allr<-c(allr,rep("-",dim(sub)[2]-length(allr)))
        #  allr<-data.frame(allr,stringsAsFactors = F)
        #  names(allr)<-names(sub)
          write.csv(rbind(as.matrix(sub),allr),fout,row.names=F)   
          print(fout)
    #  }

      }
}
start<-function(){
      fdat<-"PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed"
      flist<-"PAH/Result/Data/source/PAH-CHD-list.csv"
      #flist<-"PAH/Result/Data/source/PAH_pediatric_adult_list.csv"
      PAH_case_control <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")
      filter_dat <- formatFreq_new(PAH_case_control)
      filter_dat <- filter_allfreq_local(filter_dat,0.0001,0.0001)
      #filter_dat<-COV_filter(filter_dat)
      filter_dat<-map_filter(filter_dat)
      filter_dat<-VQSR_filter(filter_dat)
      filter_dat<- filter_dat[which(filter_dat$AF<0.001),] ## be carefule when AN<5000
      #PAH_case_control <- read.csv("PAH/Joint_calls_20170715/VCX_Control.inherited.csv_filtered_revel.csv",header = 1,sep = ",",stringsAsFactors = F,check.names = F)
      #"PAH/PAH-CHD/Data/PAH-CHD-list.csv"
      
      chd.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
      chd.list$FN<-chd.list$proband
      
      pcgc_european<-read.table("PAH/Result/Data/source/european-control.txt",stringsAsFactors = F,check.names = F,header =F )
      pcgc_european<-pcgc_european$V1
      
      pop_vcr  <-  read.table("PAH/Result/Data/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
      
      vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]
      
      
      pop <- read.csv("PAH/Result/Data/source/Phenotype_FREEZE6_pop.csv",header=1,stringsAsFactors = F,check.names = F)
      xgen_CHD_european <- intersect(chd.list$proband,pop$sampleID[which(pop$pop=="EUR")])
      
      PAH_CHD_dat<-filter_dat[which(filter_dat$ProbandName%in%c(xgen_CHD_european,vcr_chd_european)),]
      pcgc_parent<-filter_dat[which(filter_dat$ProbandName%in%c(pcgc_european)),]
      #pah_chd<-read.table("PAH/PAH-CHD/PAH-CHD.0.0001.european.REVEL.mappablity.bed",comment.char = "",header = 1,stringsAsFactors = F,check.names = F)
      
      chr_vars<-read.table("PAH/Result/Data/source/gnomad.genomes.r2.0.1.sites.all.anno.0911.bed.gz",stringsAsFactors = F,check.names = F,sep="\t",header = 1,comment.char = "")
      chr_vars<-chr_vars[which(chr_vars$FILTER=="PASS"),]
      chr_vars <- formatFreq_new(chr_vars)
      chr_vars <- filter_allfreq_local(chr_vars,0.0001,0.0001)
      
      total_case=length(xgen_CHD_european)+length(vcr_chd_european)
      total_control=7509
      
      #process(total_case,total_control,PAH_CHD_dat,"PAH-pediatric_adult",chr_vars)
      
      process(total_case,total_control,PAH_CHD_dat,"PAH-CHD",chr_vars)
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
  data<-read.csv("PAH/Result/Data/output/PAH-CHD.gnomeAD.ALL.binom.csv",header = 1,check.names = F,stringsAsFactors = T,comment.char = "")

  pvalues<-data$REVEL_LOF_Pvalue
 pvalues[which(pvalues<0)]<-0
  plot(-log10(pvalues),-log10())
  
}