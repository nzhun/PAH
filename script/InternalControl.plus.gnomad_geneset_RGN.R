setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")
source("Pipeline/Enrichment/Denovo_enrichment.R")
source("PAH/PAH_10032017/script/filter_CCHMC.R")


qcase<-function(data){
  return(dim(unique(data[,c("ProbandName","Gene.refGene")]))[1])
}

test_geneset<-function(total_case,total_control,xGen_PAH,chr_vars,genesets,correct){
  keys<-c("ExonicFunc.refGene","VarClass")
  index1<-grep("ExonicFunc.refGene",names( xGen_PAH));
  if(length(index1)<1){index1<-grep("VarClass",names(xGen_PAH))};
  index2<-grep("^Func.refGene",names(xGen_PAH));
  if(length(index2)<1){index2<-grep("VarFunc",names(xGen_PAH))};
  
  ct_index1<-grep("ExonicFunc.refGene",names( chr_vars));
  if(length(ct_index1)<1){index1<-grep("VarClass",names(chr_vars))};
  ct_index2<-grep("^Func.refGene",names(chr_vars));
  if(length(ct_index2)<1){index2<-grep("VarFunc",names(chr_vars))};
  
  types<-unique(c(xGen_PAH[,index1],chr_vars[,ct_index1]))
  
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
  
  if(is.null(subcase)||is.null(subcontrol)){return()}
  
  if(dim(subcase)[1]>0 && length(index1)>0){
    #   print(index1)
    #  print(names(subcase)[index1])
    len_case=qcase(subcase[grep("^synonymous",subcase[,index1]),c("ProbandName","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0 && length(ct_index1)>0){
    len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol[,ct_index1])])
  }
  test_syn<-Btest(len_case,total_case,len_control,total_control,correct) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-c("SYN",len_case,len_control,test_syn$p.value,test_syn$estimate)
  
  ### LGD
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[which(subcase[,index1] %in% lof_class |subcase[,index1]%in%lof_fun),c("ProbandName","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0){
    len_control=sum(subcontrol$AC_NFE[which(subcontrol[,ct_index1] %in%lof_class |subcontrol[,ct_index2]%in%lof_fun)])
  }
  test_lof<-Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("LGD",len_case,len_control,test_lof$p.value,test_lof$estimate))
  
  
  # mis
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[grep("^nonsynonymous",subcase[,index1]),c("ProbandName","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0){
    len_control=sum(subcontrol$AC_NFE[grep("^nonsynonymous",subcontrol[,ct_index1])])
  }
  test_mis<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
  rs<-rbind(rs,c("MIS",len_case,len_control,test_mis$p.value,test_mis$estimate))
  
  
  ### missense D (CADD>=25)
  coln<-grep("CADD.*ph",names(subcase),ignore.case = T)
  
  coln2<-grep("CADD.*ph",names(subcontrol),ignore.case = T)
  for(cadd in seq(5,35,by = 5)){
    if(dim(subcase)[1]>0){
      len_case=qcase(subcase[intersect(grep("^nonsynonymous",subcase[,index1]) ,  
                                       which(as.numeric(subcase[,coln])>=cadd)),c("ProbandName","Gene.refGene")]);
    }
    if(dim(subcontrol)[1]>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]), which(as.numeric(subcontrol[,coln2])>=cadd))])
    }
    test_cadd<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
    rs<-rbind(rs,c(paste("CADD",cadd,sep=""),len_case,len_control,test_cadd$p.value,test_cadd$estimate))
  }
  ## missense D +LOF
  #for(cadd in seq(15,30,by = 5)){
  cadd=25
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[c(intersect(grep("^nonsynonymous",subcase[,index1]) ,
                                       which(as.numeric(subcase[,coln])>=cadd)), which(
                                         subcase$ExonicFunc.refGene %in% lof_class |
                                           subcase$Func.refGene%in%lof_fun)
    ),c("ProbandName","Gene.refGene")]) ;
  }
  if(dim(subcontrol)[1]>0){
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol[,ct_index1]) ,
                                                  which(as.numeric(subcontrol[,coln2])>=cadd)),which( 
                                                    subcontrol$ExonicFunc.refGene%in%lof_class |
                                                      subcontrol[,ct_index2]%in%lof_fun))])
  }
  test_cadd_lof<-Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+LGD",len_case,len_control,test_cadd_lof$p.value,test_cadd_lof$estimate))
  # }
  
  
  ### missense D (CADD>=25 and metasvm ==d D)
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[intersect(grep("^nonsynonymous",subcase[,index1]),
                                     which(as.numeric(subcase[,coln])>=25  & subcase$MetaSVM_pred=="D" )),c("ProbandName","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]) ,
                                                which(as.numeric(subcontrol[,coln2])>=25 & subcontrol$MetaSVM_pred=="D"))])
  }
  test_cadd_meta<-Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+METASVM",len_case,len_control,test_cadd_meta$p.value,test_cadd_meta$estimate))
  
  ## missense D +LOF
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[c(intersect(grep("^nonsynonymous",subcase[,index1]),
                                       which((as.numeric(subcase[,coln])>=25  & subcase$MetaSVM_pred=="D"))) ,
                             which(subcase[,index1] %in% lof_class |subcase[,index2]%in%lof_fun)),c("ProbandName","Gene.refGene")])[1];
    
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),
                                                  which( as.numeric(subcontrol[,coln2])>=25   & subcontrol$MetaSVM_pred=="D"))
                                        ,which(subcontrol[,ct_index1]%in%lof_class |subcontrol[,ct_index2]%in%lof_fun))])
  }
  
  test_cadd_meta_lof<- Btest(len_case,total_case,len_control,total_control,correct)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DCADD25+METASVM+LGD",len_case,len_control,test_cadd_meta_lof$p.value,test_cadd_meta_lof$estimate))
  
  
  
  
  ## missense D( mcap>=0.05)
  for(thred in seq(0.025,0.2,by = 0.025)){
    if(length(subcase)>0){
      index<-intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$MCAP)>=thred))
      len_case<-qcase(subcase[index,]);
      #len_case=dim(unique(subcase$ProbandName[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05)),c("ProbandName","Gene.refGene")]))[1];
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$MCAP)>=thred))])
    }
    test_mcap<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-rbind(rs,c(paste("mcap>=",thred,sep=""),len_case,len_control,test_mcap$p.value,test_mcap$estimate))
  }
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
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$MCAP)>=0.05))
                                        ,which(subcontrol[,ct_index1]%in%lof_class |subcontrol[,ct_index2]%in%lof_fun))])
  }
  test_mcap_lof<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("MCAP>=0.05+LGD",len_case,len_control,test_mcap_lof$p.value,test_mcap_lof$estimate))
  
  
  ## missense D( revel>=0.75)
  ## missense D +LOF
  if(length(subcase)>0){
    
    len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$REVEL)>=0.5)),]));
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=0.5))])
  }
  #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
  #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
  test_revel<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL>=0.5",len_case,len_control,test_revel$p.value,test_revel$estimate))
  
  for(thred in seq(0.55,0.9,by = 0.05)){
    if(length(subcase)>0){
      
      len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$REVEL)>=thred)),]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=thred))])
    }
    #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
    #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    test_revel<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-rbind(rs,c(paste("DREVEL>=",thred,sep=""),len_case,len_control,test_revel$p.value,test_revel$estimate))
  }
  
  ## missense D +LOF
  #len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
  #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  if(length(subcase)>0){
    len_case=qcase((subcase[c(intersect(grep("nonsynonymous",subcase[,index1]),
                                        which((as.numeric(subcase$REVEL)>=0.5  ))) ,
                              which(subcase[,index1] %in% lof_class |subcase[,index2]%in%lof_fun)),]) );
  }
  if(length(subcontrol)>0){
    
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=0.5))
                                        ,which(subcontrol[,ct_index1]%in%lof_class |subcontrol[,ct_index2]%in%lof_fun))])
  }
  test_revel_lof<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL0.5+LGD",len_case,len_control,test_revel_lof$p.value,test_revel_lof$estimate))
  
  ## missense D(mcap>=0.75)
  #len_case=length(unique(subcase$ProbandName[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75)]));
  #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
  if(length(subcase)>0){
    len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$REVEL)>=0.75)),]));
    #  len_case=length(unique(subcase[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.75)),c("ProbandName","Gene.refGene")]));
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=0.75))])
  }
  
  test_revel2<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("DREVEL0.75",len_case,len_control,test_revel2$p.value,test_revel2$estimate))
  
  ## missense D +LOF
  # len_case=length(unique(subcase$ProbandName[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
  #  len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.75) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
  if(length(subcase)>0){
    len_case=qcase((subcase[c(intersect(grep("nonsynonymous",subcase[,index1]),
                                        which((as.numeric(subcase$REVEL)>=0.75  ))) ,
                              which(subcase[,index1] %in% lof_class |subcase[,index2]%in%lof_fun)),]) );
  }
  if(length(subcontrol)>0){
    len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=0.75))
                                        ,which(subcontrol[,ct_index1]%in%lof_class |subcontrol[,ct_index2]%in%lof_fun))])
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
  if(!is.null(sbinom)){
    #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  # correct<-as.numeric(sbinom[which(sbinom$Type=="SYN"),"OR"])
  correct<-1
  print(correct)
  genesets="ALL";
  group="ALL"
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  print(All_genes_binom)
  
  ## chanel ##
  
  #### chanel ### 
  chanelgenes<-read.table("PAH/documents/Channelopathy_gene_list.csv",stringsAsFactors = F,sep="\t")
  genesets<-chanelgenes;
  group=paste(dim(chanelgenes)[1],"chanel genes")
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  
  
  ####### 25
  hhe_genesets=sets$gene_hhe;
  group="Heart_top25"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhe_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hhe_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  hle_genesets=sets$gene_hle;
  group="lung_top25"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hle_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hle_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  #### endothelial 
  hendoth_genesets=sets$gene_high_endothelial;
  group="Endothelial_top25"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hendoth_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hendoth_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  #### endothelial 
  hendoth_genesets=sets$gene_high_endothelial_heart;
  group="Endothelial_top25_heart"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hendoth_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hendoth_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  #### endothelial 
  hendoth_genesets=sets$gene_high_endothelial_lung;
  group="Endothelial_top25_lung"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hendoth_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hendoth_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  ## lung and heart
  
  
  he_genesets=sets$gene_he;
  group="lung&Heart_top25"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,he_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(he_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  hhle_genesets=sets$gene_hhle;
  group="heart or lung_top25"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhle_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hhle_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
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
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hle_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  
  he_50_genesets=sets$gene_he_50;
  group="lung&Heart_top50"
  # All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,he_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(he_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  hhle_50_genesets=sets$gene_hhle_50;
  group="heart or lung_top50"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,hhle_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(hhle_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ##############################
  
  
  
  misz_genesets=sets$misz_sets;
  group="misZ >3"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(misz_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  pli_genesets=sets$pli_sets;
  #  lofz_1.64=sets$
  group="pLi >0.9"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,pli_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(pli_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  genesets=sets$lofz_sets;
  group="lofZ >3"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  chdgenes<-read.csv("PAH/Result/Data/source/CHD.253GeneList.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-chdgenes$Gene;
  group="253CHD"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  
  genesets=sets$lofz_1.64_sets;
  group="lofZ >1.64"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  chdgenes<-read.csv("PAH/Result/Data/source/CHD.253GeneList.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-chdgenes$Gene;
  group="253CHD"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  pahgenes<-read.csv("PAH/Result/Data/source/PAH_associated11-13.txt",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-pahgenes$Gene;
  genesets<-genesets[which(genesets!="EIF2AK4")]
  group="11 PAH genes"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  pahgenes<-read.csv("PAH/Result/Data/source/PAH_associated11-13.txt",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",strip.white = T,fill = T)
  genesets<-pahgenes$Gene;
  genesets<-genesets[which(genesets!="EIF2AK4" & genesets!="BMPR2")]
  group="10 PAH genes-BMPR2"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  genesets<- c(pahgenes$Gene,"ABCC8","SOX17");
  group="13 PAH genes"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  genesets<- c(pahgenes$Gene,"ABCC8","SOX17");
  genesets<-genesets[which(genesets!="EIF2AK4" & genesets!="BMPR2")]
  group="13 PAH genes-BMPR2"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]
  group="SOX17_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ## SOX17 and PLI >0.9
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_pli_genesets<- intersect(SOX17_genesets,pli_genesets)
  group="SOX17_pLi_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  SOX17_pli_genesets<- intersect(SOX17_genesets,sets$lofz_1.64_sets)
  group="SOX17_lof1.64_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  # SOX17_pli_genesets<- intersect(SOX17_genesets,)
  # group="SOX17_pLi_target"
  # #All_genes_binom<-c()
  # sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
  # #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  # if(!is.null(sbinom)){
  #   sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
  #   All_genes_binom<-rbind(All_genes_binom,sbinom)
  # }
  
  ## SOX17 and misZ>3
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_misz_genesets<- intersect(SOX17_genesets,misz_genesets)
  group="SOX17_misz_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_misz_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
  group="SOX17_lung25_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  ## SOX17 and heart top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_heart_genesets<- intersect(SOX17_genesets,hhe_genesets)
  group="SOX17_heart25_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_heart_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  SOX17_lung_heart_genesets<- intersect(SOX17_genesets,he_genesets)
  group="SOX17_lung_heart25_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_heart_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  ############### top 50 ###################
  
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_50_genesets<- intersect(SOX17_genesets,hle_50_genesets)
  group="SOX17_lung50_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  ## SOX17 and heart top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_heart_50_genesets<- intersect(SOX17_genesets,hhe_50_genesets)
  group="SOX17_heart50_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_heart_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  SOX17_lung_heart_50_genesets<- intersect(SOX17_genesets,he_50_genesets)
  group="SOX17_lung_heart50_target"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_heart_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  ####### Exclude SOX17 #####
  
  
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]
  SOX17_genesets<-SOX17_genesets[which(SOX17_genesets!="SOX17")]
  group="SOX17_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ## SOX17 and PLI >0.9
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_pli_genesets<- intersect(SOX17_genesets,pli_genesets)
  group="SOX17_pLi_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_pli_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_pli_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  ## SOX17 and misZ>3
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_misz_genesets<- intersect(SOX17_genesets,misz_genesets)
  group="SOX17_misz_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_misz_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_misz_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_genesets<- intersect(SOX17_genesets,hle_genesets)
  group="SOX17_lung25_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  ## SOX17 and heart top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_heart_genesets<- intersect(SOX17_genesets,hhe_genesets)
  group="SOX17_heart25_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_heart_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  SOX17_lung_heart_genesets<- intersect(SOX17_genesets,he_genesets)
  group="SOX17_lung_heart25_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_genesets,correct)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_heart_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  ############### top 50 ###################
  
  ## SOX17 and lung top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_lung_50_genesets<- intersect(SOX17_genesets,hle_50_genesets)
  group="SOX17_lung50_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  
  ## SOX17 and heart top25
  #sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
  SOX17_heart_50_genesets<- intersect(SOX17_genesets,hhe_50_genesets)
  group="SOX17_heart50_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_heart_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_heart_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  SOX17_lung_heart_50_genesets<- intersect(SOX17_genesets,he_50_genesets)
  group="SOX17_lung_heart50_target-1"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,SOX17_lung_heart_50_genesets,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(SOX17_lung_heart_50_genesets),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  abcc8<- c("ABBC8","ABBC8")
  group="abcc8"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,abcc8,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(abcc8),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  
  abcc8<- c("SOX17","SOX17")
  group="SOX17"
  #All_genes_binom<-c()
  sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,abcc8,correct)
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  if(!is.null(sbinom)){
    sbinom$Genesets<-paste(length(abcc8),group)
    All_genes_binom<-rbind(All_genes_binom,sbinom)
  }
  ## before 
  # 
  # sbinom<-test_geneset(total_case,total_control,xGen_PAH,chr_vars,genesets,1)
  # #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  # sbinom$Genesets<-paste(length(genesets),group,"_raw")
  # All_genes_binom<-rbind(All_genes_binom,sbinom)
  
  
  write.csv(All_genes_binom,file = paste("PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/",outname,".gnomeAD.genesets.binom.csv",sep=""),row.names = F)
  return (All_genes_binom)
}


plots<-function(){ 
  png("PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/rate.gene.case.control.png")  
  cmb<-c();for(g in intersect(unique(PAH_CHD_dat$Gene.refGene),unique(chr_vars$Gene.refGene))){cmb<-rbind(cmb,c(g,length(which(PAH_CHD_dat$Gene.refGene==g))/total_case, length(which(chr_vars$Gene.refGene==g))/total_control ) )}
  
  plot(-log10(as.numeric(cmb[,2])),-log10(as.numeric(cmb[,3])),xlab="cases",ylab="gnomad")
  abline(0,1)
  dev.off()
  
  
  png("PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/gene.case.control.OR.png")  
  #cmb<-c();for(g in intersect(unique(PAH_CHD_dat$Gene.refGene),unique(chr_vars$Gene.refGene))){cmb<-rbind(cmb,c(g,length(which(PAH_CHD_dat$Gene.refGene==g))/total_case, length(which(chr_vars$Gene.refGene==g))/total_control ) )}
  yy<--log10((as.numeric(cmb[,2])+0.0001)/(0.0001+as.numeric(cmb[,3])))
  plot(1:dim(cmb)[1],yy,xlab="cases",ylab="gnomad")
  abline(0,1)
  abline(h=0.5,col="yellow")
  abline(h=-0.5,col="yellow")
  
  dev.off()
  
  fgenes<-cmb[which(yy  >0.75|yy < -0.75),1]
  write.csv(PAH_CHD_dat[which(PAH_CHD_dat$Gene.refGene%in%fgenes),],file = "PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/PAH_CHD_dat.falsegenes.csv")
  write.csv(chr_vars[which(chr_vars$Gene.refGene%in%fgenes),],file = "PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/gnomad.falsepositive.csv")
}

#single_variant



#start<-function(){
    fdat="PAH/PAH_10032017/hg38/CCHMC_Freeze_One.GL.pVCF.AF0.01.anno.bravo.gnomad.intersection.target.IMAF.bed.gz" ## regeneron
    flist<-"PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped"
    ctr_file="Resources/Control/InHouse_Control/hg38/Filter/Inhouse_control.filter.european.target.laf.reannoBRAVO.bed.gz" #gnomad.2.0.2.exon.hg38.bravo.bed.gz" #"Resources/gnomad/gnomad.genomes.refgene.base1.hg38.bravo-3.Exac-3.esp-3.1kg-3.target.bed.gz" #"Resources/gnomad/gnomad.genomes.refgene.base1.hg38.bravo-4.Exac-4.esp-3.1kg-3.NFE.bed.gz"
   
   # fgnomad<-"Resources/gnomad/LA/gnomad.genomes.r2.0.2.hg38.anno.bravo.intersection.bed.gz"
    fgnomad<-"Resources/gnomad/LA/gnomad.genomes.r2.0.2.hg38.anno.bravo.intersection.bed.gz"
     
    fexclude<-"PAH/PAH_10032017/src/Excluding_CCHMC_asso.txt"
    exclude<-read.csv(fexclude,header=1,stringsAsFactors = F,check.names = F,comment.char = "")
    

    
    chd.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
    if(dim(chd.list)[2]<2){
      chd.list<-read.table(flist,header = 1,sep="\t",stringsAsFactors = F,strip.white = T,check.names = F,comment.char = "",quote = "")
    }
    chd.list$proband<-chd.list$ID
    chd.list$FN<-chd.list$proband
    chd.list<-chd.list[which(!chd.list$proband%in%exclude$`#FAM`),]
    
    pcgc_european<-read.table("Resources/Control/InHouse_Control/src/EUR_control_4597.txt",stringsAsFactors = F,check.names = F,header =F )
    pcgc_european<-pcgc_european$V1
    
    # pop_vcr  <-  read.table("PAH/Result/Data/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
    vcr_chd_european<-c()
    # vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]
    
    
    pop <- read.csv("PAH/PAH_10032017/src/PAH.pop.csv",header=1,stringsAsFactors = F,check.names = F)
    pop$sampleID<-pop$proband
    xgen_CHD_european <- intersect(chd.list$proband,pop$sampleID[which(pop$pop=="EUR")])
      
    total_case=length(xgen_CHD_european)
    total_control=length(pcgc_european)
    
    
    gnomad<-read.table(fgnomad, stringsAsFactors = F,check.names = F,sep="\t",header = 1,comment.char = "",quote="")
    ctr_dat<-read.table(ctr_file,stringsAsFactors = F,check.names = F,sep="\t",header = 1,comment.char = "",quote = "")
    PAH_case_control <- read.table(fdat,header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "",quote = "")

    if(length(grep("ProbandName",names(PAH_case_control)))<1){
      PAH_case_control$ProbandName<-unlist(lapply(PAH_case_control$proband,FUN = function(x) unlist(strsplit(x,split = "(",fixed = T))[1]  ))
    }
    
    if(length(grep("ProbandName",names(ctr_dat)))<1){
      ctr_dat$ProbandName<-unlist(lapply(ctr_dat$proband,FUN = function(x) unlist(strsplit(x,split = "(",fixed = T))[1]  ))
    }
    
    

    
    
    case_filter_dat<- PAH_case_control[which( PAH_case_control$ProbandName%in% xgen_CHD_european &
                                           PAH_case_control$ExonicFunc.refGene!="unknown" &
                                           PAH_case_control$`#CHROM`!="Y" & PAH_case_control$gnomad_COV10!="." & PAH_case_control$gnomad_COV15!="."),]
    
    filter_dat<- case_filter_dat[which( case_filter_dat$`#CHROM`!="X" ),]
    
    ctr_dat<- ctr_dat[which(ctr_dat$ExonicFunc.refGene!="unknown"),]
    
    case_dat<-filter_case_inH_merge( filter_dat,total_case)
    PAH_CHD_dat<-case_dat$data
    PAH_CHD_filter<-case_dat$filter
    print(c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
            dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
            dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
            dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) )
    
    write.csv(PAH_CHD_filter,file="PAH/PAH_10032017/Result/CCHMC.comb.filter.csv")
    
   ctr_dat<-  ctr_dat[which(ctr_dat[,1]!="X"),]
   gdat<-filter_control_merge(ctr_dat,length(pcgc_european))
   ctr_vars<-gdat$data
   ctr_filter_rs<-gdat$filter #$1==5 && $2==353888 || ($1==15 && $2==22866769
   remove<-which((ctr_vars[,1]==5 & ctr_vars[,2]==353888) |(ctr_vars[,1]==15 & ctr_vars[,2]==22866769))
   if(length(remove)>0){
    	ctr_vars<-ctr_vars[-remove,]
   }
    
   print(c(dim( ctr_vars)[1]/total_control,dim(ctr_vars[grep("^synony", ctr_vars$ExonicFunc.refGene),])[1],
           dim( ctr_vars[grep("^synony", ctr_vars$ExonicFunc.refGene),])[1]/total_control,
           dim( ctr_vars[grep("^nonsynony", ctr_vars$ExonicFunc.refGene),])[1]/total_control,
           dim(ctr_vars[grep("non",ctr_vars$ExonicFunc.refGene,invert = T),])[1]/total_control) )
   
   write.csv(ctr_filter_rs,file="PAH/PAH_10032017/Result/internalcontrol.comb.filter.csv")
   
   
   gnomad<-gnomad[which(gnomad[,1]!="X"),]
   gnomad<-gnomad[which(gnomad$Gene.refGene!="LOC100134391"),]
   gdat<-filter_gnomad_merge(gnomad)
   gnomad_dat<-gdat$data
   gfilter_rs<-gdat$filter
   total_gnomad=7509
   
   print(c(sum( gnomad_dat$AC_NFE)/total_gnomad,sum(gnomad_dat$AC_NFE[grep("^synony", gnomad_dat$ExonicFunc.refGene)]),
           sum( gnomad_dat[grep("^synony", gnomad_dat$ExonicFunc.refGene),"AC_NFE"])/total_gnomad,
           sum( gnomad_dat[grep("^nonsynony", gnomad_dat$ExonicFunc.refGene),"AC_NFE"])/total_gnomad,
           sum(gnomad_dat[grep("non",gnomad_dat$ExonicFunc.refGene,invert = T),"AC_NFE"])[1]/total_gnomad) )
   write.csv(gfilter_rs,file="PAH/PAH_10032017/Result/gnomad.comb.filter.csv")
   
   
   #exclude_snv<-read.csv("PAH/PAH_10032017/Result/Exclude.SNV.txt",header=1,stringsAsFactors = F)
  index<-which(PAH_CHD_dat[,1] ==7 & PAH_CHD_dat$FROM==48644620 &PAH_CHD_dat$REF=="C" & PAH_CHD_dat$ALT=="A");
  if(length(index)>0){PAH_CHD_dat<-PAH_CHD_dat[-index,]}
  index<-which(ctr_vars[,1] ==7 & ctr_vars$FROM==48644620 &ctr_vars$REF=="C" & ctr_vars$ALT=="A");
  if(length(index)>0){ctr_vars<-ctr_vars[-index,]}

  ctr_vars$AC_NFE=1;
  
  nms<-intersect(names(ctr_vars),names(gnomad_dat))
  chr_vars<-rbind(ctr_vars[,nms],gnomad_dat[,nms])
  
  
    index<-grep("CHR",names(chr_vars),ignore.case = T)
    auto_gnomad<-chr_vars[which(chr_vars[,index]!="X"),]
    index<-grep("CHR",names(PAH_CHD_dat),ignore.case = T)
    auto_pah<-PAH_CHD_dat[which(PAH_CHD_dat[,index]!="X" ),]
    #names(chd_dat)[grep("proband",names(chd_dat))]<-"ProbandName"
     total_cmb <- total_control+total_gnomad
    result<-process(total_case,total_cmb,PAH_CHD_dat,"PAH_Cin_hg38_all_RGN_Inhouse_plus_gnomad",chr_vars) #,"ALL","ALL")
    result<-process(total_case,total_cmb,auto_pah,"PAH_Cin_hg38_Auto_0.15_Inhouse_plus_gnomad",auto_gnomad) #,"ALL","ALL")
#}
#start()
   #plots()
   
   gnomad_dat$proband <- "."
   PAH_CHD_dat$AC_NFE <- 1
   PAH_CHD_dat$Hom_NFE <- 1
   ctr_vars$Hom_NFE <- 1
   gnomad_dat$Is_homo <- F
   gnomad_dat$Genotype <- "."
  # PAH_CHD_dat
   
   nms<-intersect(names(ctr_vars),names(gnomad_dat))
   chr_vars<-rbind(ctr_vars[,nms],gnomad_dat[,nms])
   
   nm<-intersect(names(PAH_CHD_dat),names(chr_vars))
   write.table(rbind(PAH_CHD_dat[,nm],chr_vars[,nm]),file="PAH/PAH_10032017/hg38/merged.filter.CCHMC.Internal.plus.gnomad.corrected.bed",sep="\t",row.names=F,quote=F)
