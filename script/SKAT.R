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
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  if(length(grep("genomicSuperDups",names(data)))>0){
    index<-grep("Score",data$genomicSuperDups)
    as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return(data)
}


VQSR_filter <- function(data){
  if(length(grep("FILTER",names(data)))>0){
    filter_ignore <- c("VQSRTrancheSNP99.90to100.00",
                      "VQSRTrancheSNP99.70to99.80","VQSRTrancheSNP99.80to99.90",
                       #"VQSRTrancheSNP99.60to99.70","VQSRTrancheSNP99.50to99.60","VQSRTrancheINDEL99.50to99.60","VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  return(data)
}

COV_filter <- function(data){
  if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
    data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.8 & 2*data$AC_xGen/data$AN_xGen >0.8 & 2*(data$AC_VCR/data$AN_VCR) >0.8  &data$GnomAD_Genome_cov10>0.8 ),]
  }
  if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
    data <- data[which(2*data$AC_Control/data$AN_Control >0.8 & 2*data$AC_VU/data$AN_VU >0.8 ),]
  }
  return(data)
}




assocition <- function(xGen_PAH,VCR_PAH,control_dat){
  
  ### association test fisher exac test ###
  ## only work on european
  ## filter: popfreq<10-3, cadd 20
  europeans <- unique(c(pop_vcr$V1,euro_PAH,VU_euro))
  len_case <- length(unique(PAH_case_control$ProbandName[which(PAH_case_control$ProbandName %in% europeans)]))+
    length(unique(SPH_cases$ProbandName[which(SPH_cases$ProbandName %in%  europeans)])) +
    length(unique(VU_cases$ProbandName[which(VU_cases$ProbandName %in%  europeans)]))
  #    len_case <- length(unique(xGen_PAH$ProbandName)) +length(unique(VCR_PAH$ProbandName))
  len_control <- length(unique(control_dat$ProbandName))
  list_nm <- intersect(names(PAH_case_control),intersect(names(SPH_cases),names(VU_cases)))
  PAH_case_control <- formatFreq(PAH_case_control)
  SPH_cases <- formatFreq(SPH_cases)
  VU_cases <- formatFreq(VU_cases)
  control_dat <- formatFreq(control_dat)
  
  control_dat <-  VQSR_filter(control_dat)
  control_dat <- filter_allfreq(control_dat)
  control_dat <- COV_filter(control_dat)
  
  PAH_case_control <-  VQSR_filter(PAH_case_control)
  PAH_case_control <-  filter_allfreq(PAH_case_control)
  PAH_case_control <-  COV_filter(PAH_case_control)
  
  SPH_cases <-  VQSR_filter(SPH_cases)
  SPH_cases <- filter_allfreq(SPH_cases)
  SPH_cases <- COV_filter(SPH_cases)
  
  VU_cases <- VQSR_filter(VU_cases)
  VU_cases <- filter_allfreq(VU_cases)
  VU_cases <- COV_filter(VU_cases)
  
  variants_asso <-  rbind(
    PAH_case_control[unique(c(which(as.numeric(PAH_case_control$ExACfreq)<0.001 & (as.numeric(PAH_case_control$CADDphred)>20) 
                                    & PAH_case_control$VarClass != "synonymousSNV"),get_lgd(PAH_case_control))),list_nm],
    SPH_cases[unique(c(which(as.numeric(SPH_cases$ExACfreq)<0.001 & (as.numeric(SPH_cases$CADDphred)>20)
                             &SPH_cases$VarClass != "synonymousSNV"),get_lgd(SPH_cases))),list_nm],
    VU_cases[unique(c(which(as.numeric(VU_cases$ExACfreq)<0.001 & (as.numeric(VU_cases$CADDphred)>20)
                            &VU_cases$VarClass != "synonymousSNV"),get_lgd(VU_cases))),list_nm]
    
  )
  variants_asso <- variants_asso[which(variants_asso$ProbandName %in% europeans),]
  control_asso <- control_dat[unique(c(which(as.numeric(control_dat$ExACfreq) < 0.001 & as.numeric(control_dat$CADDphred) > 20 
                                             & control_dat$VarClass != "synonymousSNV"),get_lgd(control_dat))),]
  rs <- c()
  for( gene in unique(variants_asso$GeneName)){
    case_count=length(unique(variants_asso$ProbandName[which(variants_asso$GeneName==gene)]))
    control_count=length(unique(control_asso$ProbandName[which(control_asso$GeneName==gene)]))
    fs <- fisher.test(matrix(c(case_count,len_case,control_count,len_control),nrow = 2,byrow = T),alternative = "greater")
    rs <- rbind(rs,c(gene,case_count,len_case,control_count,len_control,fs$p.value))
  }
  
  rs <- data.frame(rs,stringsAsFactors = F)
  names(rs) <- c("Gene","Ncase","Total_case","Ncontrol","Total_control","pvalue")
  rs <- rs[order(as.numeric(rs$pvalue)),]
  write.csv(rs,"PAH/PAH-CHD/case-control-assocition.csv",quote = F,row.names = F)
  png("PAH/PAH-CHD/case_control_association.png", width = 7, height = 7, units = 'in', res = 500)
  plot(-log10((1:dim(rs)[1])/dim(rs)[1]),-log10(as.numeric(as.character(rs$pvalue))),
       ylim=c(0,-log10(as.numeric(rs$pvalue[1]))*1.1),xlab="Expected",ylab="Observed",main="Case-control-European-PAH")
  abline(0,1,lty=2)
  text(-log10((1:3)/dim(rs)[1]),-log10(as.numeric(rs$pvalue[1:3])),labels = rs$Gene[1:3],pos=c(1,3,2),cex=0.75)
  dev.off()
  return(rs)
}



skat<-function(map_case,map_control,case_D,control_D){
  ## load Lookup table, get the weight
  het_lookup<-read.table("Application/psap/lookups/full.het.pCADD.gencodeV19.allsites.txt.gz",stringsAsFactors = F)
  lof_lookup<-read.table("Application/psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",stringsAsFactors = F)
  scale=seq(0,70,0.05)
  library(SKAT)
  skat_rs<-c()
  options(try.outFile = "PAH/PAH-CHD/skat.log.txt") 
  y<-c(rep(1,dim(map_case)[1]),rep(0,dim(map_control)[1]))
  obj<-SKAT_Null_Model(formula = y~1,Adjustment = T ,out_type = "D")
  for( gene in unique(case_D$GeneName)){
    # print(gene)
    gt_map<-c()
    weight<-c()
    if(length(which(het_lookup$V1==gene))>0){
      ### lof weight
      maxv<-max(lof_lookup[which(lof_lookup$V1==gene),3],(length(which(het_lookup[which(het_lookup$V1==gene),]>0))-1)*0.05,na.rm = T)
      snp_weight<-1-as.numeric(het_lookup[which(het_lookup$V1==gene),findInterval(maxv,scale)+1])
    }else{
      ### if no popscore for a gene, the weight=0.5
      snp_weight=0.5
    }
    ### cases
    for(AAchange in unique(c(case_D$AAChange[which(case_D$GeneName==gene)],control_D$AAChange[which(control_D$GeneName==gene)]))){
      arr_case<-rep(0,dim(map_case)[1])
      arr_control<-rep(0,dim(map_control)[1])
      cadd<-case_D$CADDphred[which(case_D$AAChange==AAchange)]
      if(length(cadd)==0){
        cadd<-control_D$CADDphred[which(control_D$AAChange==AAchange)]
      }
      probands<-case_D$ProbandName[which(case_D$AAChange==AAchange)]
      controls<-control_D$ProbandName[which(control_D$AAChange==AAchange)]
      if(length(probands)>0){
        index<-as.numeric(unlist(unlist(lapply(probands,FUN = function(x) unlist(map_case[which(map_case[,2]==x),1])))))
        arr_case[index]<-1;
      }
      #  controls<-control_D$ProbandName[which(control_D$GeneName==gene)]
      if(length(controls)>0){
        index<-as.numeric(unlist(lapply(controls,FUN = function(x) unlist(map_control[which(map_control[,2]==x),1]))))
        arr_control[index]<-1;
      }
      # print(cadd)
      if(length(cadd)==0){
        snp_weight<-1-as.numeric(het_lookup[which(het_lookup$V1==gene),findInterval(cadd,scale)+1])
      }
      weight<-c(weight,snp_weight)
      gt_map<-cbind(gt_map,c(arr_case,arr_control))
      
    }
    gt_map<-as.matrix(gt_map)
    #gt_map<-gt_map[,which(colSums(gt_map)!=0)]
    pvalue<-try(SKAT(gt_map,obj,weights = weight)$p.value,silent = T)
    if(is.numeric(pvalue)){
      skat_rs<-rbind(skat_rs,c(gene,pvalue,
                               sum(unlist(lapply(1:dim(map_case)[1],function(x) sum(gt_map[x,])))),
                               sum(unlist(lapply((dim(map_case)[1]+1):(dim(map_control)[1]+dim(map_case)[1]),
                                                 function(x) sum(gt_map[x,]))))))
    }else{
      print(paste("Error:",gene))
      skat_rs<-rbind(skat_rs,c(gene,'NA',
                               sum(unlist(lapply(1:dim(map_case)[1],function(x) sum(gt_map[x,])))),
                               sum(unlist(lapply((dim(map_case)[1]+1):(dim(map_control)[1]+dim(map_case)[1]),
                                                 function(x) sum(gt_map[x,]))))))
      
    }
  }
  skat_rs<-data.frame(skat_rs,stringsAsFactors = F)
  names(skat_rs)<-c("Gene","p-value","Nmutcase","Nmutcontrol")
  write.csv(skat_rs,file="PAH/PAH-CHD/PAH-CHD-skat_revel0.5.txt",row.names = F)
   pvalues<-sort(as.numeric(skat_rs$p.value[which(!is.na(as.numeric(as.character(skat_rs$p.value))))],decreasing = T))
 
  pdf("PAH/PAH-CHD/PAH-CHD-skat_revel0.5.pdf")
  x<--log10((length(pvalues):1)/length(pvalues))
  plot(x,-log10(pvalues))
  dev.off()
  
}

#### associated gene ### 
#fasso  <-  "PAH/DOC/HongJian/source/PAH_associated11-13.txt"
fasso<- "PAH/Result/Data/source/CHD.253GeneList.csv"
gene_asso  <-  read.csv(fasso,header = 1,stringsAsFactors = F,check.names = F)
gene_asso  <-  gene_asso$Gene
minor_gene_asso  <-  setdiff(gene_asso,c("BMPR2","TBX4"))
####

chd.list<-read.csv("PAH/Result/Data/source/PAH-CHD-list.csv",header = 1,stringsAsFactors = F,strip.white = T)
chd.list$FN<-chd.list$proband

#chd_solved_vars<-read.csv("PAH/PAH-CHD/Data/PAH-CHD.solved.0608.csv",header = 1,stringsAsFactors = F,check.names = F)
chd_solved_vars<-read.table("PAH/PAH-CHD/Data/PAH-CHD.solved.0809.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
#unsolved_denovo<-denovo[!which(denovo$ProbandName%in%solved_denovos),]
chd.pah.samples<-chd_solved_vars[,c("proband","Gene.refGene")]
names(chd.pah.samples)<-c("ID","GeneName")
solved_samples<-chd.pah.samples[which(chd.pah.samples$GeneName%in%gene_asso),]

pop_vcr  <-  read.table("PAH/DOC/HongJian/source/european.txt",header = F,check.names = F,stringsAsFactors = F,fill = T)
VCR_solved  <-  unique(solved_samples$ID[which(solved_samples$GeneName %in% gene_asso & solved_samples$ID%in%pop_vcr$V1)])
vcr_chd_european<-chd.list$proband[which(chd.list$proband%in%pop_vcr$V1 )]


pop <- read.csv("PAH/PAH_2017/Relatedness/Phenotype_FREEZE6_pop.csv",header=1,stringsAsFactors = F,check.names = F)
#euro_PAH <- pop$sampleID[intersect(which(pop$pop=="EUR"),grep(paste(c("_JM","_PPH","_SPH"),collapse = "|"),pop$sampleID))]
euro_xGen_PAH <- c(pop$sampleID[which(pop$pop=="EUR")]) #,pop_vcr$V1)
#all_PAH <- c(pheno_PAH$ID,pheno_SPH_PAH$ID)
xgen_CHD_european <- intersect(chd.list$proband,euro_xGen_PAH)
XGRN_solved<-xgen_CHD_european[which(xgen_CHD_european%in%solved_samples)]


control_IDs  <-  read.csv("PAH/DOC/HongJian/source/european-control.txt",header = F,check.names = F,stringsAsFactors = F)
control_IDs   <-   control_IDs$V1




### indel IGV check xGen ###  the indels failed IGV check
indel_check  <-  read.table("PAH/JointCalls/IGV.check.indel.txt",header = 1,check.names = F,stringsAsFactors = F)
indel_check  <-  indel_check[which(indel_check$Flag==0),]

case_false_indels  <-  read.table("PAH/Image/Remote.false.list",header = F,stringsAsFactors = F,check.names = F)
ctr_false_indels  <-  read.table("PAH/Image/control.false.indel.txt",header=F,stringsAsFactors = F,check.names = F)
exclude_set_cs  <-  unlist(lapply(1:dim(case_false_indels)[1],FUN = function(x) paste(case_false_indels$V1[x],case_false_indels$V2[x],case_false_indels$V3[x]) ))
exclude_set_ctr  <-  unlist(lapply(1:dim(ctr_false_indels)[1],FUN = function(x) paste(ctr_false_indels$V1[x],ctr_false_indels$V2[x],ctr_false_indels$V3[x]) ))


### start import variant #####

PAH_case_control <- read.table("PAH/Result/Data/source/PAH_control_VCX.anno.mapp.bed.hg19_multianno.txt",header = 1,sep = "\t",stringsAsFactors = F,check.names = F,comment.char = "")



control_dat <- PAH_case_control[which(PAH_case_control$ProbandName %in%  control_IDs),]
PAH_CHD_dat<-PAH_case_control[which(PAH_case_control$ProbandName%in%c(xgen_CHD_european,vcr_chd_european)),]


#VU_european_IPH_FPPH <- VU_cases[which(VU_cases$ProbandName %in% VU_euro),]
### end import variants ###

control_dat <- formatFreq_new(control_dat)
PAH_CHD_dat<-formatFreq_new(PAH_CHD_dat)
### start filter the variants ###

control_dat <- filter_allfreq_new(control_dat,0.0001,0.0001)
#VU_european_IPH_FPPH <- filter_allfreq(VU_european_IPH_FPPH)
PAH_CHD_dat<-filter_allfreq_new(PAH_CHD_dat,0.0001,0.0001)

#VU_european_IPH_FPPH <- VQSR_filter(VU_european_IPH_FPPH)
control_dat <-  VQSR_filter(control_dat)
PAH_CHD_dat <-  VQSR_filter(PAH_CHD_dat)


#control_dat <- COV_filter(control_dat)
#vcr_dat<-COV_filter(vcr_dat)
#xgen_dat <- COV_filter(xgen_dat)
#sph_dat<-  COV_filter(sph_dat)


control_dat <- control_dat[which(!paste(control_dat$ProbandName,control_dat$CHROM,control_dat$POS)  %in%  exclude_set_ctr),]
PAH_CHD_dat <- PAH_CHD_dat[which(!paste(PAH_CHD_dat$ProbandName,PAH_CHD_dat$CHROM,PAH_CHD_dat$POS)  %in%  exclude_set_cs),]
#vcr_dat <- vcr_dat[which(!paste(vcr_dat$ProbandName,vcr_dat$CHROM,vcr_dat$POS)  %in%  exclude_set_cs),]
#sph_dat <- sph_dat[which(!paste(sph_dat$ProbandName,sph_dat$CHROM,sph_dat$POS)  %in%  exclude_set_cs),]

control_dat<-control_dat[which(control_dat$ALT!="*"),]
PAH_CHD_dat <- PAH_CHD_dat[which(PAH_CHD_dat$ALT != "*"),]

PAH_CHD_dat<-format_consensus(PAH_CHD_dat)
control_dat<-format_consensus(control_dat)
#enrich_test(PAH_CHD_dat,control_dat,gene_asso,minor_gene_asso,solved_samples)


case_euro<-read.table("PAH/PAH-CHD/PAH-CHD.european.list",header = F,stringsAsFactors = F)
map_case<-cbind(1:dim(case_euro)[1],case_euro$V1);
control_euro<-read.table("PAH/PAH-CHD/PCGC_parent.european.list",header = F,stringsAsFactors = F)
map_control<-cbind(1:dim(control_euro)[1],control_euro$V1)

## only use LgD and CADD>20 missense variants
#case_D<-xGen_PAH[c(grep("^frame|^stop",xGen_PAH$VarClass),grep("splic",xGen_PAH$VarFunc),which(as.numeric(xGen_PAH$CADDphred)>20 &xGen_PAH$VarClass=="nonsynonymousSNV")),]
case_D<-PAH_CHD_dat[c(grep("^frame|^stop",PAH_CHD_dat$VarClass),grep("splic",PAH_CHD_dat$VarFunc),which(as.numeric(PAH_CHD_dat$REVEL)>0.5 &PAH_CHD_dat$VarClass=="nonsynonymous_SNV")),]
#VCR_PAH_D<-VCR_PAH[c(grep("^frame|^stop",VCR_PAH$VarClass),grep("splic",VCR_PAH$VarFunc),which(as.numeric(VCR_PAH$CADDphred)>20 &VCR_PAH$VarClass=="nonsynonymousSNV")),]
#case_D<-rbind(xGen_PAH_D[intersect(names(VCR_PAH_D),names(xGen_PAH_D))],VCR_PAH_D[,intersect(names(VCR_PAH_D),names(xGen_PAH_D))])
control_D<-control_dat[c(grep("^frame|^stop",control_dat$VarClass),grep("splic",control_dat$VarFunc),which(as.numeric(control_dat$REVEL)>0.5 &control_dat$VarClass=="nonsynonymous_SNV")),]
#xGen_PAH_D<-xGen_PAH[c(which(as.numeric(xGen_PAH$CADDphred)>20 &xGen_PAH$VarClass=="nonsynonymousSNV")),]
#VCR_PAH_D<-VCR_PAH[c(which(as.numeric(VCR_PAH$CADDphred)>20 &VCR_PAH$VarClass=="nonsynonymousSNV")),]
#case_D<-rbind(xGen_PAH_D[intersect(names(VCR_PAH_D),names(xGen_PAH_D))],VCR_PAH_D[,intersect(names(VCR_PAH_D),names(xGen_PAH_D))])
#control_D<-control_dat[c(which(as.numeric(control_dat$CADDphred)>20 &control_dat$VarClass=="nonsynonymousSNV")),]


skat(map_case,map_control,case_D,control_D)
### read gnomad file

write.table(control_dat[which(control_dat$ALT!="*"),],file = "PAH/PAH-CHD/PCGC_parent_rare.txt",sep = "\t",row.names = F,quote = F)

#write.table(xGen_PAH,file = "PAH/PAH-CHD.0.0001.european.txt",sep = "\t",row.names = F,quote = F)
case_D<-rbind(xGen_PAH_D[intersect(names(VCR_PAH),names(xGen_PAH))],VCR_PAH_D[,intersect(names(VCR_PAH),names(xGen_PAH))])

xx<-read.csv("~/server/PAH/PAH-CHD/PAH-CHD-skat_revel0.5.txt",header = 1,stringsAsFactors = F)
