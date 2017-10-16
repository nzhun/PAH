#setwd("~/server")
setwd("/home/local/ARCS/nz2274/")
source("Pipeline/NA_script/R/untils.R")

filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                    
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




fill_proband<-function(set){
  pedf<-"PAH/Result/Data/source/VCX_Control.ped"
  ped<-read.table(pedf,header=1,comment.char = "",check.names = F,stringsAsFactors = F,fill = T)
  set$proband=set$ID;
  for(i in 1:dim(set)[1]){
    index<-grep(paste(set$ID[i],"_",sep=""),ped$ID)
    if(length(index)>0){
      set$proband[i]=ped$ID[index]
    }else{
      if(length(grep("xgen",set$CAPTURE[i],ignore.case = T))>0){
        print(paste("failed",set$ID[i],set$CAPTURE[i]));
      }
    }
  }
  return (set)
}

addfamily<-function(dat,set){
  dat$age="";
  dat$type="";
  dat$gender=""
  dat$family="";
  for(f in unique(dat$proband)){
    info<-set[which(set$proband==f),]
    dat$age[which(dat$proband==f)]<-info$Age_dx
    dat$age[which(dat$proband==f)]<-info$TYPE
    dat$gender[which(dat$proband==f)]<-info$Gender
    dat$family[which(dat$proband==f)]<-info$FamType
  }
  return(dat)
}


process<-function(file,setfile,fasso,outputf){

  dat<-read.table(file,header = 1,stringsAsFactors = F,comment.char = "",check.names = F,sep="\t")
  set<-read.csv(setfile,header = 1,comment.char ="",check.names = F,stringsAsFactors = F,strip.white = T)
  dat<-formatFreq_new(dat)
  filter_dat<-filter_allfreq_local(dat,0.0001,0.0001)
  
  #set<-fill_proband(set)
  #write.csv(set,file = setfile,quote = F,row.names = F) 
 
  gene_asso  <-  read.table(fasso,header = F,stringsAsFactors = F,check.names = F,strip.white = T)
  gene_asso  <-  gene_asso[,1]
  filter_dat<-filter_dat[which(filter_dat$Gene.refGene%in%gene_asso &filter_dat$ALT!="*" & filter_dat$AF<0.01),] # &filter_dat$VQSLOD > -2.75 
  cohort<-c()
  if(length(grep("proband",names(set)))>0){cohort<-set$proband}else{cohort<-set[,1]}
  sub_dat<-filter_dat[which(filter_dat$proband%in%cohort),]
  if(length(grep("DP_ind",names(sub_dat))) >0){sub_dat<-sub_dat[which(as.numeric(sub_dat$DP_ind)>8| as.numeric(sub_dat$AD_ind)+as.numeric(sub_dat$RD_ind)>8),] }
  #if(length(grep("GQ",names(sub_dat))) >0){sub_dat<-sub_dat[which(as.numeric(sub_dat$GQ)>90),] }
  sub_dat<-addfamily(sub_dat,set)
  write.csv(sub_dat,file = outputf,row.names = F)
}


paper_check<-function(sub_dat){
  sub_dat$publish<-0;
  sub_dat$publish_id<-""
  paper<-read.csv("PAH/Result/Data/source/known_paper_adult.csv",header = 1,stringsAsFactors = F,strip.white = T)
  paper$has=0;
  for(i in 1:dim(sub_dat)[1]){
    proband<-sub_dat$proband[i]
    gene<-sub_dat$Gene.refGene[i]
    ids<-paper$ID[which(paper$Gene==gene)]
    for(id in ids){
      a<-grep(paste(id,sep=""),proband)
  
      if(length(a)>0){sub_dat$publish[i]=1; sub_dat$publish_id[i]=id; paper[which(paper$ID==id & paper$Gene==gene),"has"]<-1;  print(paste(id,proband,i));}
    }
  }
  write.csv(paper,file = "PAH/Result/Data/source/known_paper_adult.csv",row.names = F)
}
#file="PAH/Known/PAH.known.vcf.gz.2.txt"
#setfile="PAH/Result/Data/source/PAH_pediatric_adult_list.csv"
# fasso  <-  "PAH/DOC/HongJian/source/PAH_associated11-13.txt"
#outputf="PAH/PAH_known.variants.csv"  
args<-commandArgs(trailingOnly = T)
if(length(args)<4){
  stop("At least two arguments must be supplied (input file)\n\tinput format: the tab-delimited data file, cohort_ids,genesetfile,output_file ", call.=FALSE)
}else if(length(args)==3){
  args[4]="out.csv"
}
file=args[1]
setfile=args[2]
outputf=args[4]
fgene=args[3]
print("input format: the conerted_txt, cohort_ids,genesetfile,output_file")
process(file,setfile,fgene,outputf)