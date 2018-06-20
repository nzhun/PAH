# install.packages("~/Downloads/SKAT_1.3.2.1.tgz",repos =NULL)

filter_allfreq_EUR <- function(data,freq_avg,freq_max){
  
  data <- data[which(as.numeric(data$ExAC_NFE)< freq_max
                     
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     
  ),]
  
  return(data)
}

cmb_AF<-function(afs){
  #afs<-dat$BatchInfo
  AN_SPARK <- 2958
  AN_Nimble<-17460
  AN_Med<-6090
  AN_Baylor<-1200
  
   cmd_af<-unlist(lapply(afs,FUN = function(x){
      bts<-unique(unlist(strsplit(x,",")));
      dat_bat<-c()
      
      for(bt in bts){
        infos<-unlist(strsplit(bt,":"))
        #if(length(infos)<2){infos<-rep(0,4)}
        arr<-c(0,0,0,0)
        for(i in 1:length(infos)){
          if(infos[i]!="."&&infos[i]!="" && length(infos[i])>0){arr[i]<-infos[i]}
        }
        dat_bat<-rbind(dat_bat,arr)  
      }
      ac<-sum(round(as.numeric(dat_bat[,1])*as.numeric(dat_bat[,3])))
      an<-sum(as.numeric(dat_bat[,3]))
      if(length(grep("SPARK",x,ignore.case = T)) <1){an<-an+AN_SPARK;}
      if(length(grep("Baylor",x,ignore.case = T)) <1){an<-an+AN_Baylor;}
      if(length(grep("Med",x,ignore.case = T)) <1){an<-an+AN_Med;}
      if(length(grep("Nim",x,ignore.case = T)) <1){an<-an+AN_Nimble;}
      af<-ac/an;
      return(af)
  }))
   return(cmd_af)

}


filter_allfreq_local2 <- function(data,freq_avg,freq_max){
  
  data <- data[which(
    na.pass(as.numeric(data$ExAC_ALL)<= freq_avg)
   # &na.pass(as.numeric(data$ExAC_EUR)< freq_max)
    &na.pass(as.numeric(data$esp6500siv2_all)<= freq_max)
    &as.numeric(data$ExAC_NFE)<= freq_max
     & as.numeric(data$gnomAD_exome_ALL)<= freq_max
    &as.numeric(data$gnomAD_exome_NFE)<= freq_max

    & as.numeric(data$`1000g2015aug_all`) <= freq_max
    & as.numeric(data$`1000g2015aug_eur`) <= freq_max

  ),]
  # data <- data[which( as.numeric(data$AC)< 25
  #                      &as.numeric(data$AB)>0.2
  #  ),]
  
  return(data)
}



filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(
    na.pass(as.numeric(data$ExAC_ALL)<= freq_avg)
    &na.pass(as.numeric(data$ExAC_AMR)<= freq_max)
    &as.numeric(data$ExAC_AFR)<= freq_max
    &as.numeric(data$ExAC_NFE)<= freq_max
    &as.numeric(data$ExAC_FIN)<= freq_max
    &as.numeric(data$ExAC_SAS)<= freq_max
    &as.numeric(data$ExAC_EAS)<= freq_max
    &as.numeric(data$ExAC_OTH)<= freq_max
    #    &as.numeric(data$gnomAD_genome_EAS)<freq_max
    #    &as.numeric(data$gnomAD_genome_NFE)<freq_max
    #    &as.numeric(data$gnomAD_genome_FIN)<freq_max
    #   &as.numeric(data$gnomAD_genome_OTH)<freq_max
    #  &as.numeric(data$gnomAD_genome_ASJ)<freq_max
    # &as.numeric(data$gnomAD_genome_AMR)<freq_max
    #&as.numeric(data$gnomAD_genome_ALL)<freq_max
    &
      as.numeric(data$gnomAD_exome_ALL) <= freq_max
    &as.numeric(data$gnomAD_exome_EAS) <= freq_max
    &as.numeric(data$gnomAD_exome_NFE) <= freq_max
    &as.numeric(data$gnomAD_exome_FIN) <= freq_max
    &as.numeric(data$gnomAD_exome_OTH) <= freq_max
    &as.numeric(data$gnomAD_exome_ASJ) <= freq_max
    &as.numeric(data$gnomAD_exome_AMR) <= freq_max
    &as.numeric(data$gnomAD_exome_AFR) <= freq_max
    & as.numeric(data$`1000g2015aug_all`) <= freq_max
    & na.pass(as.numeric(data$esp6500siv2_all) <= freq_max)
    #      &as.numeric(data$ESPfreq)< freq2
    #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  # data <- data[which( as.numeric(data$AC)< 25
  #                      &as.numeric(data$AB)>0.2
  #  ),]
  
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
                       "VQSRTrancheSNP99.70to99.80",
                       #    "VQSRTrancheSNP99.60to99.70",#"VQSRTrancheSNP99.50to99.60",
                       #"VQSRTrancheINDEL99.50to99.60",
                       #  "VQSRTrancheINDEL99.60to99.70",                 
                       "VQSRTrancheINDEL99.70to99.80",  
                       "VQSRTrancheINDEL99.80to99.90", "VQSRTrancheINDEL99.90to100.00")
    
    data <- data[which(!data$FILTER %in%  filter_ignore |data$FILTER=="."),]
  }
  return(data)
}

map_filter<-function(data){
  index<-grep("mappability",names(data),ignore.case = T)
  if(length(index)>0){ 
    data <- data[which(as.numeric(data[,index])==1|(data[,index]==".")),]
  }
  if(length(grep("genomicSuperDups",names(data),ignore.case = T))>0){
    index<-grep("Score",data$genomicSuperDups)
    
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;|\\\\x3d|\\\\x3b"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return (data)
}




filter_case <- function(dat,total_case){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
  dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.7*max_AN )),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnomad_COV10!="." & (as.numeric(dat$gnomad_COV10)>0.9 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("90% cov10 and 75% cov15")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case)
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  
  
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.0001,0.0001) ## 
 # dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$BRAVO_AF)<0.0001 ),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$AF)<0.005 &  (as.numeric(dat$AC)/as.numeric(dat$AN)) <0.005),]
  cnt<-(c(dim(dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],dim(dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_case,
          dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case ))
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  
  gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  dat$GQ<-gq
  PAH_CHD_dat<-dat[-which(as.numeric(dat$GQ)<60),]  ## only for RGN
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("gq 30")
  print(cnt)
  filter_rs<-rbind(filter_rs,c("gq 70 ",cnt))
  PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  remove<- which( PAH_CHD_dat$AB<0.25 |PAH_CHD_dat$AD< 5|PAH_CHD_dat$DP_ind<15)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("AB 0.2 ",cnt))
  print("AB0.2")
  print(cnt)
  #PAH_CHD_dat<-PAH_CHD_dat[which(PAH_CHD_dat$FS<25),] ## only for RGN
#  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
#         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
#  filter_rs<-rbind(filter_rs,c("FS 25 ",cnt))
#  print("FS 25")
#  print(cnt)
  
  
  
 
  calt<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){ max(unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ) }) ) 
  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(calt)<4),] 
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("max alt ",cnt))
  
  
  # gq<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  is_indel<-unlist(lapply(1:dim(PAH_CHD_dat)[1],FUN = function(x){return(nchar(PAH_CHD_dat$REF[x])!=nchar(PAH_CHD_dat$ALT[x]))}))
  # PAH_CHD_dat$Is_indel<-is_indel
  # indel_fail<-intersect(which(PAH_CHD_dat$FS>200|(PAH_CHD_dat$ReadPosRankSum!="." & as.numeric(PAH_CHD_dat$ReadPosRankSum) < -20) |PAH_CHD_dat$SOR > 10 ),which(is_indel==T))
  snv_fail<-c(which(as.numeric(PAH_CHD_dat$MQ) < 40|
                                        PAH_CHD_dat$QD < 4
                                               |PAH_CHD_dat$AB  < 0.3
                                               |PAH_CHD_dat$AD_ind < 6
                                               | PAH_CHD_dat$DP_ind < 15
                                               |PAH_CHD_dat$FS > 20
                                               | PAH_CHD_dat$SOR < 4
                                                | PAH_CHD_dat$GQ < 90
                                               ), which(  as.numeric(PAH_CHD_dat$ReadPosRankSum) < -2| as.numeric(PAH_CHD_dat$MQRankSum ) < -2 )) #| PAH_CHD_dat$SOR > 3
  # # snv_gq_fail<-c() #intersect(which(gq<70),which(is_indel==F))
  remove<- intersect(unique(snv_fail),which(is_indel==F)) #(c(indel_fail,snv_fail,snv_gq_fail))
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
   cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
          dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
   filter_rs<-rbind(filter_rs,c("follow ",cnt))
  # 
  
  
  print(c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
          dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) )
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}


filter_gnomad<-function(dat){
  total_control=7509
  gfilter_rs<-c()
  #chd_dat <- read.csv("PAH/PAH-CHD/CHD/CHD.WES.253known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
  subdat<-dat[which(dat$FILTER%in%c("PASS") & dat$ExonicFunc.refGene!="unknown" & as.numeric(dat$AC_NFE)>0),] #,"RF"
  remove<- which(subdat$`#CHROM`=="X" & subdat$FROM>136251321 & subdat$FROM<136251327)
  if(length(remove)>0){
    subdat<-subdat[-remove,]
  }
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))
#  if(subdat$`#CHROM`!="X"){
    max_gan<-max(subdat$AN_NFE)
    index<-which((subdat$`#CHROM`!="X" & subdat$AN_NFE/max_gan >0.9)|subdat$`#CHROM`=="X")
    subdat<-subdat[index,]
    subdat<-subdat[which(subdat$AN_NFE/max_gan >0.75),]
    
 # }
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))

  subdat<-map_filter(subdat)
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  
  subdat <-subdat[grep("^MUC|^HLA",subdat$Gene.refGene,invert = T),]
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MUC|HLA",cnt))
  
  # subdat<-map_filter(subdat)
  # cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  # gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  subdat <- formatFreq_new(subdat)
  
  subdat <- filter_allfreq_local2(subdat,0.0001,0.0001)  ## Max gnomad exome <10^-4 ## Max gnomad exome <10^-4
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MAX pop AF",cnt))
  
  subdat<-subdat[which(subdat$BRAVO_AF=="."| as.numeric(subdat$BRAVO_AF)<0.0001),]
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("BRAVO_AF",cnt))
  
  subdat <- subdat[which(as.numeric(subdat$AF)<0.005 & as.numeric(subdat$AF_NFE)<0.005  ),] 
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("local AF",cnt))
  

  subdat <- subdat[which(as.numeric(subdat$VQSLOD) > -3.5 &   as.numeric(subdat$MQ) > 40  & as.numeric(subdat$FS) < 20 ),] 
  
 # subdat<-subdat[grep("0|0|0|0|0",subdat$AB_HIST_ALT,fixed = T,invert = F),]
#  subdat<-subdat[grep("0|0|0|0|0",subdat$AB_HIST_ALL,fixed = T,invert = F),]
  remove<-which(as.numeric(subdat$BaseQRankSum)< -2 |
                  as.numeric(subdat$ClippingRankSum) < -1 |
                  as.numeric(subdat$InbreedingCoeff) < -0.1 |
                  as.numeric(subdat$FS) > 20 |
                  as.numeric(subdat$QD) < 3 |
                  as.numeric(subdat$MQRankSum) < -2 |
                  as.numeric(subdat$ReadPosRankSum) < -2 |
                  as.numeric(subdat$SOR) > 3 |
                  as.numeric(subdat$DP/subdat$AN) < 15|
                  as.numeric(subdat$MQ) < 50
  )
  
  
  index<-intersect(grep('non',subdat$ExonicFunc.refGene,invert = T),remove)
  if(length(index)>0){
    subdat<-subdat[-index,]
  }
#  subdat <- subdat[which(as.numeric(subdat$VQSLOD) > -10  ),] 
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("VQSLOD",cnt))
  return (list(data=subdat,filter=gfilter_rs))
}

ind_var<-function(dat){
    is_homo<-unlist(lapply(dat$Genotype,FUN = function(x){
    alleles=unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ; return(alleles[1]==alleles[2]) }))
    ad<-unlist(lapply(dat$Genotype,FUN = function(x){
    alleles=unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ; 
    alleles[which(alleles==".")]=0;
    alleles<-as.numeric(alleles)
    dp4=unlist(strsplit(x,":"))[2];  
    alts=as.numeric(unlist(strsplit(dp4, ",")))
    return (alts[alleles[2]+1])
    
  }))
  dp<-unlist(lapply(dat$Genotype,FUN = function(x){
    alleles=unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ; 
    alleles[which(alleles==".")]=0;
    alleles<-as.numeric(alleles)
    dp4=unlist(strsplit(x,":"))[2];  
    alts=as.numeric(unlist(strsplit(dp4, ",")))
    return (alts[alleles[2]+1]+alts[alleles[1]+1])
    
  }))
  ab<-unlist(lapply(dat$Genotype,FUN = function(x){
    alleles=unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ; 
    alleles[which(alleles==".")]=0;
    alleles<-as.numeric(alleles)
    dp4=unlist(strsplit(x,":"))[2];  
    alts=as.numeric(unlist(strsplit(dp4, ",")))
    ab=alts[alleles[2]+1]/sum(alts);
    return(ab)
  }))  ## For RGN
  dat$AD<-ad;
  dat$AB<-ab;
  dat$DP_ind<-dp;
  dat$Is_homo<-is_homo;
  return(dat)
}


filter_control <- function(dat,total_control){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c();
  
#  dat$CMB_AF<-cmb_AF(dat$BatchInfo)
  dat<- dat[which(dat$CMB_AF<0.01),]
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("cohort 0.01 ",cnt))
  print(cnt)
  
#  dat<-ind_var(dat)
 # dat$calt<-unlist(lapply(dat$Genotype,FUN = function(x){length(unlist(strsplit(unlist(strsplit(x,split = ":"))[2],split = "," )))-1 }))
  remove<- which( (dat$AB<0.15 & dat$Is_homo==F)|(dat$AD<4 & dat$AB<0.25 ) |
                    dat$calt>3
  ) #ctr_filter_dat$DP_ind<10|ctr_filter_dat$DP_ind >150 ||ctr_filter_dat$AD<4
  if(length(remove)>0){
    dat<- dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("AB 0.15 ",cnt))
  print(cnt)

 # gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
 # dat$GQ_Ind<-gq
  dat<-dat[which(as.numeric(dat$GQ_Ind)>30),]  ## only for RGN
  
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("GQ 30 ",cnt))
  print(cnt)
  
  dat<-dat[which( dat$GnomAD_Genome_cov10!="." & (as.numeric(dat$GnomAD_Genome_cov10)>0.9) ),] #  &  as.numeric(dat$GnomAD_Genome_cov15)>0.75
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  print("90% cov10 ")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  
  
  
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_control,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_control ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.0001,0.0001) ## 
  dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  bravo_index<-grep("BRAVO_.*AF",names(dat),ignore.case = T)
  dat<-dat[which( as.numeric(dat[,bravo_index])<0.0001 ),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  

  ctr_filter_dat<-dat
 
  is_indel<-unlist(lapply(1:dim(ctr_filter_dat)[1],FUN = function(x){return(nchar(ctr_filter_dat$REF[x])!=nchar(ctr_filter_dat$ALT[x]))}))
  ctr_filter_dat$Is_indel<-is_indel
  remove<- which((as.numeric(ctr_filter_dat$DP_ind)<15 |as.numeric(ctr_filter_dat$AB)<0.25 | as.numeric(ctr_filter_dat$AD)<4) & ctr_filter_dat$Is_indel==T)
  if(length(remove)>0){
    ctr_filter_dat<-ctr_filter_dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("indel-specific ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  return(list(data=ctr_filter_dat,filter=filter_rs))
}



filter_control_merge <- function(dat,total_control){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c();
  
  #  dat$CMB_AF<-cmb_AF(dat$BatchInfo)
  dat<- dat[which(dat$CMB_AF<0.01),]
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("cohort 0.01 ",cnt))
  print(cnt)
  
  #  dat<-ind_var(dat)
  # dat$calt<-unlist(lapply(dat$Genotype,FUN = function(x){length(unlist(strsplit(unlist(strsplit(x,split = ":"))[2],split = "," )))-1 }))
  # remove<- which( (dat$AB<0.25 & dat$Is_homo==F)|(dat$AD<4 & dat$AB<0.25 ) |
  #                   dat$calt>2
  # ) #ctr_filter_dat$DP_ind<10|ctr_filter_dat$DP_ind >150 ||ctr_filter_dat$AD<4
  

  
  dat<-dat[which( dat$GnomAD_Genome_cov10!="." & (as.numeric(dat$GnomAD_Genome_cov10)>0.9) ),] #  &  as.numeric(dat$GnomAD_Genome_cov15)>0.75
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  print("90% cov10 ")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  
  
  
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_control,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_control ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.0001,0.0001) ## 
  # dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  bravo_index<-grep("BRAVO_.*AF",names(dat),ignore.case = T)
  dat<-dat[which( as.numeric(dat[,bravo_index])<0.0001 ),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  remove<- which( dat$AB<0.2 | dat$DP_ind< 6 |dat$AD<3 |  dat$calt>2 )
  if(length(remove)>0){
    dat<- dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("AB 0.15 ",cnt))
  print(cnt)
  
  # gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  # dat$GQ_Ind<-gq
  dat<-dat[which(as.numeric(dat$GQ_Ind)>30),]  ## only for RGN
  
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("GQ 60 ",cnt))
  print(cnt)
  
  ctr_filter_dat<-dat
  
  # is_indel<-unlist(lapply(1:dim(ctr_filter_dat)[1],FUN = function(x){return(nchar(ctr_filter_dat$REF[x])!=nchar(ctr_filter_dat$ALT[x]))}))
  # ctr_filter_dat$Is_indel<-is_indel
   remove<- which((as.numeric(dat$GQ_Ind) <70 | dat$DP_ind< 10|   as.numeric(ctr_filter_dat$AB)< 0.25|  as.numeric(ctr_filter_dat$AD)< 6) & ctr_filter_dat$Is_indel==T)
   if(length(remove)>0){
     ctr_filter_dat<-ctr_filter_dat[-remove,] ## only for RGN
   }
  # cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
  #        dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  # filter_rs<-rbind(filter_rs,c("indel-specific ",cnt))
  # print("Bravo 10-4")
  # print(cnt)
  
  return(list(data=ctr_filter_dat,filter=filter_rs))
}



filter_case_inH <- function(dat,total_case){
  
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
  dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.6*max_AN )),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnomad_COV10!="." & (as.numeric(dat$gnomad_COV10)>0.9 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("90% cov10 and 75% cov15")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case)
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))

  
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.0001,0.0001) ## 
  dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$BRAVO_AF)<0.0001 ),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$AF)<0.01 &  (as.numeric(dat$AC)/as.numeric(dat$AN)) <0.01),]
  cnt<-(c(dim(dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],dim(dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_case,
          dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case ))
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  
  gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  dat$GQ<-gq
  PAH_CHD_dat<-dat[-which(as.numeric(dat$GQ)<60),]  ## only for RGN
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("gq 30")
  print(cnt)
  filter_rs<-rbind(filter_rs,c("gq 70 ",cnt))
  PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  remove<- which( PAH_CHD_dat$AB<0.25 )  #PAH_CHD_dat$DP_ind<10 |
  remove<-c(remove,which( (PAH_CHD_dat$AD<5 & PAH_CHD_dat$AB<0.3)) ) #
  remove<-c(remove,which((PAH_CHD_dat$AD<4)) ) #
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("AB 0.2 ",cnt))
  print(cnt)
  PAH_CHD_dat<-PAH_CHD_dat[which(PAH_CHD_dat$FS<40),] ## only for RGN
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("FS 25 ",cnt))
  print("FS 25")
  print(cnt)
  
  
  
#  PAH_CHD_dat<-PAH_CHD_dat[which(PAH_CHD_dat$Gene.refGene!="HCFC1"),] ## only for RGN
#  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
 #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
#  filter_rs<-rbind(filter_rs,c("HCFC1 ",cnt))
  #  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(PAH_CHD_dat$QD)>1.5),] 
  #  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
  #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # filter_rs<-rbind(filter_rs,c("QD 2.5 ",cnt))
  
  calt<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){ max(unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ) }) ) 
  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(calt)<3),] 
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("max alt ",cnt))
  
  
  # gq<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
   is_indel<-unlist(lapply(1:dim(PAH_CHD_dat)[1],FUN = function(x){return(nchar(PAH_CHD_dat$REF[x])!=nchar(PAH_CHD_dat$ALT[x]))}))
  # PAH_CHD_dat$Is_indel<-is_indel
  # indel_fail<-intersect(which(PAH_CHD_dat$FS>200|(PAH_CHD_dat$ReadPosRankSum!="." & as.numeric(PAH_CHD_dat$ReadPosRankSum) < -20) |PAH_CHD_dat$SOR > 10 ),which(is_indel==T))
   snv_fail<-intersect(which(is_indel==F),which(as.numeric(PAH_CHD_dat$MQ)< 40|
                                                  (PAH_CHD_dat$ReadPosRankSum!="." &( as.numeric(PAH_CHD_dat$ReadPosRankSum) < -5 |
                                                                                        as.numeric(PAH_CHD_dat$MQRankSum )< -5 )) |PAH_CHD_dat$AB<0.25 | PAH_CHD_dat$DP_ind<15)) #| PAH_CHD_dat$SOR > 3
  # # snv_gq_fail<-c() #intersect(which(gq<70),which(is_indel==F))
   remove<- snv_fail #(c(indel_fail,snv_fail,snv_gq_fail))
   if(length(remove)>0){
     PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
   }
  # cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
  #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # filter_rs<-rbind(filter_rs,c("follow ",cnt))
  # 
  
  
  print(c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
          dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) )
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}





filter_case_inH_merge <- function(dat,total_case){
  
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
 #dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.9*max_AN )),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnomad_COV10!="." &dat$gnomad_COV15!="." & (as.numeric(dat$gnomad_COV10)>0.9 )  & (as.numeric(dat$gnomad_COV15)>0.85 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("90% cov10 and 75% cov15")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case)
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  
  
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.0001,0.0001) ## 
 # dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$BRAVO_AF)<0.0001 ),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$AF)<0.01 &  (as.numeric(dat$AC)/as.numeric(dat$AN)) <0.01),]
  cnt<-(c(dim(dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],dim(dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_case,
          dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case ))
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  dat<-dat[which(as.numeric(dat$MQ) >40),]
 # gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  #dat$GQ<-gq
  rms<-which(as.numeric(dat$GQ_Ind)<60)
  if(length(rms)>0){
    PAH_CHD_dat<-dat[-rms,]  ## only for RGN
  }else{
    PAH_CHD_dat<-dat
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("gq 30")
  print(cnt)
  filter_rs<-rbind(filter_rs,c("gq 70 ",cnt))
  PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  remove<- which( PAH_CHD_dat$AB<0.25 |PAH_CHD_dat$DP_ind < 9  | PAH_CHD_dat$FS > 45 )  #PAH_CHD_dat$DP_ind<10 |
  #remove<-c(remove,which( (PAH_CHD_dat$AD<5 & PAH_CHD_dat$AB<0.3)) ) #
  remove<-c(remove,which((PAH_CHD_dat$AD<4)) ) #
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  remove<- grep("|",PAH_CHD_dat$Genotype,fixed = T)
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,]
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("AB 0.2 ",cnt))
  print(cnt)

  calt<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){ max(unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ) }) ) 
  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(calt)<3),] 
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("max alt ",cnt))
  
  
  # gq<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  is_indel<-unlist(lapply(1:dim(PAH_CHD_dat)[1],FUN = function(x){return(nchar(PAH_CHD_dat$REF[x])!=nchar(PAH_CHD_dat$ALT[x]))}))
  # PAH_CHD_dat$Is_indel<-is_indel
  # indel_fail<-intersect(which(PAH_CHD_dat$FS>200|(PAH_CHD_dat$ReadPosRankSum!="." & as.numeric(PAH_CHD_dat$ReadPosRankSum) < -20) |PAH_CHD_dat$SOR > 10 ),which(is_indel==T))
 # remove<-which(as.numeric(PAH_CHD_dat$MQ)< 40|(PAH_CHD_dat$ReadPosRankSum!="." &( as.numeric(PAH_CHD_dat$ReadPosRankSum) < -5 |as.numeric(PAH_CHD_dat$MQRankSum )< -5 )) )
  snps<- which(PAH_CHD_dat$Is_indel==F)
  remove<-intersect(snps,which(as.numeric(PAH_CHD_dat$FS)>20|as.numeric(PAH_CHD_dat$MQ)< 40|as.numeric(PAH_CHD_dat$MQRankSum)< -2|as.numeric(PAH_CHD_dat$ReadPosRankSum)< -2 |as.numeric(PAH_CHD_dat$GQ_Ind) <90 
                               |PAH_CHD_dat$SOR > 10| PAH_CHD_dat$DP_ind >250 | PAH_CHD_dat$DP_ind< 15| PAH_CHD_dat$AD <8 |as.numeric(PAH_CHD_dat$AB )<0.3 ) )
#  remove<-intersect(which(PAH_CHD_dat$Is_indel==F),remove) #| PAH_CHD_dat$SOR > 3
 
  #|PAH_CHD_dat$AB<0.25 | PAH_CHD_dat$DP_ind<15
  # # snv_gq_fail<-c() #intersect(which(gq<70),which(is_indel==F))
  #remove<- snv_fail #(c(indel_fail,snv_fail,snv_gq_fail))
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
  
  # cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
  #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # filter_rs<-rbind(filter_rs,c("follow ",cnt))
  # 
  
  
  print(c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
          dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) )
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}


filter_control_0.001 <- function(dat,total_control){
  filter_rs<-c();
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  #  dat$CMB_AF<-cmb_AF(dat$BatchInfo)
  dat<- dat[which(dat$CMB_AF<0.01),]
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("cohort 0.01 ",cnt))
  print(cnt)
  
  #  dat<-ind_var(dat)
  # dat$calt<-unlist(lapply(dat$Genotype,FUN = function(x){length(unlist(strsplit(unlist(strsplit(x,split = ":"))[2],split = "," )))-1 }))
  remove<- which( dat$AB<0.3 ) #ctr_filter_dat$DP_ind<10|ctr_filter_dat$DP_ind >150 ||ctr_filter_dat$AD<4
  remove<- which(   dat$calt>2   ) 
  remove<- c(remove,which( dat$AD <4  ))
  
  remove<- c(remove,which( dat$DP_Ind <20  ))
  if(length(remove)>0){
    dat<- dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("AB 0.15 ",cnt))
  print(cnt)
  
  # gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  # dat$GQ_Ind<-gq
  dat<-dat[which(as.numeric(dat$GQ_Ind)>90),]  ## only for RGN
  
  cnt<-c(dim( dat)[1]/total_control,dim(dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("GQ 30 ",cnt))
  print(cnt)
  
  dat<-dat[which( dat$GnomAD_Genome_cov10!="." & (as.numeric(dat$GnomAD_Genome_cov10)>0.9) ),] #  &  as.numeric(dat$GnomAD_Genome_cov15)>0.75
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  print("90% cov10 ")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  
  
  
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_control,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_control ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.001,0.001) ## 
  dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.001 & as.numeric(dat$gnomAD_genome_NFE)<0.001),] ## 
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  bravo_index<-grep("BRAVO_.*AF",names(dat),ignore.case = T)
  dat<-dat[which( as.numeric(dat[,bravo_index])<0.001 ),]
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  
  ctr_filter_dat<-dat
  
  is_indel<-unlist(lapply(1:dim(ctr_filter_dat)[1],FUN = function(x){return(nchar(ctr_filter_dat$REF[x])!=nchar(ctr_filter_dat$ALT[x]))}))
  ctr_filter_dat$Is_indel<-is_indel
  remove<- which((as.numeric(ctr_filter_dat$DP_ind)<20 |as.numeric(ctr_filter_dat$AB)<0.3 | as.numeric(ctr_filter_dat$AD)<6 |as.numeric(ctr_filter_dat$VQSLOD)< -100) & ctr_filter_dat$Is_indel==T)
  remove<- c(remove,which((as.numeric(ctr_filter_dat$VQSLOD)< -10) & ctr_filter_dat$Is_indel==T))
  
  if(length(remove)>0){
 
       ctr_filter_dat<-ctr_filter_dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( dat)[1]/total_control,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_control,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_control) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  return(list(data=ctr_filter_dat,filter=filter_rs))
}




filter_case_inH_0.001 <- function(dat,total_case){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
  
  dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.7*max_AN )),]
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnomad_COV10!="." & (as.numeric(dat$gnomad_COV10)>0.9 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("90% cov10 and 75% cov15")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case)
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  
  
  dat <- formatFreq_new(dat)
  dat <- filter_allfreq_local2(dat,0.001,0.001) ## 
  dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.001 & as.numeric(dat$gnomAD_genome_NFE)<0.001),] ## 
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$BRAVO_AF)<0.001 ),]
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  print("Bravo 10-4")
  print(cnt)
  
  dat<-dat[which( as.numeric(dat$AF)<0.01 &  (as.numeric(dat$AC)/as.numeric(dat$AN)) <0.01),]
  cnt<-(c(dim(dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],dim(dat[grep("^synony",dat$ExonicFunc.refGene),])[1]/total_case,
          dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case ))
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  
  gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  dat$GQ<-gq
  PAH_CHD_dat<-dat[-which(as.numeric(dat$GQ)<90),]  ## only for RGN
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("gq 30")
  print(cnt)
  filter_rs<-rbind(filter_rs,c("gq 70 ",cnt))
  PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  print(cnt)
  
 
  calt<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){ max(unlist(strsplit(unlist(strsplit(x,":"))[1],"/")) ) }) ) 
  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(calt)<3),] 
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("max alt ",cnt))
  
  
  # gq<-unlist(lapply(PAH_CHD_dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  is_indel<-unlist(lapply(1:dim(PAH_CHD_dat)[1],FUN = function(x){return(nchar(PAH_CHD_dat$REF[x])!=nchar(PAH_CHD_dat$ALT[x]))}))
  PAH_CHD_dat$Is_indel<-is_indel
  indel_fail<-intersect(which(PAH_CHD_dat$AB<0.3 |PAH_CHD_dat$DP_ind<20 |PAH_CHD_dat$GQ_Ind<90 |PAH_CHD_dat$AD<5 | PAH_CHD_dat$calt>2),which(PAH_CHD_dat$Is_indel==T))
  
  remove<-c(indel_fail)
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("follow ",cnt))
  
  
  
  print(c(dim( PAH_CHD_dat)[1]/total_case,dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1],
          dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) )
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}



format_geneName<-function(data){
  key="Gene.refGene"
  id<-grep(key,names(data))[1]
  index<-grep("-|;",data[,id])
  gene2<-unlist(lapply(data[index,id],FUN = function(x){
    arr<-unlist(strsplit(x,split = "-|;",fixed = F)); 
    if(length(arr)>2){
      if(nchar(arr[1])>nchar(arr[length(arr)])) {
        return(arr[1])
      }else{
        return(arr[length(arr)])
      }
    }
    if(nchar(arr[1]) >3 &nchar(arr[2])>3){
      return(arr[1])
    }else{
      return(x)
    }
  }))
  
  data[index,id]<-gene2
  #data$Gene.refGene[which(data$Gene.refGene=="LOC101928841")]<-"ADPRHL1"
  data<-data[grep("^MUC|^HLA|LOC101928841",data$Gene.refGene,invert = T),]
  return(data)
}







skat<-function(map_case,map_control,case_D,control_D){
  ## load Lookup table, get the weight
 # het_lookup<-read.table("Application/psap/lookups/full.het.pCADD.gencodeV19.allsites.txt.gz",stringsAsFactors = F)
#  lof_lookup<-read.table("Application/psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz",stringsAsFactors = F)
 # scale=seq(0,70,0.05)
  library(SKAT)
  skat_rs<-c()
 # options(try.outFile = paste(outf,".log",sep="")) 
  y<-c(rep(1,dim(map_case)[1]),rep(0,dim(map_control)[1]))
  obj<-SKAT_Null_Model(formula = y~1,Adjustment = T ,out_type = "D",n.Resampling=10000)
  key_gene=grep("GeneName|Gene.refGene",names(case_D))
  ### only works on rare missense variants and LGD
  for( gene in unique(case_D[,key_gene])){
   #  print(gene)
    gt_map<-c()
    weight<-c() ##  setup weight
    snp_weight<-1;
    ### cases set up matrix
    key_case=names(case_D)[grep("AAChange",names(case_D))]
    key_control=names(control_D)[grep("AAChange",names(control_D))]
    
    key_gene_case=names(case_D)[grep("GeneName|Gene.refGene",names(case_D))]
    
    key_gene_control=names(control_D)[grep("GeneName|Gene.refGene",names(control_D))]
    weight<-c()
    variants<-unique(rbind(case_D[which(case_D[,key_gene_case]==gene),c("FROM","REF","ALT")],
                           control_D[which(control_D[,key_gene_control]==gene),c("FROM","REF","ALT")]))
    
    for(t in 1:dim(variants)[1]){
      v=variants[t,]
      arr_case<-rep(0,dim(map_case)[1])
      arr_control<-rep(0,dim(map_control)[1])
      # cadd<-case_D$CADDphred[which(case_D[,key_case]==AAchange)]
      # if(length(cadd)==0){
      #   cadd<-control_D$CADDphred[which(control_D[,key_control]==AAchange)]
      # }
     
      probands<-case_D$ProbandName[which(case_D$FROM==v$FROM  &case_D$REF==v$REF & case_D$ALT==v$ALT)]
      controls<-c()
      if(length(grep("ProbandName",names(control_D))) >0 ){
        controls<-control_D$ProbandName[which(control_D$FROM==v$FROM & control_D$REF==v$REF & control_D$ALT==v$ALT)]
      }else{
        ### set up gnomad simulation
        ids<-which(control_D$FROM==v$FROM & control_D$REF==v$REF & control_D$ALT==v$ALT)
        if(length(ids) >0){
          ids<-ids[1]
         # if(length(ids)>1){print(AAchange)}
          # gts<-rbinom(dim(map_control)[1],2,as.numeric(control_D$AF_NFE[ids]))
          # controls<-map_control[which(gts >0),2];
          # if(length(controls)==0){controls<-control}
        #  print("Y")
       #  print(c(ids,dim(map_control)[1],control_D$AC_NFE[ids],control_D$Hom_NFE[ids]))
          gts<-rep(0,(dim(map_control)[1]-control_D$AC_NFE[ids]+control_D$Hom_NFE[ids]))
          if(control_D$AC_NFE[ids]-2*control_D$Hom_NFE[ids] >0){
            gts<-c(gts,rep(1,(control_D$AC_NFE[ids]-2*control_D$Hom_NFE[ids])))
          }
          if(control_D$Hom_NFE[ids]>0){
            gts<-c(gts,rep(2,control_D$Hom_NFE[ids]))
          }
          if(length(gts)!=dim(map_control)[1]){print(paste("Error:",gene));next;}
          gts<-sample(gts,length(gts),replace = F)
          ctr_hit<-which(gts > 0)
          if(length(ctr_hit)>0){
            controls<-map_control[ctr_hit,2]
          }else{
            print("NO calls",control_D[ids,])
          }
          if(control_D$REVEL[ids]!="."){
            snp_weight<-as.numeric(control_D$REVEL[ids])
          }else{
            snp_weight<-1
          }
        }
      }
      ## simulate from gnomad 
      if(length(probands)>0){
        id2<-which(case_D$FROM==v$FROM  &case_D$REF==v$REF & case_D$ALT==v$ALT)[1]
        if(case_D$REVEL[id2]!="."){
          snp_weight<-as.numeric(case_D$REVEL[id2])
        }else{
          snp_weight<-1
        }
        index<-as.numeric(unlist(unlist(lapply(probands,FUN = function(x) unlist(map_case[which(map_case[,2]==x),1])))))
        arr_case[index]<-1;
      }
      #  controls<-control_D$ProbandName[which(control_D$GeneName==gene)]
      if(length(controls)>0 ){ ### need randomly generate for gnomad
        index<-as.numeric(unlist(lapply(controls,FUN = function(x) unlist(map_control[which(map_control[,2]==x),1]))))
        arr_control[index]<-1;
      }
     # weights<-c()
      # print(cadd)
      # if(length(cadd)==0){
      #   snp_weight<-1-as.numeric(het_lookup[which(het_lookup$V1==gene),findInterval(cadd,scale)+1])
      # }
      weight<-c(weight,snp_weight)
      gt_map<-cbind(gt_map,c(arr_case,arr_control))
      
    }
    gt_map<-as.matrix(gt_map)
    #gt_map<-gt_map[,which(colSums(gt_map)!=0)]
  #  weight<-rep(1,dim(gt_map)[2])
    skat<-try(SKAT(gt_map,obj,weights = weight,method = "SKATO"),silent = T)  ## need add permutation
    re<-try(Get_Resampling_Pvalue(skat),silent = T)
    repv<-NA
    if(length(re)>1){repv<-re$p.value}
    if(length(skat)>0){ #is.numeric(pvalue)
      skat_rs<-rbind(skat_rs,c(gene,skat$p.value,repv,
                               sum(unlist(lapply(1:dim(map_case)[1],function(x) sum(gt_map[x,])))),
                               sum(unlist(lapply((dim(map_case)[1]+1):(dim(map_control)[1]+dim(map_case)[1]),
                                                 function(x) sum(gt_map[x,]))))))
    }else{
      print(paste("Error:",gene))
      skat_rs<-rbind(skat_rs,c(gene,'NA','NA',
                               sum(unlist(lapply(1:dim(map_case)[1],function(x) sum(gt_map[x,])))),
                               sum(unlist(lapply((dim(map_case)[1]+1):(dim(map_control)[1]+dim(map_case)[1]),
                                                 function(x) sum(gt_map[x,]))))))
      
    }
  }
  skat_rs<-data.frame(skat_rs,stringsAsFactors = F)
  names(skat_rs)<-c("Gene","p-value","emp-pvalue","Nmutcase","Nmutcontrol")
  return(skat_rs)

}


vat<-function(dat,ctr_dat,total_case,total_control){
  
  #z(t)=sum(1..m)(sum(1..n) di*cij*(pij-mean(pi)))/sqrt(sum(1..m)(sum(1..n) (di*cij)^2))
  #m is the number of variants on gene i and n is the number of samples, cij is the variant allele copy in i, di is has or not, pij is the phenotype of j
  mph<-total_case/2+total_control/2
  for (g in unique(dat$Gene.refGene)){
   subdat<-dat[which(dat$Gene.refGene==g),]
   ctr_sub<-ctr_dat[which(ctr_dat$Gene.refGene==g),]
   thres=0.5
   subdat$REVEL[grep('non',subdat$ExonicFunc.refGene,invert = T)]<-1
   subdat$REVEL[which(subdat$REVEL==".")]<-0
   subdat$Affected=2
   ctr_sub$REVEL[grep('non',ctr_sub$ExonicFunc.refGene,invert = T)]<-1
   ctr_sub$REVEL[which(ctr_sub$REVEL==".")]<-0
   ctr_sub$Affected=1
   sub_dat$copy<-1;
   sub_dat$copy[which(sub_dat$Is_homo==T)]<-2
   ctr_sub$copy<-1;
   ctr_sub$copy[which(sub_dat$Is_homo==T)]<-2
   
   z<-c()
   for(i in c(1:dim(subdat)[1])){
     
   }
  }
}



filter_gnomad_merge<-function(dat){
  total_control=7509
  gfilter_rs<-c()
  #chd_dat <- read.csv("PAH/PAH-CHD/CHD/CHD.WES.253known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
  subdat<-dat[which(dat$FILTER%in%c("PASS") & dat$ExonicFunc.refGene!="unknown" & as.numeric(dat$AC_NFE)>0),] #,"RF"
  remove<- which(subdat$`#CHROM`=="X" & subdat$FROM>136251321 & subdat$FROM<136251327)
  if(length(remove)>0){
    subdat<-subdat[-remove,]
  }
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))
  #  if(subdat$`#CHROM`!="X"){
  max_gan<-max(subdat$AN_NFE)
#  index<-which((subdat$`#CHROM`!="X" & subdat$AN_NFE/max_gan >0.9)|subdat$`#CHROM`=="X")
  index<-which(subdat$AN_NFE/max_gan >0.9)
  
   subdat<-subdat[index,]
#  subdat<-subdat[which(subdat$AN_NFE/max_gan >0.9),]
  
  # }
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))
  
  subdat<-map_filter(subdat)
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  
  subdat <-subdat[grep("^MUC|^HLA",subdat$Gene.refGene,invert = T),]
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MUC|HLA",cnt))
  
  # subdat<-map_filter(subdat)
  # cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  # gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  subdat <- formatFreq_new(subdat)
  
  subdat <- filter_allfreq_local2(subdat,0.0001,0.0001)  ## Max gnomad exome <10^-4 ## Max gnomad exome <10^-4
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("MAX pop AF",cnt))
  
  subdat<-subdat[which(subdat$BRAVO_AF=="."| as.numeric(subdat$BRAVO_AF)<0.0001),]
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("BRAVO_AF",cnt))
  
  subdat <- subdat[which(as.numeric(subdat$AF)<0.01  ),] 
  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
  gfilter_rs<-rbind(gfilter_rs,c("local AF",cnt))
  
  
  subdat <- subdat[which(as.numeric(subdat$VQSLOD) > -5.5  & as.numeric(subdat$MQ >20) & as.numeric(subdat$FS) < 35 ),] 
  
  index<-intersect(grep('non',subdat$ExonicFunc.refGene,invert = T),which(as.numeric(subdat$FS) >20|as.numeric(subdat$MQRankSum)< -2 |
                                                                            as.numeric(subdat$VQSLOD) < -2.5 |as.numeric(subdat$MQ <40)))
  if(length(index)>0){
    subdat<-subdat[-index,]
  }
#  remove<-subdat
 #   subdat <- subdat[which(as.numeric(subdat$VQSLOD) > -10  ),] 
#  cnt<-( c(sum(subdat$AC_NFE)/total_control,sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control))
#  gfilter_rs<-rbind(gfilter_rs,c("VQSLOD",cnt))
  return (list(data=subdat,filter=gfilter_rs))
}