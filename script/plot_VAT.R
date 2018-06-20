#for file in c("REVEL/merged.filter.CCHMC.Internal.control.all.txt","MCAP/merged.filter.CCHMC.Internal.control.mcap.all.txt","CADD/merged.filter.CCHMC.Internal.control.mcap.all.txt"){
setwd("/home/local/ARCS/nz2274/")
library("metap")
source("Pipeline/NA_script/R/untils.R")
lib<-load_dataset()
exac<-load_exac()
sox17_target<-read.table("PAH/Result/Data/SOX17.target.txt",header=1,check.names = F,stringsAsFactors = F)
SOX17_genesets<- names(sox17_target)[2:dim(sox17_target)[2]]

fchd_gene="PAH/Result/Data/source/CHD.253GeneList.csv"

chd_genes<-read.csv( fchd_gene,header=1,stringsAsFactors=F,check.names=F,comment.char="",quote="")
chd_genes<-chd_genes$Gene

hchannel<-read.table("PAH/documents/Channel_highlyExps.txt",header=1,stringsAsFactors = F)
hchannel_geneset<-hchannel$GeneName

all_ion_ch<-read.table("PAH/documents/All_ion_plus_channel.txt",header=F,stringsAsFactors = F,fill=T,sep="\t")
all_ion_ch_geneset<-all_ion_ch$V1


ion_ch<-read.csv("PAH/documents/ION_Chanel.anno.csv",header=1,stringsAsFactors = F)
ion_ch_geneset<-ion_ch$GeneName

pota<-read.table("PAH/documents/Potassium.channle.txt",header=F,stringsAsFactors = F,fill=T)
pota_geneset<-pota$V1



channel<-read.csv("PAH/documents/Channelopathy_gene_list.csv",header=1,stringsAsFactors = F)
channel_geneset<-channel$GeneName

smc<-read.csv("PAH/documents/human_SMC_genelists.csv",header=1,stringsAsFactors = F)
smc_geneset<-smc$GeneName


smad_target<-read.table("Resources/GeneSets/SMAD_genesets.txt",header=F,stringsAsFactors = F)
smad_geneset<-smad_target$V1

tgf_signal<-read.table("PAH/documents/TGF-Beta.pathway.genes.list",header=F,stringsAsFactors = F)
tgf_geneset<-tgf_signal$V1

wnt_signal<-read.table("PAH/documents/Wnt-signal.pathway.genes.list",header=F,stringsAsFactors = F)
wnt_geneset<-wnt_signal$V1

bmp_signal<-read.table("PAH/documents/BMP_signaling-HSA-201451].gene.list",header=F,stringsAsFactors = F)
bmp_geneset<-bmp_signal$V1


smad4_down<-read.table("PAH/documents/Transcriptional_activity_of_SMAD2_SMAD3:SMAD4_heterotrimer_Homo_sapiens_R-HSA-2173793.gene.list",header=F,stringsAsFactors = F)
smad4_down_geneset<-smad4_down$V1



VEGF_pathway<-read.table("PAH/documents/VEGF_gene_set.txt",header=F,stringsAsFactors = F)
VEGF_geneset<-VEGF_pathway$V1



listDLOF <- read.table("PAH/documents/PAH-UK-DMIS-LOF.gene.list",stringsAsFactors = F,header = F)

list_dmis <- read.table("PAH/documents/PAH-UK-DMIS.gene.list",stringsAsFactors = F,header = F)
list_skat <- read.table("PAH/documents/PAH-UK-SKATO.gene.list",stringsAsFactors = F,header = F)
list_lof <- read.table("PAH/documents/UK-PTV-TopGenes.txt",stringsAsFactors = F,header = F)

UKs<-c(listDLOF$V1,list_dmis$V1,list_lof$V1,list_skat$V1)

#folder="PAH/PAH_10032017/hg38/VAT/"

myplot <- function(keyf,outpf,folder){

#outpfs=c("CCHMC_ALL.internal.plus.gnomad.autosome","CCHMC_APAH.internal.plus.gnomad.autosome","CCHMC_IPAH.internal.plus.gnomad.autosome")  #c("CCHMC_ALL.internal.plus.gnomad","CCHMC_APAH.internal.plus.gnomad","CCHMC_IPAH.internal.plus.gnomad") #c("gnomad","InterNalControl","APAH-internal","APAH-gnomad","IPAH-internal","IPAH-gnomad")
#keyfs<-c("CCHMC_ALL.internal.plus.gnomad.autosome.all.txt","CCHMC_APAH.internal.plus.gnomad.autosome.all.txt","CCHMC_IPAH.internal.plus.gnomad.autosome.all.txt")  #c("CCHMC_ALL.internal.plus.gnomad.all.txt","CCHMC_APAH.internal.plus.gnomad.all.txt","CCHMC_IPAH.internal.plus.gnomad.all.txt")   #c("CCHMC.internal.plus.gnomad.all.txt")  #c("merged.filter.CCHMC.gnomad.all.txt","merged.filter.CCHMC.Internal.control.all.txt","APAH.InternalControl.all.txt","APAH.gnomad.all.txt","IPAH.InternalControl.all.txt","IPAH.gnomad.all.txt")
#for(i in 1:length(keyfs)){
#eyf
      #  outpf=paste(keyf,"_out",sep="")
        #keyf="merged.filter.CCHMC.Internal.control.all.txt"
        f<-"REVEL"
        EN=200;        
        pdf(paste(folder,"/VAT.",outpf,".pdf",sep=""))
    #    for(f in types){
          raw_dat<-read.table(paste(keyf,sep=""),stringsAsFactors = F)
          names(raw_dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
        	dat<-raw_dat
        	names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
        	dat<-dat[order(dat[,"permutp"]),]
        #	dat<-dat[which(dat$V4<1),]
        	lambda<- qchisq(median(c(raw_dat[,"permutp"],rep(1,20000-dim(raw_dat)[1]))),1,lower.tail = F)/0.456
        	y<--log10(c(dat[,"permutp"],rep(1,20000-dim(dat)[1])))
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main=paste("ALL"),xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
     #    #print(f)
          abline(0,1,col="gray")
          abline(h=-log10(0.05/20000),col="gray")
          index<-which(dat[,"permutp"]< 0.05/20000)
          if(length(index)>0){
    #        #print(dat[index,])
            points(x[index],y[index],pch=20,col="red",cex=1.5)
            text(x[index],y[index],dat[index,1],pos=2,cex=1)
           # legend("bottomright",legend=dat[index,1],cex=0.75,bty='n')
          }
          
          exp_dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""),stringsAsFactors = F)
          names(exp_dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          exp_dat<-exp_dat[order(exp_dat[,"permutp"])[1:EN],]
          exp_dat$FDR<-as.numeric(p.adjust(exp_dat$permutp,method = "BH",n = 18000))
          exp_dat <- gene_exp_endothelial(exp_dat,lib$endothelial_exp)
          exp_dat <- gene_exp_Heart(exp_dat,lib$heart_exp )
          exp_dat <- gene_exp_Lung(exp_dat,lib$lung_exp )
          exp_dat <- gene_exp_SMC(exp_dat,lib$smc_exp )
          exp_dat <- gene_zscore(exp_dat,exac)
          exp_dat <- gene_exp_diaphragm(exp_dat,lib$diaph_exp )
      
          
          cchmc<- raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""),stringsAsFactors = F)
          names(cchmc)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          for(tkey in c("UK_raremissense.Binomal","UK_PAH_PTV_Binomical","UK_raremissense_PTV.Binomal")){
            uk<-read.table(paste("PAH/documents/",tkey,".txt",sep=""),stringsAsFactors = F,header=1)
            
            arr<-c()
            for(g in uk$Gene){
              index<-which(uk$Gene==g)
              cid<-which(cchmc$GeneName==g)
              if(length(cid)>0){
                #   #print(cchmc$permutp[cid])
                ps<- c(uk$p.value[index],cchmc$permutp[cid])
                if(length(ps) <2){next;}
                s<-sumlog(ps)
                arr<-rbind(arr,c(g,s$chisq,s$df,s$p,unlist(s$validp)))
              }
            }
            if(is.null(arr)|| dim(arr)[1]<1){next;}
            if(dim(arr)[2] <4){next}
            write.csv(arr[order(as.numeric(arr[,4])),],file = paste(folder,"/",outpf,".",tkey,".csv",sep=""))
            
            plot(1:dim(arr)[1],-log10(as.numeric(arr[,4])),xaxt='n',pch=20,xlab="",main=paste("Meta",tkey,sep=""),ylab="-log10(meta pvalue)",col="red",ylim=c(1,11),xlim=c(0,dim(arr)[1]))
            for(i in 1:dim(arr)[1]){
              points(c(i,i,i),c(-log10(as.numeric(arr[i,4:6]))),xaxt='n',type='l',lty=2,xlab="",ylab="-log10(meta pvalue)",col="gray")
            }
            
            points(1:dim(arr)[1],-log10(as.numeric(arr[,5])),xaxt='n',pch=3,xlab="",ylab="-log10(meta pvalue)",col="green")
            
            points(1:dim(arr)[1],-log10(as.numeric(arr[,6])),xaxt='n',pch=4,xlab="",ylab="-log10(meta pvalue)",col="blue")
            index<-which(as.numeric(arr[,4])<0.001)
            text(index,-log10(as.numeric(arr[index,4])),labels = arr[index,1],cex=1,pos=c(1,3))
            abline(h=-log10(0.05/20000))
            legend("topright",legend = c("Meta","UK","CCHMC"),col=c("red","green","blue"),pch=c(20,3,4),bty='n')
            
          } 
          
          
     #   }
   #     dev.off()
        
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%hchannel_geneset),]
          exp_dat$Hchannel_target<-0
          exp_dat$Hchannel_target[which(exp_dat[,1]%in%hchannel_geneset)]<-1
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="CHANNEL_HIGH_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
          #      #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
            #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          #  }
          #dev.off()
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%channel_geneset),]
          exp_dat$channel_target<-0
          exp_dat$channel_target[which(exp_dat[,1]%in%channel_geneset)]<-1
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="CHANNEL_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
          #      #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
            #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          #  }
          #dev.off()
          
          
          
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%smc_geneset),]
          exp_dat$SMC_target<-0
          exp_dat$SMC_target[which(exp_dat[,1]%in%smc_geneset)]<-1
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="SMC_celline_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
          #      #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
            #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          #  }
          #dev.off()
          
   
          
 
          
          pota<-read.table("PAH/documents/Potassium.channle.txt",header=F,stringsAsFactors = F,fill=T)
          pota_geneset<-pota$V1
          
          exp_dat$ALl_ION_CHANNEL<-0
          exp_dat$ALl_ION_CHANNEL[which(exp_dat[,1]%in%all_ion_ch_geneset)]<-1
          
          exp_dat$ION_CHANNEL<-0
          exp_dat$ION_CHANNEL[which(exp_dat[,1]%in%ion_ch_geneset)]<-1
          
          
          exp_dat$POTASSIUM<-0
          exp_dat$POTASSIUM[which(exp_dat[,1]%in%pota_geneset)]<-1
          
       # pdf(paste(folder,"/VAT-CHD.",outpf,".pdf",sep=""))
      #  for(f in types){
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%chd_genes),]
          exp_dat$CHD_target<-0
          exp_dat$CHD_target[which(exp_dat[,1]%in%chd_genes)]<-1
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="CHD_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
    #      #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
  #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
      #  }
        #dev.off()
        
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%VEGF_geneset),]
          exp_dat$VEGF_target<-0
          exp_dat$VEGF_target[which(exp_dat[,1]%in%VEGF_geneset)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="VEGF_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
  #        #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
#            #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%smad4_down_geneset),]
          exp_dat$SMAD4_down_target<-0
          exp_dat$SMAD4_down_target[which(exp_dat[,1]%in%smad4_down_geneset)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="SMAD4_down_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
 #         #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
    #        #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          
          
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%bmp_geneset),]
          exp_dat$BMP_target<-0
          exp_dat$BMP_target[which(exp_dat[,1]%in%bmp_geneset)]<-1
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="BMP_signal_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
   #       #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
  #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          
          
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%wnt_geneset),]
          exp_dat$WNT_target<-0
          exp_dat$WNT_target[which(exp_dat[,1]%in%wnt_geneset)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="Wnt_signal_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
   #       #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
  #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
        
          
          
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%tgf_geneset),]
          exp_dat$TGFB_target<-0
          exp_dat$TGFB_target[which(exp_dat[,1]%in%tgf_geneset)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="TGFB_signal_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
    #      #print(f)
          abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
    #        #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          
          
          
       
          
     #   pdf(paste(folder,"/VAT-SOX17.",outpf,".pdf",sep=""))
       # for(f in types){
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%SOX17_genesets),]
          exp_dat$SOX17_target<-0
          exp_dat$SOX17_target[which(exp_dat[,1]%in%SOX17_genesets)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="SOX17_targets",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
   #       #print(f)
          abline(0,1,col="gray")
          index<-which(dat[,"permutp"]<0.001)
          if(length(index)>0){
  #          #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
    #    }
     #   dev.off()
        
        #smad_geneset
        
    #    pdf(paste(folder,"/VAT-SMAD.",outpf,".pdf",sep=""))
     #   for(f in types){
          dat<-raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%smad_geneset),]
          
          exp_dat$SMAD_target<-0
          exp_dat$SMAD_target[which(exp_dat[,1]%in%smad_geneset)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="SMAD_target",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
 #         #print(f)
          abline(0,1,col="gray")
          
          abline(h=-log10(0.05/20000),col="gray")
          
          index<-which(dat[,"permutp"]<0.001)
          if(length(index)>0){
   #         #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
    #   }
   #     dev.off()
        
      #  pdf(paste(folder,"/VAT-UKs.",outpf,".pdf",sep=""))
       # for(f in types){
          dat<- raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""))
          names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
          dat<-dat[order(dat[,"permutp"]),]
          dat<-dat[which(dat[,1]%in%UKs),]
          
          
          exp_dat$UK<-0
          exp_dat$UK[which(exp_dat[,1]%in%UKs)]<-1
          
          y<--log10(dat[,"permutp"])
          x<--log10(c(1:length(y))/length(y))
          plot(x,y,main="UK_genes",xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
 #         #print(f)
          abline(0,1,col="gray")
          abline(h=-log10(0.05/20000),col="gray")
          index<-which(dat[,"permutp"]<0.001)
          if(length(index)>0){
   #         #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            legend("right",legend=dat[index,1],cex=0.75)
    #      }
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
   #     dev.off()
        
     #   pdf(paste(folder,"/VAT_Endoth_",outpf,".pdf",sep=""))
          exp_dat$hit<-unlist(lapply(1:dim(exp_dat)[1],FUN = function(x){return (sum(exp_dat[x,grep("target|UK",names(exp_dat))]))}))
          sum_arr<-unlist(lapply(grep("target|UK",names(exp_dat)),FUN = function(x){return(sum(exp_dat[,x]))}))
          exp_dat[EN+1,]<-rep(-1,dim(exp_dat)[2])
          exp_dat[EN+1,grep("target|UK",names(exp_dat))]<- sum_arr
        
     #   for(f in types){
         
          dat<-exp_dat
          # dat<-dat[,c(1:5)]
         # names(dat)<-  c(,"Endothelial")
          dat<-dat[order(dat$Endothelial),]
          y<--log10(dat[,"permutp"])
          x<-dat$Endothelial #-log10(c(1:length(y))/length(y))
          plot(x,y,main="Endothelial",xlab="Endothelial expression",ylab="-log10 permut-p",pch=20)
   #       #print(f)
         # abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
    #        #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          #  legend("right",legend=dat[index,1],cex=0.75)
          }else{
            index=1
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
          }
          
          dat<-dat[order(dat[,"permutp"]),]
          qut=66
          for(qut in c(0,33,66,90)){
            subdat<-dat[which(dat[,"Endothelial"]>qut),]
            y<--log10(subdat[,"permutp"])
            x<- -log10(c(1:length(y))/length(y))
        
            plot(x,y,main=paste(f,"\nendothelial expression top",(100-qut),length(y)),xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
            abline(0,1,col="gray")
            abline(h=-log10(0.05/20000),col="gray")
            index<-which(y > 3)
            if(length(index)>0){
       #       #print(subdat[index,])
              text(x[index],y[index],subdat[index,1],pos=1:4,cex=1)
              #  legend("right",legend=dat[index,1],cex=0.75)
            }
          }
      #  }
     #   dev.off()
        
        
        
        
        #pdf(paste(folder,"/VAT_SMC_",outpf,".pdf",sep=""))
        
        
       # for(f in types){
      #    dat<- raw_dat # read.table(paste(folder,"/",f,"/",keyf,sep=""),stringsAsFactors = F)
      #    names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
      #    dat<-dat[order(dat[,"permutp"])[1:EN],]
      #    dat <- gene_exp_SMC(dat,lib$smc_exp)
         # dat<-dat[,c(1:5)]
         # names(dat)<-  c("GeneName","BestZ","BestT","permutp","SMC")
          dat<-dat[order(dat$SMC),]
          y<--log10(dat[,"permutp"])
          x<-dat$SMC #-log10(c(1:length(y))/length(y))
          plot(x,y,main=f,xlab="Smooth muscle expression",ylab="-log10 permut-p",pch=20)
    #      #print(f)
          # abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
     #       #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            #  legend("right",legend=dat[index,1],cex=0.75)
          }
          
          dat<-dat[order(dat[,"permutp"]),]
          qut=66
          for(qut in c(0,33,66,90)){
            subdat<-dat[which(dat[,"SMC"]>qut),]
            y<--log10(subdat[,"permutp"])
            x<- -log10(c(1:length(y))/length(y))
            plot(x,y,main=paste(f,"\nSMC expression top",(100-qut),length(y)),xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
            abline(0,1,col="gray")
            abline(h=-log10(0.05/20000),col="gray")
            index<-which(y > 3)
            if(length(index)>0){
        #      #print(subdat[index,])
              text(x[index],y[index],subdat[index,1],pos=1:4,cex=1)
              #  legend("right",legend=dat[index,1],cex=0.75)
            }
          }
       # }
     #   dev.off()
        
        
        
      #  pdf(paste(folder,"/VAT_lung_",outpf,".pdf",sep=""))
        
        
      #  for(f in types){
     #     dat<- raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""),stringsAsFactors = F)
    #      names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
    #      dat<-dat[order(dat[,"permutp"])[1:EN],]
          
    #      dat <- gene_exp_Lung(dat,lib$lung_exp)
         # dat<-dat[,c(1:5)]
          dat<-dat[order(dat$LUNG_EXP),]
          y<--log10(dat[,"permutp"])
          x<-dat$LUNG_EXP #-log10(c(1:length(y))/length(y))
          plot(x,y,main=f,xlab="Lung expression", ylab="-log10(Observed VT-pvalue)",pch=20)
      #    #print(f)
          # abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
            #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            #  legend("right",legend=dat[index,1],cex=0.75)
          }
          
          dat<-dat[order(dat[,"permutp"]),]
          qut=66
          for(qut in c(0,33,66,90)){
            subdat<-dat[which(dat[,"LUNG_EXP"]>qut),]
            y<--log10(subdat[,"permutp"])
            x<- -log10(c(1:length(y))/length(y))
            plot(x,y,main=paste(f,"\n lung expression top",(100-qut),length(y)),xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
            abline(0,1,col="gray")
            abline(h=-log10(0.05/20000),col="gray")
            index<-which(y > 3)
            if(length(index)>0){
              #print(subdat[index,])
              text(x[index],y[index],subdat[index,1],pos=1:4,cex=1)
              #  legend("right",legend=dat[index,1],cex=0.75)
            }
          }
           
     #   }
      #  dev.off()
        
        
        
      #  pdf(paste(folder,"/VAT_heart",outpf,".pdf",sep=""))
        
        
      #  for(f in types){
    #      dat<- raw_dat #read.table(paste(folder,"/",f,"/",keyf,sep=""),stringsAsFactors = F)
  #        names(dat)<-c("GeneName","BestZ","BestT","permutp","Nvar")
   #       dat<-dat[order(dat[,"permutp"])[1:EN],]
          
    #      dat <- gene_exp_Heart(dat,lib$heart_exp)
         # dat<-dat[,c(1:5)]
          dat<-dat[order(dat$HEART_EXP),]
          y<--log10(dat[,"permutp"])
          x<-dat$HEART_EXP #-log10(c(1:length(y))/length(y))
          plot(x,y,main="HeartExpression",xlab="Heart expression", ylab="-log10(Observed VT-pvalue)",pch=20)
          #print(f)
          # abline(0,1,col="gray")
          index<-which(dat$permutp<0.001)
          if(length(index)>0){
            #print(dat[index,])
            text(x[index],y[index],dat[index,1],pos=1:4,cex=1)
            #  legend("right",legend=dat[index,1],cex=0.75)
          }
          
          
          dat<-dat[order(dat[,"permutp"]),]
          qut=66
          for(qut in c(0,33,66,90)){
            subdat<-dat[which(dat[,"HEART_EXP"]>qut),]
            y<--log10(subdat[,"permutp"])
            x<- -log10(c(1:length(y))/length(y))
            plot(x,y,main=paste(f,"\n heart expression top",(100-qut),length(y)),xlab="-log10 (Expected pvalues)", ylab="-log10(Observed VT-pvalue)",pch=20)
            abline(0,1,col="gray")
            abline(h=-log10(0.05/20000),col="gray")
            index<-which(y > 3)
            if(length(index)>0){
              #print(subdat[index,])
              text(x[index],y[index],subdat[index,1],pos=1:4,cex=1)
              #  legend("right",legend=dat[index,1],cex=0.75)
            }
          }
     #   }
   #     dev.off()
      #  }
          
     
          col=2
          for(col in c(2,4)){
            exp_dat[,col]<-as.numeric(format(exp_dat[,col],digits = 2,scientific = T))
          }
          for(col in 6:14){
            exp_dat[,col]<-as.numeric(format(as.numeric(exp_dat[,col]),digits = 2,scientific = F))
          }
          exp_dat<- exp_dat[,c(1:4,6:12,17:19,21:22,27:31,5,13:16,20,23:26)]
          write.csv(exp_dat,file = paste(folder,"/",outpf,".top.anno.csv",sep=""),row.names = F)   
          
}

folder="PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/VAT/REVEL/Binomial/"
keyf="PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/VAT/REVEL/Binomial/843.APAH.CCHMC.REVEL.all.txt"
out_prefix="843.APAH.CCHMC.REVEL"
args = commandArgs(trailingOnly=TRUE)
keyf=args[1]
out_prefix=args[2]
folder =args[3]
print(args) 
myplot(keyf,out_prefix,folder)
