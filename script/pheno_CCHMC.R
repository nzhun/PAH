dat<-read.csv("~/server/PAH/PAH_10032017/PAH_10032017.denovo_filtered.csv",header = 1,stringsAsFactors = F,comment.char = "")
dat$ProbandName<-unlist(lapply(dat$proband,FUN = function(x) unlist(strsplit(x,split = "(",fixed = T) )[1]))
write.table(dat[order(dat$ProbandName),c("ProbandName","CHROM","POS")],file = "~/server/PAH/PAH_10032017/Rare_variants/denovo.query.txt",row.names = F,sep="\t",quote = F)



rt<-read.table("~/server/PAH/PAH_10032017/plink/PAH_10032017.IBD.genome",header=1,stringsAsFactors = F)
sub_rt<-rt[which(rt$PI_HAT>0.05),]
pdf("~/Dropbox/PAH_Cincinatti/Relatedness.pdf")
plot(sub_rt$Z0,sub_rt$Z1,xlab="Z0",ylab="Z1")
dev.off()

sex<-read.table("~/Dropbox/PAH_Cincinatti/PAH_10032017.sexcheck.sexcheck.txt",header=1,check.names = F,stringsAsFactors = F)

#############  sex check ##########
phenotype<-read.table("PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",header=1,strip.white = T,stringsAsFactors = F,sep = "\t",comment.char = "",check.names = F)

#phenotype<-read.table("PAH/PAH_10032017/PAH_10032017.ped",header=1,strip.white = T,stringsAsFactors = F,sep = "\t",comment.char = "",check.names = F)
#bcr_phenotype<-read.csv("BreastCancer/BCR_2017/Variant_data/BCR_2017_phenotype.csv",header=1,check.names = F,stringsAsFactors = F)

######## phenotye  annotatation ###########
sex$Gender=0
sex$Affected=0
sex$AGE=0
sex$syndrome=0
for(i in 1:dim(phenotype)[1]){
  key1=phenotype$ID[i] #substr(phenotype$ID[i],nchar(phenotype$ID[i])-15,nchar(phenotype$ID[i]))
  key2=phenotype$Father[i] #substr(phenotype$Father[i],nchar(phenotype$Father[i])-15,nchar(phenotype$Father[i]))
  
  x1<-which(sex$IID==key1)
  if(length(x1)>0){
    sex$IID[x1]=phenotype$ID[i]
    sex$Gender[x1]=phenotype$SEX[i]
    sex$Affected[x1]=phenotype$Affected[i]
    sex$AGE[x1]=phenotype$`Current Age`[i]
    sex$syndrome[x1]=phenotype$PHENOTYPE[i]
  }
  
}
write.csv(sex[grep("CCH",sex$IID),],file = "~/Dropbox/PAH_Cincinatti/PAH_2017.gender.check.csv",quote = F,row.names = F)





#ped<-read.table("PAH/PAH_10032017/src/PAH_10032017.Affected.ped",sep="\t",header=1,stringsAsFactors = F,comment.char = "",check.names = F)
ped<-read.table("PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",sep="\t",header=1,stringsAsFactors = F,comment.char = "",check.names = F)
rm<-read.csv("PAH/PAH_10032017/src/Excluding_CCHMC.txt",header = 1,check.names = F,stringsAsFactors = F)
ped$SID<-unlist(lapply(ped$ID,FUN = function(x){(substr(x,nchar(x)-15,nchar(x)))}))
ibd<-read.table("PAH/PAH_10032017/plink/PAH_10032017.IBD.genome",header=1,stringsAsFactors = F)

ibd<-ibd[which(ibd$PI_HAT>0.15),]
ibd<-ibd[which(ibd$FID1%in%ped$SID | ibd$FID2%in%ped$SID ),]
ibd$Affected=1
ibd$PAH="-"
ibd$FamState="-"
ibd$Age=0;
ibd$RELAT="";
ibd$Affected[which(ibd$IID1%in%ped[which(ped$Affected==2),"SID"])]=2
ibd$PAH[which(ibd$IID1%in%ped[which(ped$Affected==2),"SID"])]="ID1"
ibd$Affected[which(ibd$IID2%in%ped[which(ped$Affected==2),"SID"])]=2
#ibd$Gender[which(ibd$IID2%in%ped[which(ped$Affected==2),"SID"])]=0
ibd$PAH[which(ibd$IID2%in%ped[which(ped$Affected==2),"SID"] & ibd$PAH!="-")]="ID1,ID2"
ibd$PAH[which(ibd$IID2%in%ped[which(ped$Affected==2),"SID"] & ibd$PAH=="-")]="ID2"

for(i in 1:dim(ped)[1]){
  key1=substr(ped$ID[i],nchar(ped$ID[i])-15,nchar(ped$ID[i]))
  key2=substr(ped$Father[i],nchar(ped$Father[i])-15,nchar(ped$Father[i]))
  key3=substr(ped$Mother[i],nchar(ped$Mother[i])-15,nchar(ped$Mother[i]))

  x1<-which(ibd$IID1==key2)
  if(length(x1)>0){ibd$IID1[x1]=ped$Father[i]}
  x2<-which(ibd$IID2==key2)
  if(length(x2)>0){ibd$IID2[x2]=ped$Father[i]}
  
  
  
  x1<-which(ibd$IID1==key3)
  if(length(x1)>0){ibd$IID1[x1]=ped$Mother[i]}
  x2<-which(ibd$IID2==key3)
  if(length(x2)>0){ibd$IID2[x2]=ped$Mother[i]}
  
  x1<-which(ibd$IID1==key1)
  if(length(x1)>0){ibd$IID1[x1]=ped$ID[i]}
  x2<-which(ibd$IID2==key1)
  if(length(x2)>0){ibd$IID2[x2]=ped$ID[i]}
  ibd$FamState[c(x2,x1)]=paste(ped$FAMID[i],ibd$FamState[c(x2,x1)],sep=" ")
  ibd$RELAT[c(x2,x1)]=paste(ped$RELATIONSHIP_TO_PROBAND[i],sep="")
  ibd$Age[c(x2,x1)]=paste(ped$`Current Age`[i],ibd$Age[c(x2,x1)],sep=" ")
  #ibd$Gender[which(ibd$IID2%in%ped[which(ped$Affected==2),"ID"])]=2
}

ibd<-ibd[which(!ibd$IID1%in%rm$`#PatientID` & !ibd$IID2%in%rm$`#PatientID` ),]
rts<-c()
arr<-c()
for(i in 1:dim(ped)[1]){
  famf=c(ped$ID[i],ped$Father[i])
  famm=c(ped$ID[i],ped$Mother[i])
  index2<-which(ibd$IID1%in% famf  & ibd$IID2 %in% famf )
  if(length(index2)>0) ibd$RT[index2]="father-child"
  else {
    if(famf[2]!="X"){rts<-c(rts,"Father -child "); arr<-rbind(arr,famf)}
  }
  index2<-which(ibd$IID1%in% famm  & ibd$IID2 %in% famm )
  if(length(index2)>0) ibd$RT[index2]="mother-child"
  else{
    if(famm[2]!="X" ){
      rts<-c(rts,"Mother -child ");
      arr<-rbind(arr,famm)
    }
  }
}

fail_TR<-cbind(arr,rts)
pdf("PAH/PAH_10032017/QC/Relatedness_PAH.pdf")
key="CCHMC"
index<-intersect(which(ibd$PI_HAT>0.1 ),intersect(grep("CCHMC",ibd$IID1),grep("CCHMC",ibd$IID2)))
dg1<-intersect(which(ibd$PI_HAT>0.4 & ibd$PI_HAT<0.6  ) ,intersect(grep("CCHMC",ibd$IID1),grep("CCHMC",ibd$IID2)))## parents
#dg11<-intersect(which(ibd$PI_HAT>0.4 & ibd$PI_HAT<0.6 &ibd$RT=="D1" ) ,intersect(grep("PAH",ibd$IID1),grep("PAH",ibd$IID2)))## parents
dg2<-intersect(which(ibd$PI_HAT>0.2 & ibd$PI_HAT < 0.4 ),intersect(grep("CCHMC",ibd$IID1),grep("CCHMC",ibd$IID2)))
dg0<-intersect(which(ibd$PI_HAT>0.6 & ibd$Z2>0.8),intersect(grep("CCHMC",ibd$IID1),grep("CCHMC",ibd$IID2)))

trueR<-which(ibd$RT!="UN")
plot(ibd$Z0[index],ibd$Z1[index],xlab="proportion of IBD=0", ylab="proportion of IBD=1")
points(ibd$Z0[dg2],ibd$Z1[dg2],col="blue")
points(ibd$Z0[dg1],ibd$Z1[dg1],col="red")
points(ibd$Z0[trueR],ibd$Z1[trueR],col="cyan",pch=3)
points(ibd$Z0[dg0],ibd$Z1[dg0],col="green")
legend("topright",legend = c("replicates","parent-child","sibling","others","In_doc"),fill=c("green","red","blue","black","cyan"),bty = 'n')

excep_trio<-c("JM238","JM725","JM990","JM145_")

dev.off()


write.csv(ibd[dg0,],file = "PAH/PAH_10032017/QC/Replicates.csv",quote = T,row.names = F)
write.csv(ibd[dg1,],file = "PAH/PAH_10032017/QC/Parent_child.csv",quote = T,row.names = F)
write.csv(ibd[dg2,],file = "PAH/PAH_10032017/QC/sibling.csv",quote = T,row.names = F)
write.csv(ibd[c(intersect(which(ibd$PI_HAT>0.4 &ibd$RT=="UN" ),c(intersect(grep("CCHMC",ibd$IID1),grep(key,ibd$IID2)))),
                intersect(which(ibd$PI_HAT<0.2 &ibd$RT!="UN" ),c(intersect(grep(key,ibd$IID1),grep(key,ibd$IID2))))
),],file = "PAH/PAH_10032017/QC/question_realtions.csv",quote = T,row.names = F)
write.csv(fail_TR,file="PAH/PAH_10032017/QC/question_realtions_fail.csv",quote = T,row.names = F)
#write.table(ibd[ ,],file = "PAH/PAH_2017/Relatedness/question_realtions.csv",append = T,quote = F,sep = "\t",row.names = F)


## inconsistent samples ####
incon<-read.table("PAH/PAH_10032017/src/PAH_10032017_inconsistent.txt",header = 1,stringsAsFactors = F)
cchmc_ped<-read.table("PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",header=1,stringsAsFactors = F,quote="",sep="\t",comment.char = "",check.names = F)
cchmc_ped<-cchmc_ped[,c(1:110)]
cchmc_ped <- merge(x=cchmc_ped,y=incon,by.x="ID",by.y="ID",all.x = T)
cchmc_ped$age_Dx[!is.na(cchmc_ped$corrected_dx_age)]<-cchmc_ped$corrected_dx_age[!is.na(cchmc_ped$corrected_dx_age)]
cchmc_ped$Subtrype[!is.na(cchmc_ped$corrected_phenotype)] <- cchmc_ped$corrected_phenotype[!is.na(cchmc_ped$corrected_phenotype)]
write.table(cchmc_ped[,c(1:110)],"PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",row.names = F,sep = "\t",quote = F)


incon<-read.table("PAH/Result/Data/source/PAH_CUMC_CCHMC_inconsistent.txt",header = 1,stringsAsFactors = F)
cchmc_ped<-read.csv("PAH/Result/Data/source/PAH-CHD-list.csv",header=1,stringsAsFactors = F,comment.char = "",check.names = F,fill = T)
cchmc_ped<-cchmc_ped[,c(1:dim(cchmc_ped)[2])]
cchmc_ped <- merge(x=cchmc_ped,y=incon,by.x="ID",by.y="ID",all.x = T)
cchmc_ped$Age_dx[!is.na(cchmc_ped$corrected_dx_age)]<-cchmc_ped$corrected_dx_age[!is.na(cchmc_ped$corrected_dx_age)]
cchmc_ped$Disease[!is.na(cchmc_ped$corrected_phenotype)] <- cchmc_ped$corrected_phenotype[!is.na(cchmc_ped$corrected_phenotype)]
write.csv(cchmc_ped[,c(1:13)],"PAH/Result/Data/source/PAH-CHD-list.csv",row.names = F,quote = F)


incon<-read.table("PAH/PAH_10032017/RiskGenes/additional.genes.txt",header = 1,stringsAsFactors = F,sep="\t",comment.char = "")
cchmc_ped<-read.table("PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",header=1,stringsAsFactors = F,quote="",sep="\t",comment.char = "",check.names = F)
pop<-read.csv("PAH/PAH_10032017/src/PAH.pop.csv",header = 1,stringsAsFactors = F,comment.char = "",check.names = F)
cchmc_ped<-merge(cchmc_ped,pop,by.x="ID",by.y = "proband")
#cchmc_ped<-cchmc_ped[,c(1:110)]
incon <- merge(y=cchmc_ped,x=incon,by.y="ID",by.x="proband",all.x = T)
cchmc_ped$age_Dx[!is.na(cchmc_ped$corrected_dx_age)]<-cchmc_ped$corrected_dx_age[!is.na(cchmc_ped$corrected_dx_age)]
cchmc_ped$Subtrype[!is.na(cchmc_ped$corrected_phenotype)] <- cchmc_ped$corrected_phenotype[!is.na(cchmc_ped$corrected_phenotype)]
write.table(incon,"PAH/PAH_10032017/RiskGenes/additional.genes.anno.txt",row.names = F,sep = "\t",quote = F)
