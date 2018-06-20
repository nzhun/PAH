setwd("~/server/")
rm_dat<-read.csv("PAH/PAH_10032017/src/Excluding_CCHMC.txt",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
rm_all<-read.csv("PAH/PAH_10032017/src/Excluding_CCHMC_all.csv",header=1,stringsAsFactors = F,check.names = F,comment.char = "")
pheno_dat<-read.table("PAH/PAH_10032017/src/PAH_10032017_1230.Affected.ped",header = 1,stringsAsFactors = F,check.names = F,comment.char = "",sep="\t",quote="")
pheno_dat<-pheno_dat[which(!pheno_dat$ID%in% rm_dat$`#PatientID`),]
pop<-read.csv("PAH/PAH_10032017/src/PAH.pop.csv",header = 1,stringsAsFactors = F,comment.char = "",check.names = F)
pheno_pop<-merge(pheno_dat,pop,by.x="ID",by.y = "proband")
pheno_pop$Subtrype[which(pheno_pop$Subtrype=="")]<-pheno_pop$PHENOTYPE[which(pheno_pop$Subtrype=="")]
pheno_pop$Subtrype[which(pheno_pop$Subtrype=="PAH")]<-pheno_pop$PHENOTYPE[which(pheno_pop$Subtrype=="PAH")]
pheno_dat<-pheno_pop
tab1<-c();
sexs<-unique(pheno_pop$SEX)
enths<-unique(pheno_pop$popback)
regions<-unique(pheno_pop$ETHNICITY)
for(type in unique(pheno_pop$Subtrype)){
   rowc<-c()
   sub_pheno<-pheno_pop[which(pheno_pop$Subtrype==type),]
   ## gender ##
   for(g in sexs){
     rowc<-c(rowc,length(which(sub_pheno$SEX==g)))
   }
   for(eh in enths){
     rowc<-c(rowc,length(which(sub_pheno$popback==eh)))
   }
   for(rg in regions){
     rowc<-c(rowc,length(which(sub_pheno$ETHNICITY==rg)))
   }
   rowc<-c(rowc,length(which(sub_pheno$age_Dx <19 & sub_pheno$age_Dx >-1)))
   
   rowc<-c(rowc,length(which(sub_pheno$age_Dx >18)))
    tab1<-cbind(tab1,rowc)  
}
tab1<-as.matrix(tab1)
colnames(tab1)<-unique(pheno_pop$Subtrype)
row.names(tab1)<-c(sexs,enths,regions,"child (dx age <19)","adult (dx age >=19)")
write.csv(tab1,"~/Dropbox (CGC)/US PAH Genetic Manuscript 2018/table_figures/table1_phenotype.csv",row.names = T)


plot_tab<-function(synd,title){
    plot(0,1,type='n',xlim=c(1,length(synd)),ylim=c(0,max(synd)),xaxt='n',xlab="",ylab="",bty='n',main=title)
    for(i in 1:length(synd)){
      lines(c(i,i),c(0,synd[i]),lwd=20)
    }
  by=1
  if(length(synd)>10){
      by=as.integer(length(synd)/10)
  }
    axis(1,at = seq(1,length(synd),by =by),labels = names(synd)[seq(1,length(synd),by=by)],cex=0.7,las=2)
  
}
synd<-table(pheno_dat$PHENOTYPE)
plot_tab(synd,"Phenotyp")

synd<-table(pheno_dat$SEX)
barplot(synd,col="white")
plot_tab(synd,"SEX")
text(x= 1.5,y=500, labels =paste("F:M=", format(synd[2]/synd[3],digits = 2)))


synd<-table(pheno_dat$SEX[which(pheno_dat$PHENOTYPE=="APAH")])
plot_tab(synd,"APAH_SEX")
text(x= 1.5,y=500, labels =paste("F:M=", format(synd[1]/synd[2],digits = 2)))

synd<-table(pheno_dat$SEX[which(pheno_dat$PHENOTYPE=="IPAH")])
plot_tab(synd,"IPAH_SEX")
text(x= 1.5,y=500, labels =paste("F:M=", format(synd[2]/synd[3],digits = 2)))

#table(pheno_dat$IsForDiagnosticRHC)
pheno_dat$Dx_date<-pheno_dat$DateRHCDone
pheno_dat$Dx_date[which(pheno_dat$IsForDiagnosticRHC=="N")]<-pheno_dat$DateNonDxRHC[which(pheno_dat$IsForDiagnosticRHC=="N")]
pheno_dat$Dx_age<-as.integer(format(as.Date.character(pheno_dat$Dx_date,format = "%m/%d/%Y"),"%Y"))-pheno_dat$YOB
ages<-pheno_dat$Dx_age;
ages[is.na(ages)]<-"-1"
barplot(table(ages))
plot_tab(table(ages),title = "Dx_age")
boxplot(ages,main="Dx_age")
ages<-2017-pheno_dat$YOB
plot_tab(table(ages),title = "age")


ages<-pheno_dat$Dx_age[which(pheno_dat$PHENOTYPE=="APAH")];
ages[is.na(ages)]<-"-1"
plot_tab(table(ages),title = "Dx_age_APAH")
ages<-2017-pheno_dat$YOB
plot_tab(table(ages),title = "age_APAH")


ages<-pheno_dat$Dx_age[which(pheno_dat$PHENOTYPE=="IPAH")];
ages[is.na(ages)]<-"-1"
plot_tab(table(ages),title = "Dx_age_IPAH")
ages<-2017-pheno_dat$YOB
plot_tab(table(ages),title = "age_IPAH")


table(pheno_dat$pop)

