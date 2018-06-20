cchmc<-read.csv("PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/PAH-Cinci-hg38-RGN-gnomad-inhouse-left.ALL.binom.csv",header=1,stringsAsFactors = F,comment.char = "")

#cchmc<-read.csv("PAH/PAH_10032017/Result/CCHMC_Internalcontrol_plus_gnomad/PAH-Cinci-hg38-RGN-gnomad-inhouse-left.ALL.binom.csv",header=1,stringsAsFactors = F,comment.char = "")
cumc<-read.csv("PAH/Result/CUMC_gnomad/PAH_CUMC_UNCCHMC_unsolved.ALL.binom.csv",header=1,stringsAsFactors = F,comment.char = "")
nms<-intersect(names(cumc),names(cchmc))
cchmc<-cchmc[,nms]
cumc<-cumc[,nms]
#names(cumc)<-names(cchmc)
v<-c()
c1<-1673
c2<-303 ##410
n=7509+5418
p<- (c1+c2)/(c1+c2+n)
all_rs<-c()
for(gene in unique(c(cchmc$Gene,cumc$Gene))){
  #print(gene)
  v1<-cchmc[which(cchmc$Gene==gene),]
  v2<-cumc[which(cumc$Gene==gene),]
  if(dim(v2)[1]==0){rs<-c(rs,v1);next}
  if(dim(v1)[1]==0){rs<-c(rs,v2);next}
  rs<-c(gene)
  for(i in seq(2,150,4)){
   # print(i)
    ncase<- as.numeric(v1[i]+v2[i])
    ncontrol<- as.numeric(v1[i+1])
    s<-Btest(ncase,c1+c2,ncontrol,n)
    rs<-c(rs,unlist(c(ncase,ncontrol,s$p.value,s$estimate)))
  }
  all_rs<-rbind(all_rs,rs)
}

all_rs<-as.data.frame(all_rs,stringsAsFactors = F)
names(all_rs)<-nms[1:153]
file<- "PAH/Result/CUMC_gnomad/merged.CUMC.CCHMC.twocontrol.binom.csv"
write.csv(all_rs,file = file,row.names = F)

intgenes<-c()
pdf(paste("PAH/Result/CUMC_gnomad/merged.CUMC.CCHMC.twocontrol.Binom.pdf",sep=""))
allp<- read.csv(file,header = T,comment.char = "",stringsAsFactors = F)
p<-allp[which(allp$Gene!="BMPR2" & allp$Gene!="HTT"& allp$Gene!="E2F1"),]
#arr<-unlist(lapply(1:dim(p)[1],FUN = function(x){sum(p[x,seq(3,54,by=4)])}))
#sickg<-p$Gene[which(arr==0)]
total<-dim(p)[2]-2
for(id in c(4,8,seq(28,total,by = 4))){
  key=names(p)[id]
  p<-p[order(p[,id]),]
  y=-log10(sort((p[,id]),decreasing = F))
  x=-log10((1:length(y))/length(y))
  z<-p[,id+1]
  plot(x,y,main="",xlab="",
       ylab="",pch=20,xlim=c(0,max(y,x)+0.1),
       ylim=c(0,max(y,x)+0.1),cex.axis=1.5,cex.lab=1.5)
  abline(0,1,lty=2,col="gray")
  abline(h=-log10(0.05/20000),lty=2,col="gray")
  title(main=key)
  title(xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",line = 2.5,cex.lab=1.5)
  index<-which(y>5 & z>1)
  if(length(index)>0){
  text(x[index],y[index],p$Gene[index],
       font=3,pos = 2,col="red",cex=1.5)
  points(x[index] ,y[index],pch=19,col="red",cex=2)
  }
  
}
dev.off()





###############  


cchmc<-read.csv("PAH/PAH_10032017/Result/CCHMC_gnomad_burden/PAH-Cinci-hg38-RGN-left.ALL.binom.csv",header=1,stringsAsFactors = F,comment.char = "")
cumc<-read.csv("PAH/Result/CUMC_gnomad/PAH_CUMC_UNCCHMC_unsolved.ALL.binom.csv",header=1,stringsAsFactors = F,comment.char = "")
#nms<-intersect(names(cumc),names(cchmc))
cchmc<-cchmc[,c(1:25,38:93)]
cumc<-cumc[,c(1:25,106:161)]
names(cumc)<-names(cchmc)
v<-c()
c1<-1673
c2<-410
n=7509
p<- (c1+c2)/(c1+c2+n)
all_rs<-c()
for(gene in unique(c(cchmc$Gene,cumc$Gene))){
  v1<-cchmc[which(cchmc$Gene==gene),]
  v2<-cumc[which(cumc$Gene==gene),]
  if(dim(v2)[1]==0){rs<-c(rs,v1);next}
  if(dim(v1)[1]==0){rs<-c(rs,v2);next}
  rs<-c(gene)
  for(i in seq(2,81,4)){
    
    ncase<- as.numeric(v1[i]+v2[i])
    ncontrol<- as.numeric(v1[i+1])
    s<-Btest(ncase,c1+c2,ncontrol,n)
    rs<-c(rs,unlist(c(ncase,ncontrol,s$p.value,s$estimate)))
  }
  all_rs<-rbind(all_rs,rs)
}

all_rs<-as.data.frame(all_rs,stringsAsFactors = F)
names(all_rs)<-names(cchmc)
file<- "PAH/Result/CUMC_gnomad/merged.CUMC.CCHMC.binom.csv"
write.csv(all_rs,file = file,row.names = F)

intgenes<-c()
pdf(paste("PAH/Result/CUMC_gnomad/merged.CUMC.CCHMC.Binom.pdf",sep=""))
allp<- read.csv(file,header = T,comment.char = "",stringsAsFactors = F)
p<-allp[which(allp$Gene!="BMPR2" & allp$Gene!="HTT"& allp$Gene!="E2F1"),]
#arr<-unlist(lapply(1:dim(p)[1],FUN = function(x){sum(p[x,seq(3,54,by=4)])}))
sickg<-p$Gene[which(arr==0)]
total<-dim(p)[2]-2
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
  if(y[1]>5){
    text(x[1],y[1],p$Gene[1],
         font=3,pos = 2,col="red",cex=1.5)
    points( max(x) ,max(y),pch=19,col="red",cex=2)
  }
  
}
dev.off()



