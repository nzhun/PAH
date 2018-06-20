path="/home/local/ARCS/nz2274/" #"~/server/"
setwd(paste(path,"PAH/Result/",sep=""))

known<-read.csv("Data/contribution.csv",header=1,check.names = F,stringsAsFactors = F)
total<-c("Total",25,135,79,178)
#known<-rbind(known,total)
known<-known[unlist(lapply(1:dim(known)[1],FUN = function(x)if(sum(as.numeric(known[x,2:5])) >0)x)),]
#index<-c(1:,12)
## solved samples ##
cols<-colors()[c(56,637,494,618,472,373,468,656,382,344,630,373,76,24, 595 ,443,41)]
pdf(paste("Data/output/Contribution_four_published.pdf",sep=""),width = 11,height=7,family ="Times")
par(mai=c(1,0,0,0),mfrow=c(2,2),oma=c(3,2,2,5),mar=c(2,2,3,2 ),xpd=NA)
#par(mfrow=c(2,2),mar=c(1, 0, 2, 0), oma = c(1, 0, 1, 0))
arr<-known[,2:5]
#arr<-as.numeric(as.matrix(arr))
indexs<-c("A.","B.","C.","D.")
for(i in 1:length(colnames(arr))){
  nm<-colnames(arr)[i]
  nm_s<- unlist(strsplit(unlist(strsplit(nm,split = "(",fixed = T))[2],split=")",fixed=T))[1]
  vec<-as.numeric(arr[,i])
  dn=2
  if(sum(vec)>100){dn=1}
  ptg<-format(vec/sum(vec)*100,digits = dn,width = 2)
  #print(ptg)
  ids<-which(ptg!=" 0" &ptg!=" 0.00"&ptg!=" 0.0")
  ptg<-paste(ptg,"%",sep="") #rownames(arr)[ids]," ",
  pie(vec[ids],col=cols[ids],labels = ptg[ids],radius = 1.05,clockwise =1,init.angle = -156, angle = 90,cex=1.2)
  if(i==1){title(main="Familial",line=2.5,cex.main=2);title(ylab="Child",line=0,cex.lab=2,font.lab=2)}
  if(i==2){title(main="Non-familial",line=2.5,cex.main=2);}
  if(i==3){title(ylab="Adult",line=0,cex.lab=2,font.lab=2)}
  
  pos="right"
 # title(main=paste(indexs[i],nm),line=1.2,cex.main=1.5)
  title(xlab=nm_s,line = 2,cex.lab=1.5)
  yx=0.2
  if(i==3){yx=0.5}
  #legend(-2,yx,legend =ptg[ids],fill=cols[ids],bty="n",cex = 1.05,pt.cex = 0.85)
  
}
nms<-known$Gene
#nms[which(nms=="Rare variants")]="Rare variants\ncandidates"

lpos=1.55
#nms[which(nms=="De Novo")]=substitute(paste(italic("De Novo"),"\ncandidates"))
#plot(0,1,xlim=c(0,5),ylim=c(0,3),type='n',bty='n',axes = F)
legend(lpos,3.5,legend =nms,fill=cols,bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
text(lpos+0.8,-0.1,"candidates",cex=1.3)
# i=9
# op <- par(font = 3)
# lpos=1.55
# legend(lpos,3.5,legend =nms[1:i],fill=cols[1:i],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
# op <- par(font = 1)
# #legend(-7.6,1,legend =expression(paste(italic("De Novo"),  " \ncandidates")),fill=cols[9],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.2,pt.cex = 1,xpd = NA)
# i=i+1
# #legend(lpos,0,legend =substitute(paste(italic("De novo"))),fill=cols[i],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
# #text(lpos+0.7,-0.45,"candidates",cex=1.3)
# legend(lpos,0.6,legend =substitute(paste(italic("De novo"))),fill=cols[i],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
# text(lpos+0.6,0.2,"candidates",cex=1.3)
# 
# i=i+1
# #legend(lpos,-0.45,legend =nms[i],fill=cols[i],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
# legend(lpos,0.25,legend =nms[i],fill=cols[i],bty="n",x.intersp = 1,y.intersp = 1.9,cex = 1.3,pt.cex = 1,xpd = NA)
#write.table(x = arr,file = "~/server/PAH_CU_VU/publish4_contribution.csv",sep = "\t")
dev.off()
