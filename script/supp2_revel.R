dat<-read.csv("~/Dropbox/CUMC/PAH-CHD/Revel_Compete/CHD_allvscontrol.csv",header=1,stringsAsFactors = F)
index<-c(3,2,5,6,1,8)
nms<-c("ALL","PolyPhen","CADD20","CADD25","MetaSVM","REVEL")
N_case=2645
N_control=1911
p=N_case/(N_case+N_control)
enrich_rate<-function(N_case,N_control,p){
  return(N_control*p/(N_case*(1-p)))
}


dat$Enich0<-0
dat$Enich_low<-0
dat$Enich_high<-0
for(i in 1:dim(dat)[1]){
  obv_case<-dat$Case[i]
  obv_control<-dat$Control[i]
  test<-binom.test(obv_case,(obv_case+obv_control),p)
  pe=as.numeric(test$estimate);
  enrich1<-enrich_rate(N_case,N_control,pe)
  conf<-as.numeric(test$conf.int)
  enrich_l<-enrich_rate(N_case,N_control,conf[1])
  enrich_h<-enrich_rate(N_case,N_control,conf[2])
  dat$Enich0[i]<-enrich1
  dat$Enich_low[i]<-enrich_l
  dat$Enich_high[i]<-enrich_h
}
tiff("~/Dropbox/PAH_CHD/Revel_competer.tiff", width=7, height=7, units="in", res=300, compression="lzw")
x<-1:length(nms)
plot(0,1,xlab = "",ylab="Enrichment",xaxt='n',xlim=c(1,6),ylim=c(1,2),type='n', bty="n")
axis(1,x,labels = nms,las=2)
j=1;
for(i in index){
  points(c(j,j,j),c(dat$Enich_low[i],dat$Enich0[i],dat$Enich_high[i]),type="l")
  points(j,dat$Enich0[i],type="p",pch=20,cex=3,col="red")
  j=j+1
}
dev.off()