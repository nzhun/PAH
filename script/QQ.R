tab<-read.csv("C:/Users/yixua/Dropbox/CUMC/PAH-CHD/PAH-CHD.gnomeAD.ALL.binom.csv", header=1,stringsAsFactors = F)
tab[which(tab$REVEL_LOF_Pvalue==-1),"REVEL_LOF_Pvalue"]=1
tab<-tab[order(tab$REVEL_LOF_Pvalue),]
tab<- tab[which(tab$REVEL_LOF_Pvalue>-1),] #

pvalue<- tab$REVEL_LOF_Pvalue #
#pvalue<-pvalue[which(pvalue<1)]
hist(pvalue,breaks = c(0,0.0001,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
xx<--log10((length(pvalue):1)/length(pvalue))
yy<- sort(-log10(pvalue)) 
pdf("C:/Users/yixua/Dropbox/CUMC/PAH-CHD/PAH-CHD-binom-QQ.pdf",width=7,height=7)
plot( xx ,yy,xlab="Expected p-value",ylab="Observed p-value",ylim=c(0,8),pch=20)
points( max(xx) ,max(yy),pch=5,col="red")
text(max(xx),max(yy),tab$Gene[1],pos = 2,col="red")

abline(0,1,lty=2)
abline(h=-log10(0.05/length(pvalue)),lty=2)
dev.off()