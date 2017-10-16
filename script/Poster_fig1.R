dat<-data.frame(pop=c("European","Hispanic","African","South Asian","East Asian","Other","Female","Male"),Pediatric=c(81,30,13,10,7,3,91,53),
                Adult=c(63,37,6,7,7,3,88,24))
dat<-dat[dim(dat)[1]:1,]
tiff("~/Dropbox/CUMC/PAH/FInal/fig1_poster.tiff",width=960,height=490,res=100)
par(mar=c(6,10,1,1))
barplot(t(dat[,c(2,3)]),beside = T,horiz = T,names.arg = dat$pop,las=2,cex.names =2,cex.axis = 2,col=c("black","white"))
#abline(h=6.5,lty=2)
title(xlab="# Patients",cex.lab=2,line=4.5)
legend("right",legend = c("Pediatric(144)","Adult(112)"),fill=c("black","white"),bty='n',cex=1.5)
dev.off()