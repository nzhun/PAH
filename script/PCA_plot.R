#uff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.gz.plus.HapMap"
setwd("/home/local/ARCS/nz2274/")
#setwd("~/server/")
#folder<-"~/server/BreastCancer/AJ/VCF/PCA_HQ/"
#f_ped<-"~/server/BreastCancer/AJ/Result/AJ_case_control.ped"
#suff<-"AJ_20151202.hardfiltered.716cases.704controls.merged.commonBed.coding.HQ.gz.plus.HapMap"
pf="Resources/1000Gneome_2013/source/integrated_call_samples_v3.20130502.plusAJ.plusDom.ALL.panel"
#pf="Resources/1000Gneome_2013/source/integrated_call_samples_v3.20130502.ALL.panel"
cont<-read.table(pf,header=1,fill = T,strip.white = T,sep="\t",check.names=F,stringsAsFactors = F)

sids<-c("WIC373",
        "WIC434",
        "WIC2376",
        "WIC2414",
        "WIC67",
        "WIC152",
        "WIC177",
        "WIC328",
        "WIC410",
        "WIC438",
        "WIC627",
        "WIC2295",
        "WIC2318",
        "WIC2344",
        "WIC2374",
        "WIC2414",
        "WIC2400",
        "WIC2568",
        "WIC3003",
        "WIC3009")
#cont<-read.table("Resources/Control/InHouse_Control/plink/temp.AJ.txt",header=1,fill = T,strip.white = T,sep="\t",check.names=F,stringsAsFactors = F)
#names(cont)<-c("sample","pop","super_pop","gender")


myplot <- function(folder,suff) {
   
    f_eval<-paste(folder,suff,".eval",sep="")
    f_evec<-paste(folder,suff,".pca.evec",sep="")
   #  ped<-read.table(f_ped)
    
    
    eval <- read.table(f_eval)
    evec1.pc <- round(eval[1,1]/sum(eval)*100,digits=2)
    evec2.pc <- round(eval[2,1]/sum(eval)*100,digits=2)
    evec3.pc <- round(eval[3,1]/sum(eval)*100,digits=2)
    evec4.pc <- round(eval[4,1]/sum(eval)*100,digits=2)
    
    evec <- read.table(f_evec)
    
   ###### hapmap ########
    inds<-as.vector(mapply(as.character(evec$V1),FUN = function(x) unlist(strsplit(x,":"))[1]))
    evec$ID<-c()
    pops<-c()
    pred_pop<-c();
    for(i in 1:length(inds)){
      id=inds[i]
      if(nchar(id)<15){  
        index<-which(cont$sample==id);
      }else{
        index<-grep(paste(id,"$",sep=""),cont$sample)
      }
      if(length(index)<1){
        pops<-c(pops,"InHouse")
        evec$ID[i]<-id
      }else{
        pops<-c(pops,cont$pop[index])
        if(length(index)>1){print (i);stop()}
        evec$ID[i]<-cont$sample[index]
      }
    }
    pred_pop<-pops  
    ##### hapmap #####
    
    
    cols<-c("gray","red","dark blue","orange","cyan","violet",colors()[c(371,417,382,376,373,417,435,430,425,466,498,551,555,574,371)])
    pchs<-c(3,4,3,1,2,c(1:6),c(1:6))
    
    ## pops
    #hapmap<-grep("HapMap",evec$V1,ignore.case = T) #which(grep("HapMap",evec$V1,ignore.case = T))
   # local<- grep("HapMap",evec$V1,ignore.case = T,invert = T) # which(grep("HapMap",evec$V1,ignore.case = T,invert = T))
   # a<-unlist(lapply(as.character(evec$V1[hapmap]),function(x)(unlist(strsplit(as.character(x),split = ":")))[1]))
  #  pops<-unlist(lapply(a,function(x) unlist(strsplit(x,".",fixed = T))[2]))
    
 
    snms<-unlist(lapply(1:dim(evec)[1],FUN = function(x){unlist(strsplit(as.character(evec$V1[x]),split = ":"))[1]}))
    
    pdf(paste(folder,"/",suff,".PCA.pdf",sep=""))
    par(mar = c(5.1, 5.1, 2.1, 2.1))
    for (e in 3:5){
      plot(evec[,e], evec[,2], 
           xlab=paste("eigenvector",e-1,"\n",round(eval[e-1,1]/sum(eval)*100,digits=2), "% of observed genetic variation", sep=""), 
           ylab=paste("eigenvector1\n",round(eval[1,1]/sum(eval)*100,digits=2), "% of observed genetic variation", sep=""), 
           main="PCA",
           col=cols[1],
           pch=pchs[1],cex=0.5 # ,type='n'
           # ,xlim=xr,ylim=yr
           # ,  xlim=c(0.005,0.015),
           # ylim=c(0.002,0.012)
      ) 
      pp<-c()
      lnm<-c("Inhouse-Case")
      i=2
      d=0;
      for(c in unique(cont$super_pop)){
        j=1;
        for(p in unique(cont$pop[which(cont$super_pop==c)])){
          index<-which(pops==p)
          pp<-c(pp,index)
          #  d=d+length(index)
          points(evec[index,e], evec[index,2], 
                 col=cols[i],
                 pch=pchs[j],
                 cex=0.7
          )
          j=j+1;
        }
        all_index<-which(pops%in% unique(cont$pop[which(cont$super_pop==c)]))
        const=1
        if(c=="AFR"){const=0.999}
        
        i=i+1
        lnm<-c(lnm,c)
      }
      
      legend("bottomleft",legend=lnm,fill=cols[1:(length(lnm))],bty='n')
      legend("bottomright",legend="WIC",col="purple",pch=3,bty='n')  
      ids<-grep("TTR|JM0002|JM1341|JM1363|WIC",evec$V1)
      if(length(ids)>0){
        points(evec[ids,e], evec[ids,2], 
               col="purple",
               pch=3,
               cex=3
        )
        qids<-which(snms%in%sids)
        text(evec[qids,e], evec[qids,2], 
             labels = evec[qids,1],
             cex=0.5
        ) 
        
      }
    }
   for (e in 2:5){
    plot(evec[,e], evec[,e+1], 
         xlab=paste("eigenvector",e-1,"\n", round(eval[e-1,1]/sum(eval)*100,digits=2), "% of observed genetic variation", sep=""), 
         ylab=paste("eigenvector",e,"\n", round(eval[e,1]/sum(eval)*100,digits=2), "% of observed genetic variation", sep=""), 
         main="PCA",
         col=cols[1],
         pch=pchs[1],cex=0.5 # ,type='n'
         # ,xlim=xr,ylim=yr
         # ,  xlim=c(0.005,0.015),
         # ylim=c(0.002,0.012)
    ) 
    pp<-c()
    lnm<-c("Inhouse-Case")
    i=2
    d=0;
    for(c in unique(cont$super_pop)){
      j=1;
      for(p in unique(cont$pop[which(cont$super_pop==c)])){
        index<-which(pops==p)
        pp<-c(pp,index)
        #  d=d+length(index)
        points(evec[index,e], evec[index,e+1], 
               col=cols[i],
               pch=pchs[j],
               cex=0.7
        )
      
        j=j+1;
      }
      all_index<-which(pops%in% unique(cont$pop[which(cont$super_pop==c)]))
      const=1
      if(c=="AFR"){const=0.999}
    
      
      i=i+1
      lnm<-c(lnm,c)
    }
    
    legend("topleft",legend=lnm,fill=cols[1:(length(lnm))],bty='n')
    legend("top",legend="WIC",col="purple",pch=3,bty='n')  
     ids<-grep("TTR|JM0002|JM1341|JM1363|WIC",evec$V1)
    if(length(ids)>0){
      points(evec[ids,e], evec[ids,e+1], 
            col="purple",
            pch=3,
            cex=3
     ) 
      qids<- which(snms%in%sids)
     text(evec[qids,e], evec[qids,e+1], 
              labels = evec[qids,1],
          cex=0.5
     ) 
    }
    # points(evec[ids,2], evec[ids,3], 
    #        col="black",
    #        pch=3,
    #        cex=0.5
    # ) 
     
   }
    

    ########
    # plot(evec[,3], evec[,4], 
    #      xlab=paste("eigenvector2\n",evec2.pc, "% of observed genetic variation", sep=""), 
    #      ylab=paste("eigenvector3\n",evec3.pc, "% of observed genetic variation", sep=""), 
    #      main="PCA",
    #      col=cols[1],
    #      pch=pchs[1],cex=1 #,type='n'
    #      # ,xlim=xr,ylim=yr
    #      # ,  xlim=c(0.005,0.015),
    #      # ylim=c(0.002,0.012)
    # ) 
    # pp<-c()
    # lnm<-c("Inhouse-Case")
    # i=2
    # d=0;
    # for(c in unique(cont$super_pop)){
    #   j=1;
    #   for(p in unique(cont$pop[which(cont$super_pop==c)])){
    #     index<-which(pops==p)
    #     pp<-c(pp,index)
    #     points(evec[index,3], evec[index,4], 
    #            col=cols[i],
    #            pch=pchs[j],
    #            cex=1.2
    #     )
    #     if(p=="FIN"){
    #       points(evec[index,3], evec[index,4], 
    #              col="purple",
    #              pch=20,
    #              cex=1.2
    #       )
    #     }
    #     j=j+1;
    #   }
    #   all_index<-which(pops%in% unique(cont$pop[which(cont$super_pop==c)]))
    #   const=1
    #   bounds<-c(sort(evec[all_index,4],decreasing = F)[as.integer(const*length(all_index))],  ### maxy
    #             sort(evec[all_index,4],decreasing = T)[as.integer(const*length(all_index))], ### miny
    #             sort(evec[all_index,3],decreasing = T)[as.integer(const*length(all_index))], ### minx
    #             sort(evec[all_index,3],decreasing = F)[as.integer(const*length(all_index))] ### maxx
    #   )
    #   #abline(h=bounds[1],col=cols[i],lty=2) #
    #   #abline(h=bounds[2],col=cols[i],lty=2)
    #   #abline(v=bounds[4],col=cols[i],lty=2)
    #   #abline(v=bounds[3],col=cols[i],lty=2)
    #   if(c=="SAS"|| c=="AMR"){ 
    #     ids<-which(evec[,4]>bounds[2]& evec[,4]<bounds[1] & evec[,3]>bounds[3] & evec[,3]<bounds[4]  & pred_pop=="InHouse");
    #     pred_pop[ids]<- c
    #     
    #     points(evec[ids,3], evec[ids,4], 
    #            col=cols[i],
    #            pch=20,
    #            cex=2
    #     ) 
    #     points(evec[ids,3], evec[ids,4], 
    #            col="yellow",
    #            pch=20,
    #            cex=0.7
    #     ) 
    #     
    #   }
    #   i=i+1
    #   lnm<-c(lnm,c)
    # }
    # 
    # legend("bottomright",legend=lnm,fill=cols[1:(length(lnm))],bty='n')
    # legend("top",legend="WIC",col="purple",pch=3,bty='n')  
    # if(length(ids)>0){
    #   points(evec[ids,3], evec[ids,4], 
    #          col="purple",
    #          pch=3,
    #          cex=3
    #   ) 
    #   qids<-which(snms%in%sids)
    #   text(evec[qids,3], evec[qids,4], 
    #        labels = evec[qids,1],
    #        cex=0.5
    #   ) 
    # }
    # 
    # 
    # ########
    # plot(evec[,5], evec[,4], 
    #      xlab=paste("eigenvector4\n",evec4.pc, "% of observed genetic variation", sep=""), 
    #      ylab=paste("eigenvector3\n",evec3.pc, "% of observed genetic variation", sep=""), 
    #      main="PCA",
    #      col=cols[1],
    #      pch=pchs[1],cex=1 #,type='n'
    # 
    # ) 
    # pp<-c()
    # lnm<-c("Inhouse-Case")
    # i=2
    # d=0;
    # for(c in unique(cont$super_pop)){
    #   j=1;
    #   for(p in unique(cont$pop[which(cont$super_pop==c)])){
    #     index<-which(pops==p)
    #     pp<-c(pp,index)
    #     #  d=d+length(index)
    #     points(evec[index,5], evec[index,4], 
    #            col=cols[i],
    #            pch=pchs[j],
    #            cex=1.2
    #     )
    #     if(p=="FIN"){
    #       points(evec[index,5], evec[index,4], 
    #              col="purple",
    #              pch=20,
    #              cex=1.2
    #       )
    #     }
    # 
    #     j=j+1;
    #   }
    #   all_index<-which(pops%in% unique(cont$pop[which(cont$super_pop==c)]))
    #   bounds<-c(sort(evec[all_index,4],decreasing = F)[as.integer(0.99*length(all_index))],  ### maxy
    #             sort(evec[all_index,4],decreasing = T)[as.integer(0.99*length(all_index))], ### miny
    #             sort(evec[all_index,5],decreasing = T)[as.integer(0.99*length(all_index))], ### minx
    #             sort(evec[all_index,5],decreasing = F)[as.integer(0.99*length(all_index))] ### maxx
    #   )
    #   #abline(h=bounds[1],col=cols[i],lty=2) #
    #   #abline(h=bounds[2],col=cols[i],lty=2)
    #   #abline(v=bounds[4],col=cols[i],lty=2)
    #   #abline(v=bounds[3],col=cols[i],lty=2)
    #   #if(c=="SAS"){ ids<-which(eve[,4]>bounds[2]& eve[,4]<bounds[1] & eve[,5]>bounds[3] & eve[,5]<bounds[4]  & pred_pop!="InHouse"); pred_pop[ids]<- c }
    #   i=i+1
    #   lnm<-c(lnm,c)
    # }
    # 
    # legend("bottomright",legend=lnm,fill=cols[1:(length(lnm))],bty='n')
    # legend("top",legend="WIC",col="purple",pch=3,bty='n')  
    # ids<-grep("TTR|JM0002|JM1341|WIC",evec$V1)
    # if(length(ids)>0){
    #   points(evec[ids,5], evec[ids,4], 
    #          col="purple",
    #          pch=3,
    #          cex=3
    #   ) 
    #   qids<-which(snms%in%sids)
    #   text(evec[qids,5], evec[qids,4], 
    #        labels = evec[qids,1],
    #        cex=0.5
    #   ) 
    # }
    # 
    # ########
    # plot(evec[,6], evec[,5], 
    #      xlab=paste("eigenvector6\n", "% of observed genetic variation", sep=""), 
    #      ylab=paste("eigenvector5\n", "% of observed genetic variation", sep=""), 
    #      main="PCA",
    #      col=cols[1],
    #      pch=pchs[1],cex=1 #,type='n'
    #      
    # ) 
    # pp<-c()
    # lnm<-c("Inhouse-Case")
    # i=2
    # d=0;
    # for(c in unique(cont$super_pop)){
    #   j=1;
    #   for(p in unique(cont$pop[which(cont$super_pop==c)])){
    #     index<-which(pops==p)
    #     pp<-c(pp,index)
    #     #  d=d+length(index)
    #     points(evec[index,6], evec[index,5], 
    #            col=cols[i],
    #            pch=pchs[j],
    #            cex=1.2
    #     )
    #     if(p=="FIN"){
    #       points(evec[index,6], evec[index,5], 
    #              col="purple",
    #              pch=20,
    #              cex=1.2
    #       )
    #     }
    #     
    #     j=j+1;
    #   }
    #   all_index<-which(pops%in% unique(cont$pop[which(cont$super_pop==c)]))
    #   bounds<-c(sort(evec[all_index,5],decreasing = F)[ceiling(0.99*length(all_index))],  ### maxy
    #             sort(evec[all_index,5],decreasing = T)[floor(0.99*length(all_index))], ### miny
    #             sort(evec[all_index,6],decreasing = T)[floor(0.99*length(all_index))], ### minx
    #             sort(evec[all_index,6],decreasing = F)[ceiling(0.99*length(all_index))] ### maxx
    #   )
    #   #abline(h=bounds[1],col=cols[i],lty=2) #
    #   #abline(h=bounds[2],col=cols[i],lty=2)
    #   #abline(v=bounds[4],col=cols[i],lty=2)
    #   #abline(v=bounds[3],col=cols[i],lty=2)
    #   #if(c=="SAS"){ ids<-which(eve[,4]>bounds[2]& eve[,4]<bounds[1] & eve[,5]>bounds[3] & eve[,5]<bounds[4]  & pred_pop!="InHouse"); pred_pop[ids]<- c }
    #   i=i+1
    #   lnm<-c(lnm,c)
    # }
    # 
    # legend("right",legend=lnm,fill=cols[1:(length(lnm))],bty='n')
    # legend("top",legend="WIC",col="purple",pch=3,bty='n')  
    # ids<-grep("TTR|JM0002|JM1341|JM1363|WIC",evec$V1)
    # if(length(ids)>0){
    #   points(evec[ids,6], evec[ids,5], 
    #          col="purple",
    #          pch=3,
    #          cex=3
    #   ) 
    #   qids<-which(snms%in%sids)
    #   text(evec[qids,6], evec[qids,5], 
    #        labels = evec[qids,1],
    #        cex=0.5
    #   ) 
    # }
    
    dev.off()
    
   dpop<-data.frame(ID=inds,pop=pops,pred_pop=pred_pop,stringsAsFactors=F)
   write.csv(dpop,paste(folder,"/",suff,"PCA.csv",sep=""),row.names = F)
   
}

args<-commandArgs(TRUE)
folder=args[1];
suff=args[2];
folder=paste(folder,"/",sep="");
myplot(folder,suff)



#t<-read.table("server/others/MEGA/MEGA.plus.refere.balance.rm.summary.txt")
t<-read.table("server/others/MEGA/MEGA.plus.refere.balance.rm.higher0.05.summary.txt")
t<-t[1:3049,]
hist(t$V3,xlab="#SNPs per sample")

t1<-t[1:401,]
hist(t1$V3,xlab="#SNPs per sample",freq = F,main="#SNP",xlim = c(1000,3000),col="red",ylim=c(0,0.006))

t2<-t[402:2905,]
hist(t2$V3,xlab="#SNPs per sample",freq = F,main="1KG",add=T,col="green")

t3<-t[2906:3049,]
hist(t3$V3,xlab="#SNPs per sample",main="MEGA",freq = F,add=T,col="blue")

legend("topleft",legend = c("1KG","AJ+DOM","MEGA"),fill=c("green","red","blue"),bty='n')
legend("top",legend="WIC",col="purple",pch=3,bty='n')  
tp<-read.table("server/others/MEGA/MEGA.plus.refere.balance.rm.snp.0.05.summary.txt")

hist(-log10(tp$V6))

hist(-log10(tp$V7))
hist(-log10(tp$V8))

