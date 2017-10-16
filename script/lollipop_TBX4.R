library(stringr)
lollipop <- function(fin,outf,domains,maxlen,type_len){
    fold=1
    folder=1 #.4 
    dev=5
    fdev=4
    dev2=2
    fdev2=2
    cfont=1.2
    vars <- read.table(fin,header = 1,stringsAsFactors =F,check.names = F,strip.white = T,sep = "\t")
    vars<-vars[order(vars$Group,decreasing = F),]
  
    
    pdf(outf,width=14,height=3)
    par(mar=c(0,0,0,0),mai=c(0,0,0,0),oma=c(0,0,0,0))
    plot(0,1,xlim=c(0,maxlen),ylim=c(-3.5,3.5),xlab="",ylab="", type='n',bty='n',axes=F)
 
    center=0
    width=0.6
    wbar=0.1
    dbar=0.25
    pcex=1.2
    pdist=0.27
    mut_types<-unique(tolower(vars$Mut_Type))
    col1<-as.vector(col2rgb("red")[,1]/255)
    col2<-as.vector(col2rgb(colors()[142])[,1]/255)
    col4<-as.vector(col2rgb(colors()[576])[,1]/255)
    fmut_types<-c("D-Mis","LGD","In-Frame","Missense")
    pointcols<-c("red","blue",colors()[c(142,576)])
    lcol<-"black"
  #  pointcols<-c(rgb(col1[1],col1[2],col1[3],alpha = 0.5),rgb(0,0,139/255,0.5),rgb(col4[1],col4[2],col4[3],0.5),rgb(col2[1],col2[2],col2[3],0.5))
    direction<- 1
    for(group in toupper(unique(vars$Group) )){
        pos<-vars[which(toupper(vars$Group)==group),]
        ## lines
        for(PRC in unique(pos$AA)){
          indexs<-which(pos$AA==PRC)
          inc=0;
          for( j in indexs){
            
            loc<-pos$AA[j]
            offset=0;
            types<-tolower(pos$Mut_Type[j])
            
            if(length(grep("missense",types,ignore.case = T))>0 && pos$MetaSVM[j]=="T"){
              next;
            }
                
              
              if(unique(loc)==0){print(paste(PRC,"lost"));next}
              if(length(indexs)>1){fold=folder}else{fold=1}
              y<-c(center,center+direction*(width*fold)*0.95) #*(1.2+(loc%%5)/8)+(loc%%3)/8+inc*0.06 *(length(indexs)/3+1)
              x<-as.numeric(loc) #+as.numeric(offset)/10
            #  if(length(which(pos$AA>x-2*dev & pos$AA<x))){
             #   offset=fdev*dev
            #  }#else if(length(which(pos$AA<x+2*dev & pos$AA>x))>0){offset=-fdev*dev }
              items<-sort(unique(pos$AA[which(pos$AA>x-dev & pos$AA<x+dev)]))
              if(length(items)>1){
                if(x==max(items)){offset=fdev}
                if(x==min(items)){offset=-fdev}
                else if(length(items)>3){offset=fdev}
                else{
                  items<-sort(unique(pos$AA[which(pos$AA>x-dev2 & pos$AA<x+dev2)]))
                  if(length(items)>1){ 
                    if(x==max(items)){offset=fdev2}
                    else if(x==min(items)){offset=-fdev2}
                    }
                }
              }
               points(c(x,x+offset),y,col=lcol,type='l',lwd=1)
               inc <- inc+1;
          }
        }
        
        ##points
        for(PRC in unique(pos$AA)){
          indexs<-which(pos$AA==PRC)
          inc=0;
          for( j in indexs){
            offset<-0
            loc<-pos$AA[j]
         #   offset<-pos$offset[j]
            types<-tolower(pos$Mut_Type[j])
            if(unique(loc)==0){print(paste(PRC,"lost"));next}
            
            if(length(grep("missense",types,ignore.case = T))>0){
              if(pos$MetaSVM[j]=="D"){
                pointcol=pointcols[1]
              }else{
                pointcol=pointcols[4]
                next;
              }
            }else if(length(grep("inframe|in_frame",types,ignore.case = T))>0){
              pointcol=pointcols[3]
            }else{
              pointcol=pointcols[2]
            }
            if(length(indexs)>1){fold=folder}else{fold=1}
            x<-as.numeric(loc) #+as.numeric(offset)/10
          #  if(length(which(pos$AA>x-2*dev & pos$AA<x))){
         #     offset=fdev*dev 
        #    }#else if(length(which(pos$AA<x+2*dev & pos$AA>x))>0){offset=-fdev*dev }
            items<-sort(unique(pos$AA[which(pos$AA>x-dev & pos$AA<x+dev)]))
            if(length(items)>1){
              if(x==max(items)){offset=fdev}
              if(x==min(items)){offset=-fdev}
              else if(length(items)>3){offset=fdev}
              else{
                items<-sort(unique(pos$AA[which(pos$AA>x-dev2 & pos$AA<x+dev2)]))
                if(length(items)>1){ 
                  if(x==max(items)){offset=fdev2}
                  if(x==min(items)){offset=-fdev2}
                }
              }
              
            }
            y<-c(center,center+(width*fold+inc*pdist )*direction) #+(loc%%3)/8+(loc%%5)/8+inc*0.07*(1.2)  (length(indexs)/3+1)
            points(x+offset,y[2],pch=19,type='p',cex=pcex,col=pointcol)
            points(x+offset,y[2],pch=1,type='p',cex=1.2*pcex,col="black")
            inc <- inc+1;
          }
        }
        
        
        
        text(7,direction*width*1.5,labels = str_to_title(group),cex=cfont)
        direction<-direction-2
    
    }
 
    rect(0,center-wbar,maxlen,center+wbar,col=border_col,border=border_col)  
    for(i in 1:dim(domains)[1]){
      domain=domains[i,]
      rect(as.numeric(domain$From),center-dbar,as.numeric(domain$To),center+dbar,col=domain_col[i],border=domain_col[i])
      text(as.numeric(domain$From)+(as.numeric(domain$To)/2-(as.numeric(domain$From))/2),center,
           labels = domain$Domain,col="black",cex=cfont)
      
      text(domain$From,center,labels = domain$From,srt=90,cex=cfont*0.85,co="black",font=2)
      text(domain$To,center,labels = domain$To,srt=90,cex=cfont*0.85,col="black",font=2)
    }
  #  axis(1,at=c(0,33,131,201,501,1038),labels = c(0,33,131,201,501,"1038 aa."),line = -3,tck=-0.02,col="gray",col.axis="gray")
    text(maxlen,y = center,labels=maxlen,srt=90,cex=cfont*0.85,col="black",font=2)
    text(maxlen+6,y = center,labels="aa.",srt=90,cex=cfont*0.85,col="black",font=2)
    legend("bottom",legend = fmut_types[1:type_len],
          horiz = T,pch=19,col=pointcols[1:type_len],cex=cfont,bty='n',pt.cex = pcex,y.intersp = 1)
    dev.off()
}


# maxlen=1038
# border_col<-"gray"
# domain_col<-c("green","#51C6E6","yellow","red","orange")
# domains<-as.data.frame(rbind(c("Activin",33,131),c("Pkinase",200,500)),stringsAsFactors = F)
# names(domains)<-c("Domain","From","To")
# fin<-"PAH/Result/Data/source/BMPR2_plot.txt"
# outf<-"~/Dropbox/CUMC/PAH/FInal/BMPR2_lolliplot.pdf"
# lollipop(fin,outf,domains,maxlen,type_len=2)



maxlen=545
border_col<-"gray"
domain_col<-c("#51C6E6","#FF8833","yellow","red","orange")
domains<-as.data.frame(rbind(c("T-BOX",66,252)),stringsAsFactors = F)
names(domains)<-c("Domain","From","To")
fin<-"PAH/Result/Data/source/TBX4_plot.txt"
outf<-"PAH/Result/Data/output/TBX4_lolliplot.pdf"
lollipop(fin,outf,domains,maxlen,type_len=3)


#library(trackViewer)
# 
# ### entropy
# library(entropy)
#   groups=unique(vars_plot$Group)
#   types<-unique(vars_plot$Type)
#   for(g in groups){
#     for(t in types){
#       scores<-c()
#       
#       for(aa in unique(vars_plot[which(vars_plot$Group==g & vars_plot$Type==t)])){
#         scores<-c(scores,length(which(vars_plot$Group==g & vars_plot$Type==t & vars_plot$AA==aa)))       
#       }
#     }
#   }
# ###