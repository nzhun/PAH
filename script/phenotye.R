setwd("~/server")

#setwd("/home/local/ARCS/nz2274/")
cohort<-read.csv("PAH/Result/Data/source/PAH-CHD-list.csv",header = 1,check.names = F,stringsAsFactors = F)
females<- grep("^F",cohort$Gender,ignore.case = T)
males<-grep("^M",cohort$Gender,ignore.case = T)
peds<-grep("^child|ped",cohort$TYPE,ignore.case = T)
adults<-grep("adult",cohort$TYPE,ignore.case = T)

ped_male<-length(intersect(males,peds))
ped_male_r<-length(intersect(males,peds))/length(peds)

ped_female=length(intersect(females,peds))
ped_female_r=length(intersect(females,peds))/length(peds)


adult_male<-length(intersect(adults,males))
adult_male_r<-adult_male/length(adults)

adult_female<-length(intersect(females,adults))
adult_female_r<-adult_female/length(adults)

n<-matrix(c(ped_male,adult_male,ped_female,adult_female),byrow = T,nrow = 2)

r<-matrix(format(c(ped_male_r,adult_male_r,ped_female_r,adult_female),digits = 2),byrow = T,nrow = 2)
c<-matrix(paste(n,r,sep=","),byrow = F,nrow = 2)

c<-(rbind(c,c(length(peds),length(adults))))

colnames(c)<-c("Pediatric","Adult")
rownames(c)<-c("Male","Female","total")
c
