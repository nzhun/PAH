setwd("~/server/")
source("Pipeline/NA_script/R/untils.R")

source("Pipeline/Enrichment/Denovo_enrichment_ZN.R")
exac<-load_exac()
sets<-load_dataset() ## embryo expression
lib <- load_genesets(sets,exac)


test_enrich<-function(dat,samples,geneset){
  ## de novo in known gene
  denovo_enrich<-c()
  denovo_adult_known<-scienceTable(length(samples),dat,as.character(geneset),F,c())
  
  denovo_adult_known$type="knowngene"
  #  denovo_adult_known$group=group
  
  ## de novo in general 
  
  denovo_adult_general<-scienceTable(length(samples),dat,"all",F,c())
  denovo_adult_general$type="allgene"
  # denovo_adult_general$group=group
  denovo_enrich<-rbind(denovo_enrich,denovo_adult_general)
  
  
  denovo_adult_other<-scienceTable(length(samples),dat,"all",F,geneset)
  denovo_adult_other$type="othergene"
  #  denovo_adult_other$group=group
  denovo_enrich<-rbind(denovo_enrich,denovo_adult_other)
  ### do the gene sets test
  
  denovo_genesets_enrich<- denovo_enrichment2(length(samples),dat,geneset,lib)
  # denovo_genesets_enrich$group<-group
  denovo_enrich<-rbind(denovo_enrich,denovo_genesets_enrich)
  return (denovo_enrich)
}


#mutrate<-load_mutrate()

process<-function(fdenovo,weslist,fasso,fsolved,fped,fout,fout2){

    solved_vars<-read.csv(fsolved,header = 1,stringsAsFactors = F,check.names = F,strip.white = T)
 
    
    denovo<-read.table(fdenovo,header =  1,sep="\t",stringsAsFactors = F,check.names = F,comment.char = "")
    denovo<-formatFreq_new(denovo)
    denovo<-filter_allfreq_new(denovo,0.0001,0.001)
    denovo<-map_filter(denovo)
    
    ped<-read.table(fped,header=1,stringsAsFactors = F,sep="\t",fill = T,check.names = F,comment.char = "")
    
    #fasso<- "PAH/PAH-CHD/source/PAH-CHD.candidate.gene.txt"
    gene_asso  <-  read.table(fasso,header = F,stringsAsFactors = F,check.names = F)
    gene_asso  <-  gene_asso$V1
    
    solved.samples<-solved_vars[,c("proband","Gene.refGene")]
    names(solved.samples)<-c("ID","GeneName")
    solved_samples<-solved.samples$ID[which(solved.samples$GeneName%in%gene_asso)]
    
   
    ## denovo enrichment
    ##trio_ped
    ped_trio<-ped[which(nchar(ped$Father)>2 & nchar(ped$Mother)>2 & ped$ID%in%weslist$proband & ped$Affected==2),]
    denovo<-denovo[which(denovo$ProbandName%in%ped_trio$ID),]
    
    knowngenes<-gene_asso #knowngenes$V1
    total_set<-dim(ped_trio)[1]
    children_list<-weslist$proband[grep("child",weslist$TYPE,ignore.case = T)]
    adult_list<-weslist$proband[grep("adult",weslist$TYPE,ignore.case = T)]
    
    pediatric<-ped_trio$ID[which(ped_trio$ID%in%children_list)]
    adult<-ped_trio$ID[which(ped_trio$ID%in%adult_list)]
    
    denovo$Type="-"
    denovo$Type[which(denovo$ProbandName%in%adult)]="adult"
    denovo$Type[which(denovo$ProbandName%in%pediatric)]="child"
    
    #denovo$group<-"PAH-others"
    #denovo$group[which(denovo$ProbandName%in%chd_cases)]<-"PAH-CHD"
    
    denovo_adult<-denovo[which(denovo$ProbandName%in%adult),]
    denovo_pediatric<-denovo[which(denovo$ProbandName%in%pediatric),]
    
    #denovo_chd<-denovo[which(denovo$ProbandName%in%chd_cases),]
    #denovo_other<-denovo[which(denovo$ProbandName%in%other_cases),]
    
    
    result<-test_enrich(denovo,unique(ped_trio$ID),knowngenes)
    result$group<-"PAH-ALL"
    
    
    resultA<-test_enrich(denovo_adult,adult,knowngenes)
    resultA$group<-"adult"
    
    resultc<-test_enrich(denovo_pediatric,pediatric,knowngenes)
    resultc$group<-"pediatric"
    
    denovo_enrich<-rbind(result,resultA,resultc)
    
    ######## unsovled 
    uresult<-test_enrich(denovo[which(!denovo$ProbandName%in%solved_samples),],unique(ped_trio$ID[which(!ped_trio$ID%in%solved_samples)]),knowngenes)
    uresult$group<-"PAH-ALL-Unsolved"
    
    uresultA<-test_enrich(denovo_adult[which(!denovo_adult$ProbandName%in%solved_samples),],adult[which(!adult%in%solved_samples)],knowngenes)
    uresultA$group<-"adult_unsolved"
    
    
    uresultc<-test_enrich(denovo_pediatric[which(!denovo_pediatric$ProbandName%in%solved_samples),],pediatric[which(!pediatric%in%solved_samples)],knowngenes)
    uresultc$group<-"pediatric_unsolved"
    
    denovo_enrich<-rbind(denovo_enrich,uresult,uresultA,uresultc)
    
    write.csv(denovo_enrich,file=fout)
    
    denovo<-gene_exp_Heart(denovo,sets$heart_exp)
    denovo<-gene_exp_Lung(denovo,sets$lung_exp)
    denovo<-gene_zscore(denovo,exac)
    denovo$solved=0;
    denovo$solved[which(denovo$ProbandName%in%solved_samples)]=1
    out_nms<-c("ProbandName","ExonicFunc.refGene","GeneName","Transcript","Exon","NucleotideChange","ProteinChange","REVEL","MCAP","MetaSVM_pred","Polyphen2_HVAR_pred","CADD13_PHRED","ExAC_ALL","pLI","lof_z","mis_z","HEART_EXP","LUNG_EXP","AAChange.refGene","solved")
    write.csv(denovo[,intersect(out_nms,names(denovo))],
              fout2,row.names = F)

}


fasso  <-  "PAH/Result/Data/source/PAH_associated11-13.txt"
fsolved<-"PAH/PAH-CHD/Data/PAH-CHD.solved.0608.csv" #PAH/Result/Data/source/PAH_ped_adult_All.known.csv
flist<-"PAH/Result/Data/source/PAH-CHD-list.csv"  #"PAH/Result/Data/source/PAH_ped_vs_adult_list.csv"
fdenovo<-"PAH/Result/Data/source/VCX_denovo_0819.txt"
fped<-"PAH/Result/Data/source/VCX_Control.ped"
fout<-"PAH/Result/Data/output/PAH_CHD_denovo.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_CHD_denovos.csv"
case.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
weslist<-case.list[grep("wes",case.list$SEQUENCING,ignore.case = T),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)


fasso  <-  "PAH/Result/Data/source/PAH_associated11-13.txt"
fsolved<-"PAH/Result/Data/source/PAH_ped_adult_All.known.csv"
flist<-"PAH/Result/Data/source/PAH_pediatric_adult_list.csv"
fdenovo<-"PAH/Result/Data/source/VCX_denovo_0819.txt"
fped<-"PAH/Result/Data/source/VCX_Control.ped"
fout<-"PAH/Result/Data/output/PAH_ped_adult.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_denovos.csv"
case.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
weslist<-case.list[grep("wes",case.list$SEQUENCING,ignore.case = T),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)




fasso  <-  "PAH/Result/Data/source/PAH_associated11-13.txt"
fsolved<-"PAH/Result/Data/source/PAH_ped_adult_All.known.csv"
flist<-"PAH/Result/Data/source/PAH_pediatric_adult_list.csv"
fdenovo<-"PAH/Result/Data/source/VCX_denovo_0819.txt"
fped<-"PAH/Result/Data/source/VCX_Control.ped"

case.list<-read.csv(flist,header = 1,stringsAsFactors = F,strip.white = T)
case.list<-case.list[which(case.list$Disease=="IPAH"),]


fout<-"PAH/Result/Data/output/PAH_ped_adult_IPAH.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_IPAH_denovos.csv"
weslist<-case.list[grep("wes",case.list$SEQUENCING,ignore.case = T),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)

fout<-"PAH/Result/Data/output/PAH_ped_adult_let5.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_let5_denovos.csv"
weslist<-case.list[intersect(grep("wes",case.list$SEQUENCING,ignore.case = T),which(case.list$Age_FT <=5 & case.list$Age_FT>-0.1) ),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)



fout<-"PAH/Result/Data/output/PAH_ped_adult_gt5.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_gt5_denovos.csv"
weslist<-case.list[intersect(grep("wes",case.list$SEQUENCING,ignore.case = T),which(case.list$Age_FT > 5 & case.list$TYPE=="child PAH") ),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)


fout<-"PAH/Result/Data/output/PAH_ped_adul_female_ped.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_female_ped_denovos.csv"
weslist<-case.list[intersect(grep("wes",case.list$SEQUENCING,ignore.case = T),which(case.list$TYPE=="child PAH" & case.list$Gender=="F") ),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)


fout<-"PAH/Result/Data/output/PAH_ped_adul_male_ped.enrich.csv"
fout2<-"PAH/Result/Data/output/PAH_ped_adult_male_ped_denovos.csv"
weslist<-case.list[intersect(grep("wes",case.list$SEQUENCING,ignore.case = T),which(case.list$TYPE=="child PAH" & case.list$Gender=="M") ),]
process(fdenovo,weslist,fasso,fsolved,fped,fout,fout2)

