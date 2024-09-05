library(TCGAbiolinks)
projects <- TCGAbiolinks::getGDCprojects()$project_id 
projects <- projects[grepl("^TCGA", projects, perl = TRUE)]
projects

query <- GDCquery(
  project =projects,
  data.category = "Simple Nucleotide Variation",
  data.type="Masked Somatic Mutation",
  access="open"
)

GDCdownload(query, method = "api", files.per.chunk = 100) 

GDCdownload(query)

library(data.table)
library(maftools)
for( i in 1:33){
  projects[i]
  dir.path<- paste("GDCdata/",projects[i],sep="")
  all.maf<-list.files(path=dir.path,pattern=".gz",
                      full.names = T,recursive = T)
  maf.list<-lapply(all.maf,fread,
                   skip=7,
                   sep='\t',
                   header=T)
  maf.merge<-do.call(rbind,maf.list)
  maf2<-maf.merge[,c(5,6,16)]
  maf2$HGVsc<-maf.merge$HGVSc
  maf2$Variant_Type<-maf.merge$Variant_Type
  maf2$ref<-maf.merge$Reference_Allele #新的ref？
  table(maf2$Variant_Type)
  
  maf2<-filter(maf2,Variant_Type=="SNP")
  library(stringr)
  maf2$HGVsc<-str_sub(maf2$HGVsc,-3)
  #maf2$ref<-word(maf2$HGVsc,1,sep=fixed('>'))
  maf2$mut<-word(maf2$HGVsc,2,sep=fixed('>'))
  maf2<-maf2[,c(3,1,2,6,7)]
  colnames(maf2)<-c("sampleID","chr","pos","ref","mut")
  maf2$chr<-str_replace(maf2$chr,"chr","")
  maf2$chr<-as.numeric(maf2$chr)
  maf2<-na.omit(maf2)
  
  #hg38-hg19 convert
  maf2_output<-maf2[,c(2,3,3,1,4,5)]
  maf2_output$chr<-paste("chr",maf2_output$chr,sep="")
  maf2_output<-na.omit(maf2_output)
  write.table(maf2_output,paste(projects[i],"-maf_hg38.bed",sep=""),sep='\t',row.names=F,col.names = F,quote = F)
  
}

###

total_maftools<-NA
for( i in 24:33){
  projects[i]
  dir.path<- paste("GDCdata/",projects[i],sep="")
  all.maf<-list.files(path=dir.path,pattern=".gz",
                      full.names = T,recursive = T)
  maf.list<-lapply(all.maf,fread,
                   skip=7,
                   sep='\t',
                   header=T)
  maf.merge<-do.call(rbind,maf.list)
  maf1<-read.maf(maf.merge)
  
  hnsc.sig=oncodrive(maf1,AACol="HGVSp_Short",minMut = 5,pvalMethod = "zscore")
  hnsc.sig_FER<-filter(hnsc.sig,Hugo_Symbol=="FER")
  hnsc.sig_FER$Type<-projects[i]
  total_maftools<-rbind(total_maftools,hnsc.sig_FER,fill=T)
}
#liftOver  COAD-maf_hg38.bed hg38ToHg19.over.chain COAD-maf_hg19.bed unmapped.bed
#ls *bed |while read id;do liftOver $id hg38ToHg19.over.chain $(basename ${id} hg38.bed)hg19.bed unmapped.bed;done

write.csv(total_maftools,"maftools_oncodrive.FER.csv")

# step2 dndSNV-------------------------------------------------------------------
#library(devtools)
#install_github("im3sanger/dndscv")
#library(devtools)
#install_github("im3sanger/dndscv")
library(dndscv)
total<-NA
for ( i in 2:31){
  all.bed<-list.files(path=".",pattern="hg19.bed",
                      full.names = F,recursive = T)
  bed<-read.table(all.bed[i])
  colnames(bed)<-c("chr","pos","end","sampleID","ref",'mut')
  bed<-bed[-3]
  bed<-bed[,c(3,1,2,4,5)]
  bed$chr<-str_replace(bed$chr,"chr","")
  bed$chr<-as.numeric(bed$chr)
  dndsout = dndscv(bed)
  sel_cv = dndsout$sel_cv
  sel_cv_FER<-filter(sel_cv,gene_name=="FER")
  sel_cv_FER$Type<-str_replace(all.bed[i],"-maf_hg19.bed","")
  total<-rbind(total,sel_cv_FER)
}
write.csv(total,"dndscv_FER.result.csv")




