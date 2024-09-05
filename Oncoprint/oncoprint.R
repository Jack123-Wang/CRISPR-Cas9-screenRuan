# Organize the mutation matrix ------------------------------------------------------------------
data<-read.table("alterations_across_samples.tsv",sep='\t',header=T)
genelist<-colnames(data)[c(5:45)]
library(dplyr)
library(stringr)
k00<-NA
for(i in 1:1){
  a<-genelist[i]
  b<-cbind("Sample.ID",paste(a,"..MUT",sep=""),paste(a,"..FUSION",sep=""),paste(a,"..AMP",sep=""))
  b<-as.character(b)
  c<-select(data,b)
  c[,2]<-str_replace(c[,2],"no alteration","")
  c[,2]<-str_replace(c[,2],"not profiled","")
  c[,2]<-str_replace(c[,2],"[A-Z]","MUT;")
  c[,2]<-str_sub(c[,2],start=1,end=4)
  
  c[,3]<-str_replace(c[,3],"no alteration","")
  c[,3]<-str_replace(c[,3],"not profiled","")
  c[,3]<-str_replace(c[,3],"[A-Z]","Fusion;")
  c[,3]<-str_sub(c[,3],start=1,end=7)
  
  c[,4]<-str_replace(c[,4],"no alteration","")
  c[,4]<-str_replace(c[,4],"not profiled","")
  c[,4]<-str_replace(c[,4],"AMP","AMP;")
  c[,4]<-str_replace(c[,4],"HOMDEL","HOMDEL;")
  
  c$NF2..total<-""
  k0<-NA
  for(j in 1:nrow(c)){
    k<-c[j,]
    k[1,5]<-paste(k[1,2],k[1,3],k[1,4],"",sep="")
    k0<-rbind(k0,k)
  }
  k0<-na.omit(k0)
  k0<-k0[c(1,5)]
  colnames(k0)[2]<-a
  k00<-cbind(k00,k0)
}


for(i in 3:length(genelist)){
  a<-genelist[i]
  b<-cbind("Sample.ID",paste(a,"..MUT",sep=""),paste(a,"..FUSION",sep=""),paste(a,"..AMP",sep=""))
  b<-as.character(b)
  c<-select(data,b)
  c[,2]<-str_replace(c[,2],"no alteration","")
  c[,2]<-str_replace(c[,2],"not profiled","")
  c[,2]<-str_replace(c[,2],"[A-Z]","MUT;")
  c[,2]<-str_sub(c[,2],start=1,end=4)
  
  c[,3]<-str_replace(c[,3],"no alteration","")
  c[,3]<-str_replace(c[,3],"not profiled","")
  c[,3]<-str_replace(c[,3],"[A-Z]","Fusion;")
  c[,3]<-str_sub(c[,3],start=1,end=7)
  
  c[,4]<-str_replace(c[,4],"no alteration","")
  c[,4]<-str_replace(c[,4],"not profiled","")
  c[,4]<-str_replace(c[,4],"AMP","AMP;")
  c[,4]<-str_replace(c[,4],"HOMDEL","HOMDEL;")
  
  c$NF2..total<-""
  k0<-NA
  for(j in 1:nrow(c)){
    k<-c[j,]
    k[1,5]<-paste(k[1,2],k[1,3],k[1,4],"",sep="")
    k0<-rbind(k0,k)
  }
  k0<-na.omit(k0)
  k0<-k0[c(1,5)]
  colnames(k0)[2]<-a
  k00<-inner_join(k00,k0,by="Sample.ID")
}

#write.csv(k00,"mat.csv")
#excel 去除"(driver)"


# Determine which are change groups and which are non change groups --------------------------------------------------------

alter<-read.table("altered_samples.txt")
alter$V2<-"alter"
non_alter<-read.table("unaltered_samples.txt")
non_alter$V2<-"non_alter"

total<-rbind(alter,non_alter)
colnames(total)<-c("Sample.ID","Type")

# Compare the differences between the top 20 genes in two groups ---------------------------------------------------------

mat<-read.csv("mat.csv")
mat<-inner_join(mat,total,by="Sample.ID")
mat_alter<-filter(mat,Type=="alter")
mat_nonalter<-filter(mat,Type=="non_alter")

#TP53 alter 47.09%
test1<-as.data.frame(table(mat_alter$TP53!=""))
test1$Per<-test1$Freq/sum(test1$Freq)
test1

#TP53 non_alter 30.35%
test2<-as.data.frame(table(mat_nonalter$TP53!=""))
test2$Per<-test2$Freq/sum(test2$Freq)
test2

#TTN alter 42.94%
test1<-as.data.frame(table(mat_alter$TTN!=""))
test1$Per<-test1$Freq/sum(test1$Freq)
test1

#TTN non_alter 22.05%
test2<-as.data.frame(table(mat_nonalter$TTN!=""))
test2$Per<-test2$Freq/sum(test2$Freq)
test2

#MUC16 alter 28.69%
test1<-as.data.frame(table(mat_alter$MUC16!=""))
test1$Per<-test1$Freq/sum(test1$Freq)
test1

#MUC16 non_alter 14.72%
test2<-as.data.frame(table(mat_nonalter$MUC16!=""))
test2$Per<-test2$Freq/sum(test2$Freq)
test2


# Survival --------------------------------------------------------------

alter<-read.table("altered_samples.txt")
alter$V2<-"alter"
non_alter<-read.table("unaltered_samples.txt")
non_alter$V2<-"non_alter"

total<-rbind(alter,non_alter)
colnames(total)<-c("Sample.ID","Type")

survival<-read.table("TCGAPanCancer-Survival.tsv",sep='\t',header=T)
survival<-inner_join(survival,total,by="Sample.ID")

library(survival)
library(survminer)
library(dplyr)


fit<-survfit(Surv(OS.time, OS) ~ Type, data = survival)
print(fit)
ggsurvplot(fit,pval = TRUE, 
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           risk.table.heght=1,
           #surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           #legend.labs=c("A","B","C"),
           palette = c("#BC4835", "#4479C1"),
           data=survival)


# Compare TMB/Neoantigens/Heterogeneity between two groups ------------------------------------------------------------
meta<-read.csv("PATIENT_DATA_oncoprint.csv",header=T)
colnames(meta)[1]<-"Sample.ID1"
total$Sample.ID1<-str_sub(total$Sample.ID,start = 1,end=12)
meta<-inner_join(meta,total,by="Sample.ID1")  

meta_alter<-filter(meta,Type=="alter")
meta_nonalter<-filter(meta,Type=="non_alter")

#TMB alter 14.63498
meta_alter_TMB<-meta_alter$TMB
meta_alter_TMB<-na.omit(meta_alter_TMB)
mean(meta_alter_TMB)
meta_alter_TMB<-as.data.frame(meta_alter_TMB)
meta_alter_TMB$type<-"alter"
colnames(meta_alter_TMB)<-c("TMB","type")
#TMB non_alter 3.25467
meta_nonalter_TMB<-meta_nonalter$TMB
meta_nonalter_TMB<-na.omit(meta_nonalter_TMB)
mean(meta_nonalter_TMB)
meta_nonalter_TMB<-as.data.frame(meta_nonalter_TMB)
meta_nonalter_TMB$type<-"non-alter"
colnames(meta_nonalter_TMB)<-c("TMB","type")
#

meta_TMB<-rbind(meta_alter_TMB,meta_nonalter_TMB)
meta_TMB$TMB<-log10(meta_TMB$TMB+0.01)
library(ggplot2)
library(ggsignif)
p<-ggplot(meta_TMB,aes(x=type,y=TMB,fill=type))
p+geom_violin()+
  geom_boxplot(width=0.25)+
  theme_classic()+#
  theme(panel.grid.major = element_line(color = "gray60", size = .5,linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray60", size = .5,linetype = "dashed"))+
  geom_jitter(width = 0.2,size=0.01,alpha=0.1)+
  scale_fill_manual(values=c("#A84D3A","#4C72AF"))+
  labs(x="")+
  theme(axis.line.x=element_line(linetype=,color="black",size=1),#
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 22),  #
        axis.ticks.y=element_line(color="black",size=1,lineend = 22))+
  geom_signif(comparisons = list(c("alter","non-alter")),test="t.test")

#Neoantigens alter 185.5896
meta_alter_Neoantigens<-meta_alter$SNV.Neoantigens
meta_alter_Neoantigens<-na.omit(meta_alter_Neoantigens)
mean(meta_alter_Neoantigens)
meta_alter_Neoantigens<-as.data.frame(meta_alter_Neoantigens)
meta_alter_Neoantigens$type<-"alter"
colnames(meta_alter_Neoantigens)<-c("Neoantigens","type")
#Neoantigens non_alter 36.8935
meta_nonalter_Neoantigens<-meta_nonalter$SNV.Neoantigens
meta_nonalter_Neoantigens<-na.omit(meta_nonalter_Neoantigens)
mean(meta_nonalter_Neoantigens)
meta_nonalter_Neoantigens<-as.data.frame(meta_nonalter_Neoantigens)
meta_nonalter_Neoantigens$type<-"non-alter"
colnames(meta_nonalter_Neoantigens)<-c("Neoantigens","type")

meta_Neoantigens<-rbind(meta_alter_Neoantigens,meta_nonalter_Neoantigens)
meta_Neoantigens$Neoantigens<-log10(meta_Neoantigens$Neoantigens+0.01)
library(ggplot2)
library(ggsignif)
p<-ggplot(meta_Neoantigens,aes(x=type,y=Neoantigens,fill=type))
p+geom_violin()+
  geom_boxplot(width=0.25)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "gray60", size = .5,linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray60", size = .5,linetype = "dashed"))+
  geom_jitter(width = 0.2,size=0.01,alpha=0.1)+
  scale_fill_manual(values=c("#A84D3A","#4C72AF"))+
  labs(x="")+
  theme(axis.line.x=element_line(linetype=,color="black",size=1),#
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 22),  #
        axis.ticks.y=element_line(color="black",size=1,lineend = 22))+
  geom_signif(comparisons = list(c("alter","non-alter")),test="t.test")

#Neoantigens alter 0.1561835
meta_alter_Heterogeneity<-meta_alter$Intratumor.Heterogeneity
meta_alter_Heterogeneity<-na.omit(meta_alter_Heterogeneity)
mean(meta_alter_Heterogeneity)
meta_alter_Heterogeneity<-as.data.frame(meta_alter_Heterogeneity)
meta_alter_Heterogeneity$type<-"alter"
colnames(meta_alter_Heterogeneity)<-c("Heterogeneity","type")

#Neoantigens non_alter 0.1264779
meta_nonalter_Heterogeneity<-meta_nonalter$Intratumor.Heterogeneity
meta_nonalter_Heterogeneity<-na.omit(meta_nonalter_Heterogeneity)
mean(meta_nonalter_Heterogeneity)
meta_nonalter_Heterogeneity<-as.data.frame(meta_nonalter_Heterogeneity)
meta_nonalter_Heterogeneity$type<-"non-alter"
colnames(meta_nonalter_Heterogeneity)<-c("Heterogeneity","type")

meta_Heterogeneity<-rbind(meta_alter_Heterogeneity,meta_nonalter_Heterogeneity)

library(ggplot2)
library(ggsignif)
p<-ggplot(meta_Heterogeneity,aes(x=type,y=Heterogeneity,fill=type))
p+geom_violin()+
  geom_boxplot(width=0.25)+
  theme_classic()+
  theme(panel.grid.major = element_line(color = "gray60", size = .5,linetype = "dashed"),
        panel.grid.minor = element_line(color = "gray60", size = .5,linetype = "dashed"))+
  geom_jitter(width = 0.2,size=0.01,alpha=0.1)+
  scale_fill_manual(values=c("#A84D3A","#4C72AF"))+
  labs(x="")+
  theme(axis.line.x=element_line(linetype=,color="black",size=1),
        axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 22),  
        axis.ticks.y=element_line(color="black",size=1,lineend = 22))+
  geom_signif(comparisons = list(c("alter","non-alter")),test="t.test")

# oncoprint ---------------------------------------------------------------
library(devtools)
#install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
mat<-read.csv("mat1.csv")

mat<-t(as.matrix(mat))
colnames(mat)<-mat[1,]
mat<-mat[-1,]
mat[1:10, 1:10]
mat1<-mat[c(1:20),]

col = c("FUSION"="#FDF7AF", "AMP" = "#C22E1E", "MUT" = "#31723F","HOMDEL" = "#2C689D")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.0002, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#E6E6E5", col = NA))
  },
  FUSION = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.0002, "pt"), h*0.33, 
              gp = gpar(fill = col["FUSION"], col = NA))
  },

  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.0002, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.0002, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["MUT"], col = NA),
              )
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.0002, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  }
)


# just for demonstration
alter_fun = list(
  background = alter_graphic("rect", fill = "#E6E6E5"),   
  FUSION = alter_graphic("rect", height = 0.33,fill = col["FUSION"]),
  AMP = alter_graphic("rect", fill = col["AMP"]),
  MUT = alter_graphic("rect", fill = col["MUT"]),
  HOMDEL = alter_graphic("rect", fill = col["HOMDEL"])
)
column_title = "OncoPrint for TCGA PanCancer, genes in 20 candidate genes"
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL", "AMP", "MUT","FUSION"), 
                            labels = c("Deep deletion", "Amplification", "Mutation","Fusion"))


oncoPrint(mat1,
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            cbar = anno_oncoprint_barplot()          ),
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          remove_empty_columns = F, remove_empty_rows = T)


# Other Top20 genes -------------------------------------------------------
mat<-read.csv("mat1.csv")

mat<-t(as.matrix(mat))
colnames(mat)<-mat[1,]
mat<-mat[-1,]
mat[1:10, 1:10]
mat1<-mat[c(21:40),]

#metadata
meta<-read.csv("PATIENT_DATA_oncoprint.csv",header=T)
colnames(meta)[1]<-"Sample.ID1"

alter<-read.table("altered_samples.txt")
alter$V2<-"alter"
non_alter<-read.table("unaltered_samples.txt")
non_alter$V2<-"non_alter"

total<-rbind(alter,non_alter)
colnames(total)<-c("Sample.ID","Type")
library(stringr)
total$Sample.ID1<-str_sub(total$Sample.ID,start=1,end=12)

library(dplyr)
meta<-inner_join(meta,total,by="Sample.ID1")

mat_col<-as.data.frame(colnames(mat1))
mat_col$order<-c(1:10967)
colnames(mat_col)[1]<-"Sample.ID"
mat_col<-inner_join(mat_col,meta,by="Sample.ID")

group<-mat_col[c(1,14)]
group_alter<-filter(group,Type=="alter")
group_nonalter<-filter(group,Type=="non_alter")

mat1_alter<-subset(mat1,select=group_alter$Sample.ID)
mat1_nonalter<-subset(mat1,select=group_nonalter$Sample.ID)

mat_col_alter<-filter(mat_col,Type=="alter")
mat_col_nonalter<-filter(mat_col,Type=="non_alter")

library(circlize)
library(ComplexHeatmap)

oncoPrint(mat1_alter,
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            SampleType=mat_col_alter$Sample.Type, 
            Sex=mat_col_alter$Sex, 
            Age=mat_col_alter$Diagnosis.Age, 
            OS=mat_col_alter$Overall.Survival.Status, 
            OS.Time=anno_barplot(mat_col_alter$Overall.Survival..Months.,border = F, height = unit(20, "pt")),
            TMB=anno_barplot(mat_col_alter$TMB,border = F, height = unit(20, "pt")),
            Neoantigens=anno_barplot(mat_col_alter$SNV.Neoantigens,border = F, height = unit(20, "pt")),
            Intratumor.Heterogeneity=anno_barplot(mat_col_alter$Intratumor.Heterogeneity,border = F, height = unit(20, "pt")),
            na_col = "black",
            simple_anno_size = unit(7.5, "pt"), 
            col = list(
              Sex = c("Female" = "#F8DCCA", "Male" = "#E9A888"),
              SampleType = c("Primary"="#63C187","Mixed" = "#4FA06B", "Metastasis" = "#387A47", "Recurrence"="#2D5A3A"),
              OS=c("1:DECEASED"="#BABABA","0:LIVING"="#E0E0E0"),
              Age = colorRamp2(c(0, 100), c("white", "#62E18F"))
            )
          ),
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          remove_empty_columns = F, remove_empty_rows = F)

oncoPrint(mat1_nonalter,
          alter_fun = alter_fun, col = col, 
          top_annotation = HeatmapAnnotation(
            SampleType=mat_col_nonalter$Sample.Type, 
            Sex=mat_col_nonalter$Sex, 
            Age=mat_col_nonalter$Diagnosis.Age, 
            OS=mat_col_nonalter$Overall.Survival.Status, 
            OS.Time=anno_barplot(mat_col_nonalter$Overall.Survival..Months.,border = F, height = unit(20, "pt")),
            TMB=anno_barplot(mat_col_nonalter$TMB,border = F, height = unit(20, "pt")),
            Neoantigens=anno_barplot(mat_col_nonalter$SNV.Neoantigens,border = F, height = unit(20, "pt")),
            Intratumor.Heterogeneity=anno_barplot(mat_col_nonalter$Intratumor.Heterogeneity,border = F, height = unit(20, "pt")),
            na_col = "black",
            simple_anno_size = unit(7.5, "pt"), 
            col = list(
              Sex = c("Female" = "#F8DCCA", "Male" = "#E9A888"),
              SampleType = c("Primary"="#63C187","Mixed" = "#4FA06B", "Metastasis" = "#387A47", "Recurrence"="#2D5A3A"),
              OS=c("1:DECEASED"="#BABABA","0:LIVING"="#E0E0E0"),
              Age = colorRamp2(c(0, 100), c("white", "#62E18F"))
            )
          ),
          column_title = column_title, heatmap_legend_param = heatmap_legend_param,
          remove_empty_columns = F, remove_empty_rows = F)


