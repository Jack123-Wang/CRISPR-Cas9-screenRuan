
# step1:reading files --------------------------------------------------------------
library(stringr)
files<-list.files(pattern="tsv")
data<-NA
for (i in 1:33){
  a<-read.table(files[i],sep='\t',header=T)
  a$type<-word(files[i],1,sep=fixed('.'))
  data<-rbind(data,a)
}

data<-unique(data)
data<-na.omit(data)
data<-as.data.frame(data)
# step2:Univariate  analysis-33 cancers --------------------------------------------------
library(dplyr)
library("survival")
library("survminer")
total<-NA
for (i in 1:length(unique(data$type))){
  a<-unique(data$type)[i]
  b<-filter(data,type==a)
  res.cox <- coxph(Surv(OS.time, OS) ~ FER, data = b)
  x1<-summary(res.cox) 
  
  out_multi1 <- data.frame()
  out_multi1 <- cbind(
    coef=x1$coefficients[,"coef"],
    HR=x1$conf.int[,"exp(coef)"],
    HR.95L=x1$conf.int[,"lower .95"],
    HR.95H=x1$conf.int[,"upper .95"],
    pvalue=x1$coefficients[,"Pr(>|z|)"],
    type=a)
  total<-rbind(total,out_multi1)
}


total<-as.data.frame(total)
total<-na.omit(total)

total$pvalue<-as.numeric(total$pvalue)

#p<0.05
total_fil<-filter(total,pvalue<0.05)

#write.csv(total,"TCGA-FER-univariate.cox.results.csv")

# step3:Figure ---------------------------------------------------------------

out_multi<-read.csv("TCGA-FER-univariate.cox.results.csv")

out_multi
hz <- paste(round(out_multi$HR,3),
            "(",round(out_multi$HR.95L,3),
            "-",round(out_multi$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Characteristics",out_multi$id),
                   c(NA,"Coefficient",round(out_multi$coef,3)),
                   c(NA,"P value",ifelse(out_multi$pvalue<0.001,"P < 0.001",round(out_multi$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz),
                   c(NA,"Cancer Type",out_multi$type))


library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_multi$HR),
           lower=c(NA,NA,out_multi$HR.95L), 
           upper=c(NA,NA,out_multi$HR.95H),
           boxsize=0.3,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, #
           zero=1,lwd.zero=1,      
           colgap=unit(5,"mm"),   
           xticks = c(0,1,3),
           is.summary=c(T,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,T),
           lwd.xaxis=1,           
           lineheight = unit(0.9,"cm"), 
           graphwidth = unit(.1,"npc"), #
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #
                           "36" = gpar(lwd=2, col="black")),#
           mar=unit(rep(0.5, times = 4), "cm"),#
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")


