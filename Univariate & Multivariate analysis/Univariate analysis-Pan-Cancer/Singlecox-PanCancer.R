
library("survival")
library("survminer")

library(readxl)


data<-read.table("input1.txt",header=T,sep='\t')

res.cox1 <- coxph(Surv(OS.time, OS) ~ FER, data = data)
x1<-summary(res.cox1) 


out_multi1 <- data.frame()
out_multi1 <- cbind(
  coef=x1$coefficients[,"coef"],
  HR=x1$conf.int[,"exp(coef)"],
  HR.95L=x1$conf.int[,"lower .95"],
  HR.95H=x1$conf.int[,"upper .95"],
  pvalue=x1$coefficients[,"Pr(>|z|)"])

out_multi1 <- as.data.frame(cbind(id="FER",out_multi1)) 
out_multi1


#write.csv(out_multi1,"out_signle.output.csv")
out_single<-read.csv("out_signle.output.csv")
hz <- paste(round(out_single$HR,3),
            "(",round(out_single$HR.95L,3),
            "-",round(out_single$HR.95H,3),")",sep = "")

tabletext <- cbind(c(NA,"Characteristics",out_single$id),
                   c(NA,"Coefficient",round(out_single$coef,3)),
                   c(NA,"P value",ifelse(out_single$pvalue<0.001,"P < 0.001",round(out_single$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_single$HR),
           lower=c(NA,NA,out_single$HR.95L), 
           upper=c(NA,NA,out_single$HR.95H), #
           boxsize=0.1,lwd.ci=2,  
           ci.vertices.height = 0.08,ci.vertices=TRUE, #
           zero=1,lwd.zero=1,      
           colgap=unit(5,"mm"),    
           xticks = c(0,0.5,1), 
           is.summary=c(T,T,F,F),
           lwd.xaxis=1,            
           lineheight = unit(0.9,"cm"), #
           graphwidth = unit(.1,"npc"), #
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #
                           "4" = gpar(lwd=2, col="black")),#
           mar=unit(rep(0.5, times = 4), "cm"),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")



