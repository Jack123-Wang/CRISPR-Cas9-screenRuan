# step1：Multicox（cancer type+age+FER） ----------------------------------------------------------------
library("survival")
library("survminer")

###
data<-read.table("input1.txt",header=T,sep='\t')
multi_res.cox <- coxph(Surv(OS.time, OS) ~ FER+type+age, data = data)

x <- summary(multi_res.cox)


out_multi <- data.frame()
out_multi <- cbind(
  coef=x$coefficients[,"coef"],
  HR=x$conf.int[,"exp(coef)"],
  HR.95L=x$conf.int[,"lower .95"],
  HR.95H=x$conf.int[,"upper .95"],
  pvalue=x$coefficients[,"Pr(>|z|)"])

out_multi <- as.data.frame(cbind(id=row.names(out_multi),out_multi)) 
out_multi

out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
#write.csv(out_multi,"out_multi.output_cancer+age.csv")

# step2：Multicox（cancer type+FER） -------------------------------------------
multi_res.cox <- coxph(Surv(OS.time, OS) ~ FER+type, data = data)

x <- summary(multi_res.cox)

out_multi <- data.frame()
out_multi <- cbind(
  coef=x$coefficients[,"coef"],
  HR=x$conf.int[,"exp(coef)"],
  HR.95L=x$conf.int[,"lower .95"],
  HR.95H=x$conf.int[,"upper .95"],
  pvalue=x$coefficients[,"Pr(>|z|)"])

out_multi <- as.data.frame(cbind(id=row.names(out_multi),out_multi)) 
out_multi


out_multi[,2:ncol(out_multi)] <- as.numeric(unlist(out_multi[,2:ncol(out_multi)]))
#write.csv(out_multi,"out_multi.output_type.csv")


# step3：Figure(cancer type+age+FER) ---------------------------------------------------------------

out_multi<-read.csv("out_multi.output_cancer+age.csv")

out_multi
hz <- paste(round(out_multi$HR,3),
            "(",round(out_multi$HR.95L,3),
            "-",round(out_multi$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Characteristics",out_multi$id),
                   c(NA,"Coefficient",round(out_multi$coef,3)),
                   c(NA,"P value",ifelse(out_multi$pvalue<0.001,"P < 0.001",round(out_multi$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_multi$HR),
           lower=c(NA,NA,out_multi$HR.95L), 
           upper=c(NA,NA,out_multi$HR.95H), 
           boxsize=0.3,lwd.ci=2,   #
           ci.vertices.height = 0.08,ci.vertices=TRUE, #
           zero=1,lwd.zero=1,      #
           colgap=unit(5,"mm"),    #
           xticks = c(0,1,3,10), #
           is.summary=c(T,T,T,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,T),
           lwd.xaxis=1,            #
           lineheight = unit(0.9,"cm"), 
           graphwidth = unit(.1,"npc"), #
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), #
                           "40" = gpar(lwd=2, col="black")),#
           mar=unit(rep(0.5, times = 4), "cm"),#
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")



# step4:Figure（cancer type+FER） --------------------------------------------------

out_multi<-read.csv("out_multi.output_type.csv")

out_multi
hz <- paste(round(out_multi$HR,3),
            "(",round(out_multi$HR.95L,3),
            "-",round(out_multi$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Characteristics",out_multi$id),
                   c(NA,"Coefficient",round(out_multi$coef,3)),
                   c(NA,"P value",ifelse(out_multi$pvalue<0.001,"P < 0.001",round(out_multi$pvalue,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))


library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=3,  
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,out_multi$HR),
           lower=c(NA,NA,out_multi$HR.95L), 
           upper=c(NA,NA,out_multi$HR.95H), 
           boxsize=0.3,lwd.ci=2,   
           ci.vertices.height = 0.08,ci.vertices=TRUE, 
           zero=1,lwd.zero=1,     
           colgap=unit(5,"mm"),    
           xticks = c(0,1,3,10), 
           is.summary=c(T,T,T,F,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,F,T),
           lwd.xaxis=1,         
           lineheight = unit(0.9,"cm"), 
           graphwidth = unit(.1,"npc"), 
           cex=0.9, fn.ci_norm = fpDrawCircleCI, 
           hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                           "3" = gpar(lwd=2, col="black"), 
                           "38" = gpar(lwd=2, col="black")),
           mar=unit(rep(0.5, times = 4), "cm"),
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio")
