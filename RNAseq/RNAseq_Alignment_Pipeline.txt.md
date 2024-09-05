

```bash
#step1:Check the quality of fastq files with FASTQC
fastqc -t 4 *.gz
mulitqc ./

#step2:Low-quality bases of FASTQ were removed by Trim_Galore
ls *1.fq.gz >1
ls *2.fq.gz >2
paste 1 2 > config
cat config | while read id; do arr=(${id});fq1=${arr[0]};fq2=${arr[1]};trim_galore -q 25 --phred33 --length 60 --stringency 15 --paired -o ./clean $fq1 $fq2;done #The example analysis shown
 
#step3:FASTQ files were aligned to the mouse genome mm10 using HISAT2
ls *gz |while read id; do hisat2 -p 8 -x /home/data_HD_1/reference/index/hisat2/human/hg19/genome -1 ${id}_1.fq.gz -2 ${id}_2.fq.gz -S ${id}_RNA-Seq.hisat2.sam; done #The example analysis shown 

ls *.sam | while read id; do (samtools sort -O bam -@ 12 -o $(basename ${id} .sam).bam ${id}); done
ls *.bam | while read id; do (samtools sort ${id} -o $(basename ${id} .bam)_sorted.bam); done
ls *_sorted.bam | while read id; do samtools index ${id}; done

#step4:FeatureCount software was used to calculate gene expression from aligned bam files
gtf="/home/wangfengsheng/data_HD_1/reference/gtf/human/gencode.v35lift37.annotation.gff3"
featureCounts -T 8 -p -t exon -g gene_id -a $gtf -o all.id.txt *sorted.bam 1> counts.id.log 2>&1 #The example analysis shown

```

DESeq2

```R
data<-read.table("all.id.txt",header=T,sep="\t")
data<-data[c(1,3:6)]
data<-na.omit(data)
library(dplyr)

rownames(data)<-data$Geneid
data<-data[-1]

colnames(data)
condition <- factor(c(rep("treat",2),rep("control",2)), levels = c("treat","control"))
condition

colData=data.frame(condition,row.names = colnames(data))

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c('condition', levels(colData$condition)))

res <- data.frame(res)
library(dplyr)
res_DEG<-filter(res,pvalue<0.05)
res_DEG_Up<-filter(res,pvalue<0.05 &log2FoldChange>0)
res_DEG_Down<-filter(res,pvalue<0.05 &log2FoldChange<0)


res$ENSG<-rownames(res)
Symbol<-read.csv("gProfiler_hsapiens_2022-10-29 下午7-50-05.csv")
res<-inner_join(res,Symbol,by="ENSG")
write.csv(res,"DEG_All.csv")
```

