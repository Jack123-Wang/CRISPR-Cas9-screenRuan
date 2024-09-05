library(Seurat)
expression_matrix_1<-Read10X('../singlecell/NL1/')
expression_matrix_2<-Read10X('../singlecell/NL2/')
expression_matrix_3<-Read10X('../singlecell/NL3/')

cds1<-CreateSeuratObject(counts = expression_matrix_1, min.cells = 5, min.features = 168, project = 'N1')
cds2<-CreateSeuratObject(counts = expression_matrix_2, min.cells = 8, min.features = 168, project = 'N2')
cds3<-CreateSeuratObject(counts = expression_matrix_3, min.cells = 8, min.features = 168, project = 'N3')

cds1@meta.data$percent.mt<-PercentageFeatureSet(cds1, pattern = "^MT-")
cds2@meta.data$percent.mt<-PercentageFeatureSet(cds2, pattern = "^MT-")
cds3@meta.data$percent.mt<-PercentageFeatureSet(cds3, pattern = "^MT-")

VlnPlot(cds1, group.by = NULL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(cds2, group.by = NULL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(cds3, group.by = NULL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(cds1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cds1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(cds2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(cds2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
plot5 <- FeatureScatter(cds3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(cds3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot5 + plot6

cds1 <- subset(cds1, subset = nFeature_RNA > 0 & nFeature_RNA < 5000 & percent.mt < 5)
cds2 <- subset(cds2, subset = nFeature_RNA > 0 & nFeature_RNA < 6000 & percent.mt < 5)
cds3 <- subset(cds3, subset = nFeature_RNA > 0 & nFeature_RNA < 6000 & percent.mt < 5)

cds1 <- NormalizeData(cds1, normalization.method = "LogNormalize", scale.factor = 10000)
cds2 <- NormalizeData(cds2, normalization.method = "LogNormalize", scale.factor = 10000)
cds3 <- NormalizeData(cds3, normalization.method = "LogNormalize", scale.factor = 10000)

cds.combined <- merge(cds1, y = c(cds2,cds3), add.cell.ids = c("N1", "N2",'N3'), project = "N1+N2+N3")
cds.combined <- FindVariableFeatures(cds.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cds.combined)
cds.combined <- ScaleData(cds.combined, features = all.genes)
cds.combined <- RunPCA(cds.combined, features = VariableFeatures(object = cds.combined))
VizDimLoadings(cds.combined, dims = 1:2, nfeatures = 30, reduction = "pca")
PCAPlot(cds.combined, reduction = "pca")
ElbowPlot(cds.combined, ndims = 30, reduction = "pca")
#cds.combined <- FindNeighbors(cds.combined, dims = 1:10, reduction = "pca")
#cds.combined <- FindClusters(cds.combined, resolution = 0.1)
cds.combined <- FindNeighbors(cds.combined, dims = 1:12, reduction = "pca")
cds.combined <- FindClusters(cds.combined, resolution = 0.5)
cds.combined <- RunUMAP(cds.combined, features = VariableFeatures(object = cds.combined))
DimPlot(cds.combined, 
        group.by = 'seurat_clusters',
        label = F,
        cols = c('0'='#F68282','1'='#ff9a36','2'='#2FF18B','3'='#D4D915','4'='#31C53F','5'='#CCB1F1','6'='#B95FBB','7'='#28CECA',
                 '8'='#1FA195','9'='#25aff5','10'='#aeadb3','11'='#8080C0','12'='#A4DFF2','13'='#AC8F14',
                 '14'='#E6C122','15'='#4B4BF7','16'='#84C1FF','17'='#B8B8DC'),
        label.size = 1.5, 
        label.color = "black",
        label.box = F,
        repel = FALSE,
        pt.size = 1.0, 
        reduction = "umap")
FeaturePlot(cds.combined, features = c("FER"), reduction = "umap")

cds.combined_1 <- cds.combined[,cds.combined@meta.data$seurat_clusters %in% c(0:8)]
cds.combined_2 <- cds.combined[,cds.combined@meta.data$seurat_clusters %in% c(9:17)]
VlnPlot(cds.combined, group.by = 'seurat_clusters', features = c("FER"))

####Differential analysis####
exprset <- cds.combined@assays$RNA
exprset <- as.data.frame(exprset@data)
FER_high <- colnames(exprset)[exprset[rownames(exprset)=='FER',]>0]
FER_low <- colnames(exprset)[exprset[rownames(exprset)=='FER',]==0]
cds.combined@meta.data$FER_type <- ifelse(rownames(cds.combined@meta.data) %in% FER_high,'high','low')
Idents(cds.combined)="FER_type"
cluster.markers_FER <- FindMarkers(cds.combined, 
                                   ident.1 = c('high'), 
                                   ident.2 = c('low'), 
                                   min.pct = 0.1, 
                                   logfc.threshold = 0)
cluster.markers_FER$gene <- rownames(cluster.markers_FER)
cds.combined_diff <- cds.combined[rownames(cds.combined) %in% rownames(cluster.markers_FER),]

av <-AverageExpression(cds.combined_diff, assays = "RNA",group.by = "FER_type") 
results <- as.data.frame(av$RNA)
results$gene <- rownames(results)
merge <- merge(results,cluster.markers_FER,'gene')
merge <- merge[,c(1:4)]
write.csv(merge,'Differential analysis.csv',quote = F,row.names = T)

DimPlot(cds.combined, 
        group.by = 'FER_type',
        label = T, 
        label.size = 3, 
        label.color = "black",
        label.box = T,
        repel = FALSE,
        pt.size = 1.0, 
        reduction = "umap")


