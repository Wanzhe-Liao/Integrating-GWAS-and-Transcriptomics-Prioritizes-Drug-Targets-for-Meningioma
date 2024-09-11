library(harmony)
library(data.table)
#设定工作目录
setwd('./scRNA/')
setwd("C:/Users/84759/Desktop/胆管炎cholangitis")
load('Meningioma.RData')
scRNA <- readRDS("C:/Users/84759/Desktop/胆管炎cholangitis/Pipeline/scRNA_anno1.RDS")
#设定样本,CT为对照组，PR为疾病组
samples <- c('CT','PSC')

library(Seurat)
library(tidyverse)
library(vroom)


#读取数据
PSC1 <- Read10X(data.dir = './PSC1/')
PSC2 <- Read10X(data.dir = './PSC2/')
PSC3 <- Read10X(data.dir = './PSC3/')
PSC4 <- Read10X(data.dir = './PSC4/')
PSC5 <- Read10X(data.dir = './PSC5/')
PSC6 <- Read10X(data.dir = './PSC6/')
PSC7 <- Read10X(data.dir = './PSC7/')
PSC8 <- Read10X(data.dir = './PSC8/')


for (i in 11:22)
  {rownames(MSCi) <- gsub("hg19_","",rownames(MSCi))}
  

CT1 <- Read10X(data.dir = './CT1/')
CT2 <- Read10X(data.dir = './CT2/')
CT3 <- Read10X(data.dir = './CT3/')
CT4 <- Read10X(data.dir = './CT4/')
CT5 <- Read10X(data.dir = './CT5/')
CT6 <- Read10X(data.dir = './CT6/')
CT7 <- Read10X(data.dir = './CT7/')
CT8 <- Read10X(data.dir = './CT8/')
CT9 <- Read10X(data.dir = './CT9/')
CT10 <- Read10X(data.dir = './CT10/')
CT11 <- Read10X(data.dir = './CT11/')
CT12 <- Read10X(data.dir = './CT12/')
CT13 <- Read10X(data.dir = './CT13/')
CT14 <- Read10X(data.dir = './CT14/')
CT15 <- Read10X(data.dir = './CT15/')
CT16 <- Read10X(data.dir = './CT16/')
CT17 <- Read10X(data.dir = './CT17/')
CT18 <- Read10X(data.dir = './CT18/')
CT19 <- Read10X(data.dir = './CT19/')
CT20 <- Read10X(data.dir = './CT20/')
CT21 <- Read10X(data.dir = './CT21/')

#分别创建SeuratObject
scRNA1 <- CreateSeuratObject(counts = CT1, project = "CT")
scRNA2 <- CreateSeuratObject(counts = CT2, project = "CT")
scRNA3 <- CreateSeuratObject(counts = CT3, project = "CT")
scRNA4 <- CreateSeuratObject(counts = CT4, project = "CT")
scRNA5 <- CreateSeuratObject(counts = CT5, project = "CT")
scRNA6 <- CreateSeuratObject(counts = CT6, project = "CT")
scRNA7 <- CreateSeuratObject(counts = CT7, project = "CT")
scRNA8 <- CreateSeuratObject(counts = CT8, project = "CT")
scRNA9 <- CreateSeuratObject(counts = CT9, project = "CT")
scRNA10 <- CreateSeuratObject(counts = CT10, project = "CT")
scRNA11 <- CreateSeuratObject(counts = CT11, project = "CT")
scRNA12 <- CreateSeuratObject(counts = CT12, project = "CT")
scRNA13 <- CreateSeuratObject(counts = CT13, project = "CT")
scRNA14 <- CreateSeuratObject(counts = CT14, project = "CT")
scRNA15 <- CreateSeuratObject(counts = CT15, project = "CT")
scRNA16 <- CreateSeuratObject(counts = CT16, project = "CT")
scRNA17 <- CreateSeuratObject(counts = CT17, project = "CT")
scRNA18 <- CreateSeuratObject(counts = CT18, project = "CT")
scRNA19 <- CreateSeuratObject(counts = CT19, project = "CT")
scRNA20 <- CreateSeuratObject(counts = CT20, project = "CT")
scRNA21 <- CreateSeuratObject(counts = CT21, project = "CT")
scRNA22 <- CreateSeuratObject(counts = PSC1, project = "PSC")
scRNA23 <- CreateSeuratObject(counts = PSC2, project = "PSC")
scRNA24 <- CreateSeuratObject(counts = PSC3, project = "PSC")
scRNA25 <- CreateSeuratObject(counts = PSC4, project = "PSC")
scRNA26 <- CreateSeuratObject(counts = PSC5, project = "PSC")
scRNA27 <- CreateSeuratObject(counts = PSC6, project = "PSC")
scRNA28 <- CreateSeuratObject(counts = PSC7, project = "PSC")
scRNA29 <- CreateSeuratObject(counts = PSC8, project = "PSC")



scRNA10 <- CreateSeuratObject(counts = MSC10, project = "CT")
scRNA11 <- CreateSeuratObject(counts = MSC11, project = "BM")
scRNA12 <- CreateSeuratObject(counts = MSC12, project = "BM")
scRNA13 <- CreateSeuratObject(counts = MSC13, project = "BM")
scRNA14 <- CreateSeuratObject(counts = MSC14, project = "BM")
scRNA15 <- CreateSeuratObject(counts = MSC15, project = "BM")
scRNA16 <- CreateSeuratObject(counts = MSC16, project = "BM")
scRNA17 <- CreateSeuratObject(counts = MSC17, project = "BM")
scRNA18 <- CreateSeuratObject(counts = MSC18, project = "BM")
scRNA19 <- CreateSeuratObject(counts = MSC19, project = "BM")
scRNA20 <- CreateSeuratObject(counts = MSC20, project = "BM")
scRNA21 <- CreateSeuratObject(counts = MSC21, project = "BM")
scRNA22 <- CreateSeuratObject(counts = MSC22, project = "BM")

scRNA99 <- CreateSeuratObject(counts = MSC99, project = "BM")
scRNA100 <- CreateSeuratObject(counts = MSC100, project = "BM")




#利用merge函数把两个样本合在一起，并且给每个样本都加上各自的样本名，方便后续分析
sce.all = merge(scRNA1, y = c(scRNA2,scRNA3,scRNA4,scRNA5,scRNA6,scRNA7,scRNA8,scRNA9,scRNA10,scRNA11,scRNA12,scRNA13,scRNA14,scRNA15,scRNA16,scRNA17,scRNA18,scRNA19,scRNA20,scRNA21,scRNA22,scRNA23,scRNA24,scRNA25,scRNA26,scRNA27,scRNA28,scRNA29), 
                add.cell.ids = c("CT1", "CT2","CT3","CT4","CT5","CT6","CT7","CT8","CT9","CT10","CT11","CT12","CT13","CT14","CT15","CT16","CT17","CT18","CT19","CT20","CT21","PSC1","PSC2","PSC3","PSC4","PSC5","PSC6","PSC7","PSC8"),
                project = 'PSC', merge.data = TRUE)


sce.all = merge(scRNA1, y = c(scRNA2,scRNA3,scRNA4,scRNA22,scRNA23,scRNA24,scRNA25), 
                add.cell.ids = c("CT1", "CT2","CT3","CT4","PSC1","PSC2","PSC3","PSC4"),
                project = 'PSC', merge.data = TRUE)



sce.all = merge(scRNA11, y = c(scRNA12,scRNA13,scRNA14,scRNA15,scRNA16,scRNA17,scRNA18,scRNA19,scRNA20,scRNA21,scRNA22), 
                add.cell.ids = c("BM1", "BM2","BM3","BM4","BM5","BM6","BM7","BM8","BM9","BM10","BM11","BM12"),
                project = 'BM', merge.data = TRUE)

#metadata为样本信息，我们需要定义分组
sce.all@meta.data$patient=sce.all@meta.data$orig.ident
sce.all@meta.data$orig.ident=stringr::str_remove(sce.all@meta.data$orig.ident,'[0-9]')

# 开始流程
scRNA=sce.all
dir.create('Pipeline')
setwd('./Pipeline/')
dir.create('QC')
##计算质控指标
#计算细胞中核糖体基因比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(scRNA@active.ident))
#质控前
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
pearplot
ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 


##设置质控标准，按照你给我的文献来定，如果发现后面结果不太满意，可以调整一下阈值，但是影响不大的
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=200", "maxGene=2500", "pctMT=10",'pctHB=3'))
minGene=500
maxGene=20000
pctMT=10
pctHB=3

##1.数据质控
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)

col.num <- length(levels(scRNA@active.ident))
violin <-VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
pearplot
ggsave("QC/pearplot_after_qc.pdf", plot = pearplot, width = 12, height = 5) 
# 标准化
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
#2.降维和聚类#########################
library(Seurat)
library(tidyverse)
library(patchwork)
dir.create("cluster")
#高变基因3000个
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
plot

ggsave("cluster/VariableFeatures.pdf", plot = plot, width = 8, height = 6) 

#如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(scRNA)
#gc()
#scRNA <- ScaleData(scRNA, features = scale.genes)
##如果内存不够，可以只对高变基因进行标准化
#scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=40, reduction="pca") 
plotc <- plot1+plot2
plotc
ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 

# 选取平缓的elbow，不用更改
pc.num=1:10

#整合去批次
scRNA <- RunHarmony(scRNA, reduction = 'pca',group.by.vars = "orig.ident",reduction.save = "harmony",plot_converge = T)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pc.num,reduction.name = "umap")

library(ggsci)
p1 <- DimPlot(object = scRNA, reduction = "harmony", pt.size = .1, group.by = "orig.ident") + 
  scale_color_npg()+
  NoLegend()

p2 <- VlnPlot(object = scRNA, features = "harmony_1", group.by = "orig.ident", pt.size = .1) + 
  scale_fill_npg()+
  NoLegend()
p1 +p2 


obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = F)
obj <- ScaleData(obj, features = rownames(obj), assay = 'FLUX', verbose = F)
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 10, reduction.name = 'pca.flux', verbose = F)
obj <- FindNeighbors(obj, dims = 1:2, verbose = F)
obj <- FindClusters(obj, resolution = 0.5, verbose = F,graph.name = "RNA_snn")
obj <- RunTSNE(obj, dims = 1:2, assay = 'FLUX', reduction.name = "tsne.flux", verbose = F)
DimPlot(obj, reduction = "tsne.flux") + ggtitle('tSNE of Flux')




#reduction = "harmony" 
scRNA <- FindNeighbors(scRNA, dims = pc.num) 

# 聚类
scRNA <- FindClusters(scRNA)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)

#tSNE
scRNA <- RunTSNE(scRNA,dims = pc.num,dim.embed=2)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne.csv') 
plot1 = DimPlot(scRNA, reduction = "tsne",label = T) 
plot1
ggsave("cluster/tsne.pdf", plot = plot1, width = 8, height = 7)


#UMAP可视化
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap",label = T) 
plot2
ggsave("cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
library(patchwork)
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
plotc
ggsave("cluster/tSNE_UMAP.pdf", plot = plotc, width = 8, height = 4)

###3.细胞类型鉴定（SingleR）
### 安装SingleR   BiocManager::install('SingleR')
dir.create('cell_identify')
#手动注释
#habermann_imm <- c('CD274',"CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3",'GATA3', "KIT", "MKI67", "CDK1", "EPCAM")

#habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")

#DotPlot(scRNA, features = habermann_imm,group.by = "seurat_clusters") + coord_flip()
#DotPlot(scRNA, features = habermann_oth,group.by = "seurat_clusters") + coord_flip()

### 注释
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)#这里的n=10也即前10位的高变基因
#可视化marker基因
DoHeatmap(scRNA,features = top10$gene)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)
#celltype <- c('B_cells','T_cells','T_cells','B_cells',
#              'Fibroblasts','Endothelial_cells','Epithelial_cells','Epithelial_cells',
#              'Endothelial_cells','B_cells','Myeloids','Endothelial_cells',
#              'Mast_cells') 

library(SingleR)
#refdata <- celldex::HumanPrimaryCellAtlasData()
load("C:/Users/84759/Desktop/Meningioma/refdata_singleR.Rdata")
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata,
                    labels =refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

write.csv(celltype,"cell_identify/celltype_singleR.csv",row.names = F)
#celltype <- read.csv('cell_identify/celltype_singleR.csv')
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=4.5, reduction='umap')
p2
ggsave("cell_identify/UMAP_celltype.pdf", p2, width=8 ,height=7)
Idents(scRNA)=scRNA$celltype 

#三维
scRNA <- RunTSNE(scRNA,dims = pc.num,dim.embed=3)
data.combined <- scRNA
tmp.tsne.3<-Embeddings(object = data.combined[["umap"]])
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")
#笔者注，一般分的clusters较多时，常用的颜色配置函数配出来的颜色区分度不明显，笔者在网上搜到上述50个颜色分类，感觉比较好用，就作为自己常用的颜色条使用。
cb_palette.use <- cb_palette[1:length(unique(data.combined$celltype))]
col_match <- data.frame(cluster=unique(data.combined$celltype),col=cb_palette.use)
col_draw<- col_match[match(data.combined$celltype,col_match[,1]),2]
#library(rgl)
#plot3d(tmp.tsne.3,col = col_draw,type = 'p', radius = .001,axes=T,box=F)

library(plotly)
tmp.tsne.3 <- as.data.frame(tmp.tsne.3)
fig <- plot_ly(tmp.tsne.3, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color =data.combined$celltype, colors = cb_palette.use,size=2)
fig


genes <- c("COL7A1","ABCB9","TRIM10")
load('rf_lasso_genes.Rdata')
Idents(scRNA)=scRNA$orig.ident

scRNA1 <- subset(scRNA, TRIM10 > 0 |ABCB9 > 0|COL7A1 > 0)

VlnPlot(scRNA1,features = genes,group.by = 'celltype',pt.size = 0.2)
VlnPlot(scRNA1,features = genes,group.by = 'orig.ident',pt.size = 0.2)
FeaturePlot(scRNA1,features = genes,split.by = 'orig.ident')

markers <- FindMarkers(scRNA,features = genes,logfc.threshold = 0,min.pct = 0,min.cells.feature = 0,min.cells.group = 0 ,ident.1 = "BM", ident.2 = "CT") 
markers

saveRDS(scRNA,file ='scRNA_anno1.RDS')
#############################################
## 重要热图
table(scRNA@meta.data$celltype)
# 加入tissue_cell列
tissue_cell=paste0(scRNA@meta.data$orig.ident," ",scRNA@meta.data$celltype)
scRNA <- AddMetaData(scRNA,tissue_cell,'tissue_cell')
table(scRNA@meta.data$tissue_cell)

##生成一个参考矩阵，并不重要，后面都要替换里面的数字
ap <- AverageExpression(scRNA,features = genes,group.by ='tissue_cell' ,slot = 'data')[[1]]
## 循环
rowgroup=genes
colgroup=colnames(ap)
## 循环的逻辑是：先取每种细胞，再取每个基因，计算该细胞种类中阳性(表达>0)的细胞比例占所有细胞个数的百分比
## 下面的循环的巧妙之处在于提取自己sobj，然后在sobj中巧妙地提取S4对象中的RNA矩阵
## 稍微等几十秒
for( i in rowgroup){
  for(j in colgroup){
    sobj <- subset(scRNA,tissue_cell==j)
    ap[i,j] <- as.numeric(table(sobj@assays$RNA@data[i,]>0)[2])/length(sobj@assays$RNA@data[i,])*100
  }}
## NA是0，替换
ap[is.na(ap)] <- 0

## 热图
### ap:average proportion
ap=ap[,c('CT_Astrocytes','BM_Astrocytes','CT_Macrophages','BM_Macrophages','CT_Microglials','BM_Microglials','CT_Neural stem cells','BM_Neural stem cells','CT_Oligodendrocytes','BM_Oligodendrocytes')]
p=pheatmap::pheatmap(ap,display_numbers = T,  
                     color = colorRampPalette(c(rep("white",1), rep("#23BAC5",1)))(100),
                     cluster_rows = F,
                     cluster_cols = F,angle_col = 45,main = 'Proportion')
p
ggsave(filename = 'pheatmap_proportion.pdf',plot = p,width = 10,height = 5)


#############表达量热图##########################
#################################################
scRNA@meta.data$tissue_cell=paste0(scRNA@meta.data$orig.ident,' ',scRNA@meta.data$celltype)
library(Seurat)

ae=AverageExpression(scRNA,assays = 'RNA',group.by = 'tissue_cell',features = genes)
ae=as.data.frame(ae$RNA)
ae=log2(ae+1)
ae=na.omit(ae)
ae=ae[,c('CT_Chondrocytes','PR_Chondrocytes','CT_Endothelial_cells','PR_Endothelial_cells','CT_Macrophage','PR_Macrophage','CT_NK_cell','PR_NK_cell','CT_Monocyte','PR_Monocyte','CT_Tissue_stem_cells','PR_Tissue_stem_cells')]
p=pheatmap::pheatmap(ae,display_numbers = T,  
                     color = colorRampPalette(c(rep("white",1), rep("#FD763F",1)))(100),
                     cluster_rows = F,
                     cluster_cols = F,angle_col = 45 ,main='Expression')
p
ggsave(filename = 'pheatmap_expression.pdf',plot = p,width = 10,height = 5)

#画多个基因【共表达在一张图（scCustomize 实现）】
#scCustomize是一个单细胞转录组数据可视化的R包，里面集合了一些常用的数据可视化方法，可以与Seurat包进行联用。我们用Plot_Density_Joint_Only函数进行多基因联合密度图的绘制。
#install.packages("scCustomize")
#install.packages("Nebulosa")
library(Nebulosa)
library(scCustomize)
pbmc <- pbmc3k.final
p_density <- Plot_Density_Joint_Only(seurat_object = scRNA, 
                                     features = genes,
                                     custom_palette = BlueAndRed())
p_density +  DimPlot(scRNA,label = TRUE) 


p1=FeaturePlot(scRNA,genes,blend = T)

p2=FeaturePlot(scRNA, features = c("XBP1", "TTC28"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
               split.by = 'orig.ident',label = T) 
p2
ggsave(filename = 'XBP1-TTC28.pdf',plot = p2,width = 10,height = 5)


p3=FeaturePlot(scRNA, features = c("TRIM10"), cols =c("lightgrey", "orange"),pt.size = .9, order=T,
               split.by = 'orig.ident',label = T) 
p3
ggsave(filename = 'IGKC.pdf',plot = p3,width = 8,height = 4)

load("C:/Users/84759/Desktop/Meningioma/GS.Rdata")
scRNA <-AddModuleScore(scRNA, features= list,name = names(list))
names(x = scRNA[[]])

## 选择感兴趣的代谢通路
##提取细胞子集(DC)
Cells.sub <- subset(scRNA@meta.data, celltype=='B_cell')
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$orig.ident
DC_CT <- subset(scRNAsub, idents = c("CT"))
DC_PR <- subset(scRNAsub, idents = c("PSC"))

VlnPlot(scRNAsub, features="M5945_HALLMARK_HEME_METABOLISM41", group.by = "orig.ident", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Oxidative score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(DC_CT$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, DC_PR$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, alternative = "two.sided") #p-value < 2.2e-16

## 选择感兴趣的代谢通路
##提取细胞子集（表皮）
table(scRNA$celltype)
Cells.sub <- subset(scRNA@meta.data, celltype=='Monocyte')
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$orig.ident
Keratinocytes_CT <- subset(scRNAsub, idents = c("CT"))
Keratinocytes_PR <- subset(scRNAsub, idents = c("EC"))

VlnPlot(scRNAsub, features="M5945_HALLMARK_HEME_METABOLISM41", group.by = "orig.ident", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Oxidative score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(Keratinocytes_CT$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, Keratinocytes_PR$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, alternative = "two.sided") #p-value < 2.2e-16

scRNA$Oxidative_score=scRNA$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34
scRNA$Complement_score=scRNA$M5921_HALLMARK_COMPLEMENT23
p=FeaturePlot(scRNA, features = c("IGHG4", "Complement_score"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
              split.by = 'orig.ident',label = T) 

p
ggsave(filename = 'Oxidative_score.pdf',plot = p,width = 10,height = 5)


#########细胞通讯全局分析
# 先安装cellchat,4.1版本 的R为宜
#devtools::::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)
# 选取第一个样本，如果内存不够，再进一步减少细胞数，例如随机抽1000个
# 内存够，则不挑选直接上
#scRNA_chat <- subset(scRNA, orig.ident=='BM')
scRNA_chat <-subset(scRNA, tissue_type =='Tumor')
scRNA_chat <- scRNA
meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data

data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

cellchat <- createCellChat(object = data_input, meta = meta, group.by = "celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

write.csv(df.net,file ='cellchat.csv',quote=F)
#returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) 

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
dev.off()

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()

#####------------------------
# 关键热图

library(Seurat)
gene=read.table("C:/Users/84759/Desktop/预习单细胞/REACTOME_AGGREPHAGY.v2023.2.Hs.gmt")
gene=gene[,3:ncol(gene)]
gene=t(gene)
DoHeatmap(subset(scRNAsub,downsample=50,),features = genes,group.by = 'celltype',assay='RNA',slot = 'data',
          group.colors =c('#313c63','#b42e20','#ebc03e','#377b4c',
                          '#7bc7cd','#5d84a4'),lines.width = 10)+
  scale_fill_gradientn(colors=c('white','firebrick3'),na.value = 'white')
metadata=scRNA@meta.data


VlnPlot(scRNA,features = genes)
VlnPlot(scRNA,features = 'CDC34')
FeaturePlot(scRNA,features = genes,reduction = 'tsne',pt.size = 1)
