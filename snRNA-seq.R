rm(list=ls())
library(Seurat)
setwd("E://scRNA/GSE157827/")
samples <- list.files("./AD4_NC7")

seurat_list <- list()
for (sample in samples) {
  data.path <- paste0("./AD4_NC7/", sample)
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}

scedata <- merge(seurat_list[[1]], 
                 y = seurat_list[-1],
                 add.cell.ids = samples)

scedata[["percent.mt"]] <- PercentageFeatureSet(scedata,pattern = "^MT-")
preQC_scedata <- VlnPlot(scedata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident",                         pt.size = 0)
preQC_scedata

scedata <- subset(scedata, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
postQC_scedata <- VlnPlot(scedata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident",                          pt.size = 0)
postQC_scedata

scedata <- NormalizeData(scedata)
scedata <- FindVariableFeatures(scedata, nfeatures = 4000)

top10 <- head(VariableFeatures(scedata), 10)
plot1 <- VariableFeaturePlot(scedata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(scedata)
scedata <- ScaleData(scedata, features = all.genes)
scedata <- RunPCA(scedata, features = VariableFeatures(object = scedata))

VizDimLoadings(scedata, dims = 1:2, reduction = "pca")
DimPlot(scedata, reduction = "pca") + NoLegend()
DimHeatmap(scedata, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scedata, dims = 1:15, cells = 500, balanced = TRUE)
library(harmony)
scedata <- RunHarmony(scedata, group.by.vars = "orig.ident")

b=DimPlot(scedata, reduction = "harmony", group.by = "orig.ident")
b

ElbowPlot(scedata, reduction = "harmony")
scedata <- FindNeighbors(scedata, reduction = "harmony", dims = 1:13)
scedata <- FindClusters(scedata, resolution = 0.5, reduction = "harmony")
scedata <- RunUMAP(scedata, reduction = "harmony", dims = 1:13)

DimPlot(scedata, reduction = "umap", label = T)
scedata.markers <- FindAllMarkers(scedata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(scedata.markers,file = 'scedata.markers.csv')

library("tidyverse")
sig.markers<- scedata.markers %>% select(gene,everything()) %>%
  subset(p_val_adj<0.05 & abs(scedata.markers$avg_log2FC)>1)
dim(sig.markers)

write.csv(sig.markers, file = ' sig.markers.csv')

new.cluster.ids <- c("0"="Oligodendrocytes", "1"="Astrocytes", "2"="Neuron", "3"="Neuron", "4"="Oligodendrocytes", "5"="Microglia",  "6"="Neuron", "7"="Oligodendrocyte progenitor", "8"="Endothelial", "9"="Neuron", "10"="Neuron", "11"="Neuron", "12"="Neuron", "13"="Endothelial", "14"="Neuron")
scedata <- RenameIdents(scedata, new.cluster.ids)
scedata$celltype <- scedata@active.ident
DimPlot(scedata, group.by = "celltype", label = T)

Neuron <- subset(scedata, celltype=="Neuron")
table(Neuron@meta.data$orig.ident)

diff_Neuron <- FindMarkers(Neuron, min.pct = 0.25, 
                           logfc.threshold = 0.25,
                           group.by = "orig.ident",
                           ident.1 ="AD",
                           ident.2="Normal")
write.csv(diff_Neuron, file = "diff_Neuron.csv")

log2FC = 0.4
padj = 0.05
diff_Neuron$threshold="ns"
diff_Neuron[which(diff_Neuron$avg_log2FC  > log2FC & diff_Neuron$p_val_adj < padj),]$threshold="up"
diff_Neuron[which(diff_Neuron$avg_log2FC  < (-log2FC) & diff_Neuron$p_val_adj < padj),]$threshold="down"
diff_Neuron$threshold=factor(diff_Neuron$threshold, levels=c('down','ns','up'))

p <- ggplot(data= diff_Neuron, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(diff_Neuron$threshold)[[1]],')'),'ns',
                                 paste0("up(",table(diff_Neuron$threshold)[[3]],')' )),
                     values=c("blue", "grey","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
p

selected_rows <- diff_Neuron[row.names(diff_Neuron) %in% "IDH3G", ]

p + geom_point(size = 3, shape = 1, data = selected_rows) +
  ggrepel::geom_label_repel(
    aes(label = label),
    data = selected_rows,
    color="black"
  )

VlnPlot(scedata, 
        features = c("IDH3G"),
        pt.size = 0,
        ncol = 1,
        split.by = "orig.ident")

singlecell_gene_test <- function(SerautObj, 
                                 genes.use, 
                                 group.by=NULL, 
                                 assay = "RNA", 
                                 comp = NULL, 
                                 alpha_start = .05, 
                                 Bonferroni = T,
                                 only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
  else{
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname]
      group1_exp <- group1_exp[which(group1_exp>0)] 
      
      
      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      group2_exp <- group2_exp[which(group2_exp>0)] 
      
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
    
  }
  
  if (Bonferroni == T){
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)
    
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
    
  }
  
  else{
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out)
    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))
  }
  
  return(dfOUT)
}

library(ggsignif)
A <- singlecell_gene_test(Neuron, 
                          genes.use = c('IDH3G'),
                          group.by = 'orig.ident', 
                          comp = c("AD", "Normal"))
A1 <- singlecell_gene_test(Neuron, 
                           genes.use = c('IDH3G'),
                           group.by = 'orig.ident', 
                           comp = c("AD", "Normal"),
                           only_postive = T)

A$p_val
A$sig

anno_pvalue <- format(A$p_val, scientific = T,digits = 3)
anno_sig <- A$sig
plots_violins <- VlnPlot(Neuron, 
                         cols = c("limegreen", "navy"),
                         pt.size = 0,
                         group.by = "orig.ident",
                         features = c('IDH3G'), 
                         ncol = 1, 
                         log = FALSE,
                         combine = FALSE)
plots_violins

for(i in 1:length(plots_violins)) {
  data <- plots_violins[[i]]$data
  colnames(data)[1] <- 'gene'
  plots_violins[[i]] <- plots_violins[[i]] + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y= element_text(size=12,color="black"),
          axis.title.x = element_blank(),
          legend.position='none')+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
    scale_x_discrete(labels = c("AD","Normal"))+
    geom_signif(annotations = anno_pvalue[i],
                y_position = max(data$gene)+0.5,
                xmin = 1,
                xmax = 2,
                tip_length = 0)
}
CombinePlots(plots_violins)
