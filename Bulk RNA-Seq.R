setwd("D://Rstudio_Data/GSE5281/code")
library(GEOquery)
gset = getGEO('GSE5281', destdir = '.', getGPL = F, AnnotGPL = T)
gset = gset[[1]]
expr = exprs(gset)
pdata = pData(gset)
View(pdata)

gset@annotation
probe = read.table(file = 'GPL570-55999.txt', sep = '\t', quote = '', comment.char = '#', header = T, fill = T, stringsAsFactors = F)
ids = probe[probe$Gene.Symbol != '', c(1,11)]
View(ids)

library(dplyr)
class(expr)
View(expr)
expr = as.data.frame(expr)
expr$ID = rownames(expr)
ids = ids[-grep('///',ids$Gene.Symbol),]
exprSet = inner_join(ids,expr,by = 'ID')

library(limma)
exprSet= avereps(exprSet[,-c(1,2)], ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
pdf (file = 'rowbox.pdf')
p <- boxplot(exprSet,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F)
title(main = list('Before normalization',cex = 2 ,font = 2), xlab = list('Sample list',cex = 1.5,font = 2), ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

library(limma)
normalized_expr = normalizeBetweenArrays(exprSet)
pdf(file = 'normalized_box.pdf')
p1 <- boxplot(normalized_expr,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F)
title(main = list('Normalization',cex = 2 ,font = 2),
xlab = list('Sample list',cex = 1.5,font = 2),
ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

normalized_expr<- log2(normalized_expr)
write.csv(normalized_expr, file = "normalized_expr.csv")
library(gplots)
library(edgeR)
library(dplyr)
library("edgeR")
data <- read.csv('normalized_expr.csv',header = T,row.names = 1)
marker <- 'IDH3G'
data1 <- data[which(rownames(data) %in% marker),]
data1 <- as.data.frame(t(data1))
View(data1)
data1$group <- ifelse(data1$IDH3G > median(data1[,'IDH3G']), 'high','low')

group_list = pdata$title
control = normalized_expr[,grep('control',group_list)]
AD = normalized_expr[,-grep('control',group_list)]
exprSet1 = cbind(AD, control)
group_list = c(rep('AD', ncol(AD)),
rep('normal',ncol(control)))
table(group_list)

exprSet1 = log2(exprSet1)
data = exprSet1
group_list = factor(group_list)
design <- model.matrix( ~0 + group_list)
colnames(design ) = levels(group_list)
rownames(design ) = colnames(data)
write.csv (design, file = "ClinicalTraits.csv")

contrast.matrix <- makeContrasts("AD-normal", levels = design)
fit <- lmFit( data, design )
fit2 <- contrasts.fit( fit, contrast.matrix )
fit2 <- eBayes( fit2 )
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="alldiff.xls",sep="\t",quote=F)

allLimma=allDiff
allLimma=allLimma[order(allLimma$logFC),]
allLimma=rbind(Gene=colnames(allLimma),allLimma)
write.table(allLimma,file="GSE5281_limmaTab.txt",sep="\t",quote=F,col.names=F)

logFoldChange=1
adjustP=0.05
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
write.csv(diffSig,file="diff.csv")
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
write.csv(diffUp,file="up.csv")
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)
write.csv(diffDown,file="down.csv")

hmExp=data[rownames(diffSig),]
write.csv(hmExp, file="hmExp.csv")
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))

library(ggplot2)
allDiff$change <- ifelse(allDiff$adj.P.Val < 0.05 & abs(allDiff$logFC) > 1,
ifelse(allDiff$logFC > 1,'UP','DOWN'),
'STABLE')
table(allDiff$change)
pdf(file = 'volcano.pdf')
ggplot(data= allDiff, aes(x = -log10(adj.P.Val), y = logFC, color = change)) +
geom_point(alpha=0.8, size = 1) +
theme_bw(base_size = 15) +
theme(plot.title=element_text(hjust=0.5),
panel.grid.minor = element_blank(),
panel.grid.major = element_blank()) +
geom_hline(yintercept= 0 ,linetype= 2 ) +
scale_color_manual(name = "",
values = c("red", "green", "black"),
limits = c("UP", "DOWN", "STABLE")) +
xlim(0,xMax) +
ylim(-yMax,yMax) +
labs(title = 'Volcano', x = '-Log10(adj.P.Val)', y = 'LogFC')
dev.off()

top100exp = exprSet1[rownames(diffSig)[1:100],]
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)
library(pheatmap)
pdf(file = 'heatmap.pdf')
pheatmap(top100exp,annotation_col = annotation_col,
color = colorRampPalette(c("green", "black", "red"))(50),
fontsize  = 5)
dev.off()

top50exp_Up = exprSet1[rownames(diffUp)[1:50],]
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)
library(pheatmap)
pdf(file = 'heatmap_top50_Up.pdf')
pheatmap(top50exp_Up,annotation_col = annotation_col,
color = colorRampPalette(c("green", "black", "red"))(50),
fontsize  = 5)
dev.off()

top25exp_Down = exprSet1[rownames(diffDown)[1:25],]
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)
library(pheatmap)
pdf(file = 'heatmap_top25_Down1.pdf')
pheatmap(top25exp_Down,annotation_col = annotation_col, cluster_rows = TRUE, cluster_cols = FALSE, color = colorRampPalette(c("green", "black", "red"))(50), fontsize  = 8)
dev.off()
