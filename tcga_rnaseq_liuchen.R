rm(list = ls());gc()
setwd('/media/jun/Sync/TCGA_data/ESCA/rnaseq')
options(stringsAsFactors = F)

# 文件准备：1.path : 待整合原始数据（建议使用绝对路径）
#           2.gdc_sheet ; gdc_sample_sheet文件
#           3.annot : 注释文件（检查数据的缺失值及重复值）
path <- '/media/jun/Sync/TCGA_data/ESCA/rnaseq/gdc_download_20181103/'
gdc_sheet <- read.csv('../esca_gdc_sheet_pathology.csv')
#
annot <- read.csv('../hgnc_complete_subset_complete.csv')
annot <- dplyr::arrange(annot, ensembl_gene_id)
# 
annot <- read.csv('../esemble_symbol.csv')
annot <- dplyr::arrange(annot, gene_id)

annot <- read.csv('../')

# 检查文件：1. path 为 gdc下载后解压的各文件夹的上一级目录
#           2. gdc_sheet,annot 需和函数里的各参数一致（避免TCGA更新后各文件内容调整而出错）

# ------------------------------------------------
# 1.TCGA数据预处理（矩阵读入，注释，对应临床样本）
# ------------------------------------------------

# ------ 1.1 批量读入原始下载的数据，读入成表达矩阵

## 自定义函数1：批量读入gdc下载的数据（rnaseq_counts,fpkm）
load_tcga_gz <- function(lujing, guize = 'counts.gz$') {
  ## path为tcga下载的原始文件解压后的目录（其下一级目录为gz文件所在的文件夹）
  ## guize 为数据文件的提取规则
  ## 例如：data/文件夹/gz
  
  # 根据正则表达式 获取FPKM文件的相对路径
  file_gz <- list.files(path = lujing, pattern = guize, recursive = T)
  # 获取其绝对路径
  ab_file_path <- paste(lujing, file_gz, sep = '/')
  
  # 批量读入匹配文件
  mdt <- lapply(ab_file_path, function(x) read.table(gzfile(x)))
  # 获取第一个list的第一列作为整合数据的行名
  row_name <- mdt[[1]][, 1]
  # 只保留所有的mdt元素的第二列
  mdt2 <- lapply(mdt, function(x) x[, 2])
  # 获取所有的文件夹名作为mydata所有元素的名称
  names(mdt2) <- sapply(strsplit(file_gz, '/'), function(x) x[[1]])
  
  # 合并数据
  df_mdt <- do.call(cbind, mdt2)
  # 简化ENSG并重命名行
  rownames(df_mdt) <- sapply(strsplit(row_name, '\\.'), function(x) x[[1]])
  # 用sample_id替换列名
  colnames(df_mdt) <- gdc_sheet[match(colnames(df_mdt), gdc_sheet$File.ID), 'Sample.ID']
  # 注释
  return(df_mdt)
}

## 使用自定义函数1 批量读入数据（默认读入count数据，其他数据需要改guize参数）
exprs <- load_tcga_gz(lujing = path)
exprs <- data.frame(exprs, check.names = F)
exprs[1:3,1:3]
save(exprs, file = './result/exprs.RData')

# ------ 1.2 筛选鳞癌数据,构建分类表
library(dplyr)
sub_sheet <- filter(gdc_sheet, pathology == 'ESCC'| Sample.Type == 'Solid Tissue Normal' )
sub_sheet <- sub_sheet[, c(1,8,9)]

sub_sheet[grepl('Normal', sub_sheet$Sample.Type), 3] <- 'Normal'
sub_sheet[!grepl('Normal', sub_sheet$Sample.Type), 3] <- 'ESCC'
# sub_sheet$Sample.ID <- gsub('-', '\\.', sub_sheet$Sample.ID)

escc <- select(exprs, one_of(sub_sheet$Sample.ID))
save(escc, file = './result/exprs_escc.RData')

group <- data.frame(id = names(escc))
group$class <- sub_sheet[match(group$id, sub_sheet$Sample.ID), 3]
group$class <- factor(group$class)
save(group, file = './result/exprs_escc_pheno.RData')

# ----------------------------------------------
# 2.差异分析DESeq2
# ----------------------------------------------
rm(list = ls());gc()
library(DMwR)
library(DESeq2)
library(RColorBrewer)
library(BiocParallel)
detectCores()
register(MulticoreParam(11))

load('./result/exprs_escc.RData')
load('./result/exprs_escc_pheno.RData')

nrow(escc)
escc <- escc[rowSums(escc)>5,]
head(escc)[1:3]
escc[, c(c('TCGA-V5-A7RC-06A','TCGA-LN-A5U6-01A'),'TCGA-XP-A8T6-01A')]

# 2.1 构建DESeq2分类表(colData)
group_list <- group$class
colData <- data.frame(row.names = group$id, group_list = group_list)
# 2.2 构建dds对象
dds <- DESeqDataSetFromMatrix(countData = escc,
                              colData = colData,
                              design = ~ group_list)
# 2.3 差异分析
dds <- DESeq(dds, parallel = T)
## 2.3.1 返回标准化的数据
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)[,1:3]
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts, file="./result/exprs_escc_DESeq2_normalized.csv",
          quote = F)
## 2.3.2 log转换后的结果
rld <- vst(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.csv(rlogMat, file="./result/exprs_escc_DESeq2_normalized_vst.csv",
            quote=F)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # 生成颜色
pearson_cor <- as.matrix(cor(rlogMat, method="pearson")) # 计算相关性pearson correlation
library(amap)
hc <- hcluster(t(rlogMat), method="pearson") # 层级聚类
# 热图绘制
pdf("./result/exprs_escc_DESeq2_normalized_vst.pdf", pointsize=10)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, 
          trace="none", col=hmcol, margins=c(11,11), 
          main="The pearson correlation of each sample")
dev.off()
pca_data <- plotPCA(rld, intgroup=c("conditions"), returnData=T, ntop=5000)

## 2.3.4  差异分析
res_d <- results(dds, contrast=c("group_list", "ESCC", "Normal"))
res_d <- res_d[order(res_d$pvalue),]
head(res_d)
summary(res_d)
dim(res_d)
res_d <- na.omit(res_d)
# 所有结果先进行输出
res_d[1:3,1:3]
write.csv(res_d, file="./result/diffAll_DESeq2.csv")

table(res_d$padj<0.05)
diffSig_deseq2 <-subset(res_d, padj < 0.05 & abs(log2FoldChange) > 1)
dim(diffSig_deseq2)
head(diffSig_deseq2)
summary(diffSig_deseq2)
sum(complete.cases(diffSig_deseq2))
write.csv(diffSig_deseq2, file= "./result/diffSig_DESeq2.csv")

# -----------------------------
# 3. 绘图.1
# -----------------------------
#volcano
pdf(file="./result/vol.pdf")
xMax=max(-log10(res_d$padj))+1
yMax=14
plot(-1*log10(res_d$padj),res_d$log2FoldChange, xlab="-log10(FDR)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=res_d[res_d$padj<0.05 & res_d$log2FoldChange>1,]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="red",cex=0.4)
diffSub=res_d[res_d$padj<0.05 & res_d$log2FoldChange<(-1),]
points(-log10(diffSub$padj), diffSub$log2FoldChange, pch=20, col="green",cex=0.4)
abline(h=0, lty=2, lwd=3)
dev.off()

# heatmap
res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))

red_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)
red_de_top20

r_expr <- normalized_counts[rownames(normalized_counts) %in% 
                              diffSig_deseq2_df$EnID[1:100],]
r_expr <- na.omit(r_expr)
library(gplots)
pdf(file="./result/heatmap.pdf",width=90,height=60)
par(oma=c(10,3,3,7))
hmExp=log10(r_expr+0.001)
hmMat=as.matrix(hmExp)
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()

# 注释
annot[1:3,1:3]
escc[1:3,1:3]
head(diffSig_deseq2)
diffSig_deseq2$symbol <- annot[match(rownames(diffSig_deseq2), annot[,1]), 2]
sum(is.na(diffSig_deseq2$symbol))
diffSig_deseq2 <- diffSig_deseq2[!is.na(diffSig_deseq2$symbol), ]

sum(duplicated(diffSig_deseq2$symbol))
diffSig_deseq2_df <- data.frame(diffSig_deseq2, row.names = diffSig_deseq2@rownames)
diffSig_deseq2_df$EnID <- rownames(diffSig_deseq2_df)
diffSig_deseq2_df <- dplyr::arrange(diffSig_deseq2_df, desc(abs(log2FoldChange)))
diffSig_deseq2_df <- diffSig_deseq2_df[!duplicated(diffSig_deseq2_df$symbol), ]
write.csv(diffSig_deseq2_df, file= "./result/diffSig_DESeq2_annot.csv", 
          row.names = F, quote = F)
length(intersect(annot[,1], row.names(outTab2)))

# ----------------------
# 4. 绘图.2
# ----------------------
oncogone <- read.csv('../ongene_human.csv')
TSG <- read.table('../Human_TSGs.txt', header = T, sep = '\t', fill = T)
diffSig_deseq2_df <- read.csv("./result/diffSig_DESeq2_annot.csv")
methy <- read.table("../methylation/result/Methy_GenediffSig.txt",
                    sep="\t", row.names=1, header = T)

outTab2 <- read.table('../methylation/result/Methy_Genediff.txt', 
                      sep = '\t', header = T, row.names = 1)
outTab2$betaDiff <- outTab2$tumorGeneMeans - outTab2$normalGeneMeans
write.csv(outTab2, file="../methylation/result/Methy_Genediff_betadiff.csv",
            row.names=T)

index <- abs(outTab2$logFC) > 1 & outTab2$pvalue < 0.05
outDiff=outTab2[index, ]
outDiff$type <- ifelse(outDiff$logFC > 0, 'UP', 'Down')
outDiffUP <- dplyr::filter(outDiff, type == 'UP')
outDiffDown <- dplyr::filter(outDiff, type == 'Down')
write.csv(outDiff, file="../methylation/result/Methy_GenediffSig_diff0.2.csv",
          row.names = T)

diffSig_deseq2_df <- dplyr::arrange(diffSig_deseq2_df, desc(log2FoldChange))
diffSig_deseq2_df$type <- ifelse(diffSig_deseq2_df$log2FoldChange > 0, 'UP', 'Down')
diffSigUP <- dplyr::filter(diffSig_deseq2_df, type == 'UP')
diffSigDown <- dplyr::filter(diffSig_deseq2_df, type == 'Down')

library(Vennerable)
venn <- Venn(list(diffSigUP$symbol, oncogone$OncogeneName, rownames(outDiffDown)))
venn2 <- Venn(list(diffSigDown$symbol, TSG$GeneSymbol, rownames(outDiffUP)))

intersect(diffSigUP$symbol, outDiffDown)
intersect(diffSigDown$symbol, TSG$GeneSymbol)

# ------------------------
# 二 甲基化数据处理
# ------------------------
rm(list = ls());gc()
setwd('/media/jun/Sync/TCGA_data/ESCA/methylation/download_test')
options(stringsAsFactors = F)

# 文件准备：
methysheet <- read.table('./gdc_sample_sheet.2018-11-14.tsv', 
                         header = T, sep = '\t', fill = T)
rna_sheet <- read.csv('/media/jun/Sync/TCGA_data/ESCA/esca_gdc_sheet_pathology.csv')
clinical <- read.table('../clinical.tsv', header = T, sep = '\t', fill = T)

# 添加病理类型
methysheet <- dplyr::arrange(methysheet, Sample.ID)
subRNAsheet <- rna_sheet[, c(1,8,9)]
subClicnical <- clinical[, c(1, 2, 3)]

methysheet$Case.ID %in% subClicnical$submitter_id
sum(duplicated(subClicnical$submitter_id))
methysheet$primaryType <- clinical[match(methysheet$Case.ID, 
                                         clinical$submitter_id), 11]
write.csv(methysheet, file = '../gdc_sample_sheet.2018-11-14_primarysite.csv', 
          row.names = F, quote = F)
# 提取鳞癌数据
index <- grepl('quamous', methysheet$primaryType) | grepl('Normal', methysheet$Sample.Type)
subMETHYsheet <- methysheet[index, ]
subMETHYsheet <- dplyr::arrange(subMETHYsheet, Sample.Type)
write.csv(subMETHYsheet, row.names = F, quote = F,
          file = '../gdc_sample_sheet.2018-11-14_primarysite_squamous.csv')


methyarry <- read.table('./geneMethy.txt', header = T, 
                        sep = '\t', row.names = 1, check.names = F) 
methyarry[1:5,1:3]
names(methyarry) <- substr(names(methyarry), 1, 16)
squasheet <- read.csv('../gdc_sample_sheet.2018-11-14_primarysite_squamous.csv')
table(squasheet$Sample.Type)
squasheet <- dplyr::arrange(squasheet, desc(Sample.Type))

submethy <- methyarry[, squasheet$Sample.ID]
substr(names(submethy),14,16)

# 差异分析
library(limma)
normalNum=16          #正常样品的数目
tumorNum=96           #癌症样品的数目
#读取数据
grade <- c(rep(1,normalNum),rep(2,tumorNum))
Type <- c(rep("Normal",normalNum),rep("Tumor",tumorNum))

rt <- submethy
rt[1:3,1:3]
rt <- as.matrix(rt)
data <- rt[rowMeans(rt)>0,]

#矫正数据
data <- normalizeBetweenArrays(data)
write.table(data,file="../result/normalizeMethy.txt",sep="\t",row.names=T,quote=F)

#差异分析
testl_ll <- apply(data, 1, function(x) {
  rt <- rbind(expression=x,grade=grade)
  rt <- as.matrix(t(rt))
  wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
  
  normalGeneMeans=mean(x[1:normalNum])
  tumorGeneMeans=mean(x[(normalNum+1):ncol(data)])
  logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)
  
  normalMed=median(x[1:normalNum])
  tumorMed=median(x[(normalNum+1):ncol(data)])
  diffMed=tumorMed-normalMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    return(c(normalGeneMeans, tumorGeneMeans, logFC, wilcoxTest$p.value))
  }
})
outTab2 <- data.frame(t(do.call(cbind, testl_ll)))
names(outTab2) <- c('normalGeneMeans', 'tumorGeneMeans', 'logFC', 'pvalue')

all(row.names(outTab) == row.names(outTab2))
all(round(outTab[,4],3) == round(outTab2[,4],3))
#对p值进行矫正
pValue=outTab2[,'pvalue']
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab2=cbind(outTab2, FDR=fdr)
#输出所有基因的甲基化差异情况
write.table(outTab2, file="../result/Methy_Genediff.txt", 
            sep="\t", row.names=T, quote=F)

#输出差异甲基化的基因
index <- abs(outTab2$logFC) > 1 & outTab2$FDR < 0.05
outDiff=outTab2[index, ]
write.table(outDiff, file="../result/Methy_GenediffSig.txt",
            sep="\t", row.names=T, quote=F)

#输出热图数据文件
heatmap <- data[rownames(data) %in% rownames(outDiff),]
write.table(heatmap,file="../result/Methy_heatmap.txt", 
            sep="\t", row.names=T, quote=F)
library(pheatmap)
Type <- c(rep("normal",16),rep("tumor",96))    #修改正常和癌症样品数目
names(Type) <- colnames(heatmap)
Type <- as.data.frame(Type)

tiff(file="../result/Methy_heatmap.tiff", width = 45, height =70,
     units ="cm", compression="lzw", bg="white", res=300)
pheatmap(heatmap, annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, fontsize_row=5, fontsize_col=4)
dev.off()

# 

<<<<<<< HEAD
=======























>>>>>>> e45fec05c38065d0e784263257e64dbee7d9d689
