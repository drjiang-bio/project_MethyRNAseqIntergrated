rm(list = ls());gc()
library(magrittr)
library(dplyr)
options(stringsAsFactors = F)
projectwd <- getwd() %>% print
RNAseqd <- 'origin/rnaseq/'

# 文件准备：1.path : 待整合原始数据（建议使用绝对路径）
#           2.gdc_sheet ; gdc_sample_sheet文件
#           3.annot : 注释文件（检查数据的缺失值及重复值）
fpath <- paste0(projectwd, '/', RNAseqd, 'gdc_download_20181103/') %>% print

gdc_sheet <- read.csv('./origin/esca_gdc_sheet_pathology.csv')

# 检查文件：1. path 为 gdc下载后解压的各文件夹的上一级目录
#           2. gdc_sheet,annot 需和函数里的各参数一致（避免TCGA更新后各文件内容调整而出错）

# ------------------------------------------------
# 1.TCGA数据预处理（矩阵读入，注释，对应临床样本）
# ------------------------------------------------

# ------ 1.1 批量读入原始下载的数据，读入成表达矩阵

## 自定义函数1：批量读入gdc下载的数据（rnaseq_counts,fpkm）
load_tcga_gz <- function(lujing, guize = 'counts.gz$') {
  ## lujing为tcga下载的原始文件解压后的目录（其下一级目录为gz文件所在的文件夹）
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
exprs <- load_tcga_gz(lujing = fpath)
exprs <- data.frame(exprs, check.names = F)
exprs[1:3,1:3]
save(exprs, file = './result_RNAseq/exprs.RData')
write.csv(exprs, file = './result_RNAseq/exprs.csv')

# ------ 1.2 筛选鳞癌数据,构建分类表
sub_sheet <- filter(gdc_sheet, pathology == 'ESCC'| Sample.Type == 'Solid Tissue Normal' )
sub_sheet <- sub_sheet[, c(1,8,9)]
sub_sheet[1:3,1:3]
sub_sheet$Sample.Type <- ifelse(grepl('Normal', sub_sheet$Sample.Type), 
                                'Normal', 'ESCC')
# sub_sheet$Sample.ID <- gsub('-', '\\.', sub_sheet$Sample.ID)
escc <- select(exprs, one_of(sub_sheet$Sample.ID))

save(escc, file = './result_RNAseq/exprs_ESCC.RData')
write.csv(escc, file = './result_RNAseq/exprs_ESCC.csv')

group <- data.frame(id = names(escc))
group$class <- sub_sheet[match(group$id, sub_sheet$Sample.ID), 3]
group$class <- factor(group$class)

save(group, file = './result_RNAseq/exprs_ESCC_pheno.RData')
write.csv(group, file = './result_RNAseq/exprs_ESCC_pheno.csv')

# ----------------------------------------------
# 2.差异分析DESeq2
# ----------------------------------------------
rm(list = ls());gc()
options(stringsAsFactors = F)
library(DESeq2)
library(BiocParallel)
detectCores()
register(MulticoreParam(6))
library(amap)
library(gplots)
library(RColorBrewer)

load('./result_RNAseq/exprs_ESCC.RData')
load('./result_RNAseq/exprs_ESCC_pheno.RData')

nrow(escc)
escc <- escc[rowSums(escc)>5,]
head(escc)[1:3]

# 2.1 构建DESeq2分类表(colData)
group_list <- group$class
colData <- data.frame(row.names = group$id, group_list = group_list)
# 2.2 构建dds对象
dds <- DESeqDataSetFromMatrix(countData = escc,
                              colData = colData,
                              design = ~ group_list)
# 2.3 差异分析
dds <- DESeq(dds, parallel = T)
save(dds, file = './result_RNAseq/dds_ESCC.RData')

## 2.3.1  DESeq2 差异分析
res_d <- results(dds, contrast=c("group_list", "ESCC", "Normal"))
res_d <- res_d[order(res_d$pvalue),]
head(res_d)
summary(res_d)
dim(res_d)
# 结果输出
res_d[1:3,1:3]
save(res_d, file = './result_RNAseq/diffAll_DESeq2.RData')

## 2.3.2 edgeR 差异分析
library(edgeR)
load('./result_RNAseq/exprs_ESCC.RData')
# 设置分组信息，去除低表达量的gene以及做TMM标准化
nrow(escc)
escc <- escc[rowSums(cpm(escc) > 1) >= 2,]
names(escc)
group_list_e <- factor(c(rep("ESCC",81), rep("Normal",11)))
res_d_edgeR <- DGEList(counts = escc, group = group_list_e)
res_d_edgeR <- calcNormFactors(res_d_edgeR) 
res_d_edgeR <- estimateCommonDisp(res_d_edgeR)
res_d_edgeR <- estimateTagwiseDisp(res_d_edgeR)
# 结果输出
et <- exactTest(res_d_edgeR)
tTag <- topTags(et, n=nrow(res_d_edgeR))
tTag <- as.data.frame(tTag)
tTagSig <- tTag[tTag$PValue < 0.05 | tTag$FDR < 0.05,]
table(tTag$FDR < 0.05 & abs(tTag$logFC) > 1)
write.csv(tTagSig,file = "./result_RNAseq/diffSig_p_fdr_edgeR.csv")

## 2.3.3 返回标准化的数据(DESeq2)
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)[,1:3]
# 根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前。
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts_mad[1:4]
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
write.csv(normalized_counts, file="./result_RNAseq/exprs_ESCC_DESeq2_normalized.csv")

## 2.3.2 log转换后的结果
rld <- vst(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.csv(rlogMat, file="./result_RNAseq/exprs_ESCC_DESeq2_normalized_vst.csv")

# 热图绘制
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) # 生成颜色
pearson_cor <- as.matrix(cor(rlogMat, method="pearson")) # 计算相关性pearson correlation
hc <- hcluster(t(rlogMat), method="pearson") # 层级聚类
# 
pdf("./result_RNAseq/exprs_ESCC_DESeq2_normalized_vst.pdf", pointsize=10)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, 
          trace="none", col=hmcol, margins=c(11,11), 
          main="The pearson correlation of each sample")
dev.off()

pca_data <- plotPCA(rld, intgroup=c("group_list"), returnData=T, ntop=5000)

# -----------------------------
# 3 绘图
# -----------------------------
# 3.1 volcano
library(ggplot2)
library(DMwR)
vol_df <- data.frame(na.omit(res_d))
vol_df[complete.cases(vol_df),] %>% nrow 

vol_df$threshold <- as.factor(
  ifelse(vol_df$padj < 0.05 & abs(vol_df$log2FoldChange) >= 1,
         ifelse(vol_df$log2FoldChange > 1 ,'Up','Down'),'Not'))

pdf(file = './result_RNAseq/volcano.pdf', width = 5, height = 6.47)
ggplot(data = vol_df, aes(x = log2FoldChange, y = -log10(padj),
                              colour=threshold, fill=threshold)) +
  geom_point(alpha=0.4, size=1.2) +
  scale_color_manual(values=c("blue", "black","red"))+
  scale_x_continuous(expand = c(0,0)) +
  theme_bw(base_size = 12, base_family = "Times") +
  xlim(c(-9,9)) + ylim(c(0, 60)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.7)+
  geom_hline(yintercept =-log10(0.05) , lty=4,col="grey",lwd=0.7)+
  theme(legend.position=c(0.99,0.99), legend.justification = c(1,1),
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)", y="-log10 (p-value)",title = 'volcano')
dev.off()

# 3.2 heatmap
r_expr <- normalized_counts[rownames(normalized_counts) %in% 
                              c(diffSig_deseq2$ENSL[1:100], tail(diffSig_deseq2$ENSL,100)),]
r_expr <- na.omit(r_expr)
pdf(file="./result/heatmap.pdf",width=90,height=60)
par(oma=c(10,3,3,7))
hmMat=as.matrix(log10(r_expr+0.001))
heatmap.2(hmMat,col='greenred',trace="none")
dev.off()

library(pheatmap)
colnames(r_expr)
length(colnames(r_expr)) - sum(grepl('11A', colnames(r_expr)))
Type <- c(rep("tumor",81), rep("normal",11))    #修改正常和癌症样品数目
names(Type) <- colnames(r_expr)
Type <- as.data.frame(Type)

pheatmap(hmMat, annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, fontsize_row=5, fontsize_col=4)



heatmap <- data[rownames(data) %in% rownames(outDiff),]
write.table(heatmap,file="../result/Methy_heatmap.txt", 
            sep="\t", row.names=T, quote=F)

Type <- c(rep("normal",16),rep("tumor",96))    #修改正常和癌症样品数目
names(Type) <- colnames(heatmap)
Type <- as.data.frame(Type)

tiff(file="../result/Methy_heatmap.tiff", width = 45, height =70,
     units ="cm", compression="lzw", bg="white", res=300)
pheatmap(heatmap, annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, fontsize_row=5, fontsize_col=4)
dev.off()



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

# ---------------------------------------
# # 注释
# ---------------------------------------
rm(list = ls());gc()
diffSig_edgeR <- read.csv("./result_RNAseq/diffSig_p_fdr_edgeR.csv",
                          check.names = F, row.names = 1)
annot_EG <- read.csv('./origin/esemble_symbol.csv')
annot_EG[1:3,1:3]
index <- abs(diffSig_edgeR$logFC) > 1 & diffSig_edgeR$FDR < 0.05
difsig_fc <- diffSig_edgeR[index, ]

difsig_fc$symbol <- annot_EG[match(rownames(difsig_fc), annot_EG[,1]), 2]

duplicated(difsig_fc$symbol) %>% sum
difsig_fc[duplicated(difsig_fc$symbol),]
is.na(difsig_fc$symbol) %>% sum

difsig_fc <- na.omit(difsig_fc)
difsig_fc2 <- data.frame(mclapply(difsig_fc, function(x) {
  dt <- tapply(x, difsig_fc[,5], median, na.rm = T)
  }, mc.cores = 6), check.names = F)
difsig_fc2$status <- sapply(difsig_fc2$logFC, function(x) {
  ifelse(x > 0, 'UP', 'DOWN')
})
difsig_fc2 <- dplyr::arrange(difsig_fc2, logFC, status)

write.csv(difsig_fc2, file= "./result_interg/diffSig_p_fdr_edgeR_fc2_annot.csv")

###
dif_methy <- read.csv("./result_Methylation/diff_gene_methylation_fdr&fc2.csv",
                      row.names = 1, check.names = F)
dif_methy$symbol <- rownames(dif_methy)
inter_d <- intersect(rownames(dif_methy), difsig_fc2$symbol)
inter_me_rna <- merge(difsig_fc2, dif_methy, by='symbol')
write.csv(inter_me_rna, file = './result_interg/interg_methyRNAseq_diff.csv')


















