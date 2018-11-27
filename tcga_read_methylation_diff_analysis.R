# Read TCGA methylation data(txt type)
# -----------------------------
#  一 甲基化数据读入与整合
# -----------------------------
rm(list = ls());gc()
setwd('/media/jun/Sync/TCGA_data/ESCA/methylation/download_test')
library(data.table)
library(parallel)
options(stringsAsFactors = F)

# data process
squasheet <- read.csv('../gdc_sample_sheet.2018-11-14_primarysite_squamous.csv')
file_gz <- list.files(path = './', pattern = 'hg38.txt$', recursive = T)
subfile <- sapply(strsplit(file_gz, '/'), function(x) x[[2]])
subfilegz <- file_gz[subfile %in% squasheet$File.Name]
# self function

read_methy <- function(file) {
  # file 为文件名，包括上一级目录名，最好上二级目录设置为过工作目录
  res <- c()
  ii = 0
  for (f in file) {
    n <- gregexpr("TCGA", f)[[1]][1]
    nam <- substr(f, n, n+15)
    re <- fread(f, sep = '\t', fill = T, header = T, 
                check.names = F)[, 2]
    res[nam] <- re
    ii = ii + 1
    print(ii)
  }
  re <- data.frame(fread(file[1], sep = '\t', fill = T, header = T)[, c(1,6)])
  res <- cbind(re, res)
  return(res)
}

methy_pre <- function(methy, method = 2) {
  # methy 为数据框，第一列为位点，第二列为symbol，余下为beta值
  # method = 2 : 取所有symbol; method = 1 : 取symbol列第一个
  # unique the same symbol-duplicated(still have multi symbol per row)
  if (method == 2) {
    symbol_uni <- sapply(strsplit(methy[,2],';'), unique)
    all_name <- unlist(symbol_uni)
    counts <- data.frame(table(all_name), check.names = F)
    # rename the symbol col by the unique symbol(still have multi symbol per row)
    methy[,2] <- sapply(symbol_uni, function(x) paste(x, collapse = ';'))
    
    # calculate the duplicate frep per different row
    rep_r_num <- sapply(methy[,2] , function(x) length(strsplit(x, ';')[[1]]))
    
    # duplicate row and rename the row
    rep_data <- methy[rep(1:nrow(methy), rep_r_num),]
    rep_data[,2] <- all_name
  } else {
    # abstact the fisrt symbol per row
    symbol_uni <- sapply(strsplit(methy[,2],';'), function(x) x[1])
    all_name <- unlist(symbol_uni)
    counts <- data.frame(table(all_name), check.names = F)
    # rename the symbol col by the unique symbol(still have multi symbol per row)
    methy[,2] <- sapply(symbol_uni, function(x) paste(x, collapse = ';'))
    counts <- data.frame(table(symbol_uni), check.names = F)
    rep_data <- methy
  }
  
  # calculate sum
  sum_data <- data.frame(mclapply(rep_data[, -c(1,2)], function(x) {
    dt <- tapply(x, rep_data[,2], sum, na.rm = T)
  }, mc.cores = 6), check.names = F)
  sum_data$counts <- counts[match(rownames(sum_data), counts[,1]),2]
  # calculate mean
  result <- data.frame(lapply(sum_data[, 1:(ncol(sum_data)-1)], 
                              function(x) x/sum_data$counts), check.names = F)
  rownames(result) <- rownames(sum_data)
  return(result)
}

sfmethy <- read_methy(subfilegz)
methyarraydf <- methy_pre(sfmethy)
methyarraydf[1:5, 1:3]

save(sfmethy, file = '../result3/methylation_read_orginial.RData')
write.csv(methyarraydf, file = '../result3/methyarray_pre.csv', row.names = T)

# ------------------------
# 二 甲基化数据处理
# ------------------------
# 差异分析
library(limma)
library(DMwR)
# 文件准备：
methyarraydf[1:5,1:3]
rt <- methyarraydf[,order(substr(names(methyarraydf),14, 16), decreasing = T)]
sum(grepl('11', substr(names(rt),14, 16)))
#
normalNum <- sum(grepl('11', substr(names(rt),14, 16)))          #正常样品的数目
tumorNum <- ncol(rt) - 16           #癌症样品的数目
grade <- c(rep(1,normalNum),rep(2,tumorNum))
Type <- c(rep("Normal",normalNum), rep("Tumor",tumorNum))

rt[1:3,1:3]
rt <- rt[-c(1, 2), ]
nrow(rt)
data <- rt[rowMeans(rt)>0,]
nrow(data)

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

#对p值进行矫正
fdr <- p.adjust(as.numeric(as.vector(outTab2[,'pvalue'])), method="fdr")
outTab2 <- cbind(outTab2, FDR=fdr)
#输出所有基因的甲基化差异情况
write.csv(outTab2, file="../result3/Methy_Genediff.txt", row.names=T)

#输出差异甲基化的基因
index <- abs(outTab2$logFC) > 1 & outTab2$FDR < 0.05
diffmethy_gene <- outTab2[index, ]
write.csv(diffmethy_gene, file="../result3/Methy_GenediffSig.csv", row.names=T)

#输出热图数据文件
heatmap <- data[rownames(data) %in% rownames(diffmethy_gene),]
write.csv(heatmap,file="../result3/Methy_heatmap.csv", row.names=T)
library(pheatmap)   
names(Type) <- colnames(heatmap)
Type <- as.data.frame(Type)

tiff(file="../result3/Methy_heatmap.tiff", width = 45, height =70,
     units ="cm", compression="lzw", bg="white", res=300)
pheatmap(heatmap, annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F, fontsize_row=5, fontsize_col=4)
dev.off()

# ---------------------
# 注释
# ---------------------
library(magrittr)
outTab2 <- read.csv(file="../result3/Methy_Genediff.txt", row.names=1)
diffmethy_gene <- read.csv(file="../result3/Methy_GenediffSig.csv", row.names=1)
methyWD <-getwd()
rnaseqWD <- '/media/jun/Sync/TCGA_data/ESCA/rnaseq'
setwd(rnaseqWD)

annotHGNC <- read.csv('../hgnc_complete_subset_complete.csv')
intersect(rownames(outTab2), annotHGNC$symbol) %>% length()
symbolOmit <- setdiff(annotHGNC$symbol, row.names(outTab2))
# 
annotENSM <- read.csv('../esemble_symbol.csv')
intersect(rownames(outTab2), annotENSM$gene_name) %>% length()

annotGPL <- fread('../methylation/GPL13534-11288.txt')
annotGPL[1:3,1:3]
annotGPL <- dplyr::arrange(annotGPL, ID)

methy <- fread('../methylation/download_test/02313093-32fa-4cda-a0ab-ff4432900a20/jhu-usc.edu_ESCA.HumanMethylation450.10.lvl-3.TCGA-VR-A8ET-01A-11D-A409-05.gdc_hg38.txt')
methy[1:3,1:3]
methy <- dplyr::arrange(methy, `Composite Element REF`)
all(annotGPL$ID == methy$`Composite Element REF`)
names(annotGPL)
annotGPL$ID %>% {head(., 20)}
annotGPL$UCSC_RefGene_Name %>% (function(x) head(x, 20))
annotGPL$UCSC_RefGene_Name %>% {head(., 20)}

methy$`Composite Element REF` %>% {head(., 20)}
methy$Gene_Symbol %>% {head(., 20)}

sapply(strsplit(methy$Gene_Symbol, ';'), unique) %>% {head(., 20)}
sapply(strsplit(annotGPL$UCSC_RefGene_Name, ';'), unique) %>% {head(., 20)}

