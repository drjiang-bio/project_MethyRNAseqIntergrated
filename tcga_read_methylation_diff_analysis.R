# Read TCGA methylation data(txt type)
# -----------------------------
#  一 甲基化数据读入与整合
# -----------------------------
rm(list = ls());gc()
library(data.table)
library(magrittr)
library(parallel)
library(limma)
library(DMwR)
options(stringsAsFactors = F)

projectwd <- getwd() %>% print
methyd <- 'origin/methylation/'

# data process
squasheet <- read.csv('./origin/methylation/gdc_sample_sheet.2018-11-14_primarysite_squamous.csv')
fpath <- paste0(projectwd, '/', methyd, 'download/') %>% print
squasheet$File.Name

read_methy <- function(lujing) {
  # ## lujing为tcga下载的原始文件解压后的目录（其下一级目录为txt文件所在的文件夹）
  # 此函数依赖函数外的sheet文件(squasheet)
  # 本函数返回值即为下一个函数的输入值, 
  # 返回一个数据框，第一列为位点，第二列为symbol，余下为beta值
  file_gz <- list.files(path = lujing, pattern = 'hg38.txt$', recursive = T)
  subfile <- sapply(strsplit(file_gz, '/'), function(x) x[[2]])
  # 依据sheet表中提取对应文件
  subfilegz <- file_gz[subfile %in% squasheet$File.Name]
  file <- paste0(lujing, subfilegz) # 获取其绝对路径
  
  res <- c()
  ii = 0
  for (f in file) {
    n <- gregexpr("TCGA-", f)[[1]][1]
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
  # 本函数返回值即为处理好的甲基化数据
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

sfmethy <- read_methy(lujing = fpath)
sfmethy[1:5,1:4]
methyarraydf <- methy_pre(sfmethy)
methyarraydf[1:5, 1:3]
methyarraydf1 <- methy_pre(sfmethy, method = 1)
methyarraydf1[1:5, 1:3]

save(sfmethy, file = './result_Methylation/methylation_read_orginial.RData')
write.csv(methyarraydf, file = './result_Methylation/methylation_pre_all.csv')
write.csv(methyarraydf1, file = './result_Methylation/methylation_pre_1.csv')

# ------------------------
# 二 甲基化数据处理
# ------------------------
# 差异分析
rm(list = ls());gc()
library(limma)

# 文件准备：
methyarraydf <- read.csv('./result_Methylation/methylation_pre_all.csv', 
                         check.names = F, row.names = 1)
methyarraydf <- methyarraydf[-1,]
methyarraydf[1:5,1:3]
# 事先排序
rt <- methyarraydf[,order(substr(names(methyarraydf),14, 16), decreasing = T)]
substr(names(rt),14, 16)
# 事先计算
sum(grepl('11', substr(names(rt),14, 16)))
normalNum <- sum(grepl('11', substr(names(rt),14, 16)))   # 正常样品的数目
tumorNum <- ncol(rt) - normalNum          #癌症样品的数目

Type <- c(rep("Normal",normalNum), rep("Tumor",tumorNum))
ID <- names(rt)
Type <- data.frame(ID = ID, class = Type)

rt[1:3,1:3]
nrow(rt)
data <- rt[rowMeans(rt)>0,]
nrow(data)
data[1:3,1:3]
2^log[1:3,1:3]

#矫正数据
data <- normalizeBetweenArrays(data)
log <- log2(data)
write.csv(data, file="./result_Methylation/methylation_pre_all_normalize.csv")

data <- read.csv("./result_Methylation/methylation_pre_all_normalize.csv", 
                 row.names = 1)

#差异分析1
methy_diff <- function(data_methy, normalNum, tumorNum) {
  # 本函数需对beta矩阵 事先排序, 将正常样本至于矩阵前列
  # 本函数需 事先计算 对照（正常）与实验组（癌症）的样本数目
  # dara_methy:处理好的基因总体平均甲基化beta值矩阵，列为样本，行为基因;
  # normalNum: 正常样品的数目; tumorNum: 癌症样品的数目
  # grade: 样本分类信息，使用数字向量(1,1,2,2,2)，对应矩阵样本，须正常样本在前列
  grade <- c(rep(1, normalNum), rep(2,tumorNum))
  res_l <- apply(data_methy, 1, function(x) {
    rt <- rbind(expression=x, grade=grade)
    rt <- as.matrix(t(rt))
    wilcoxTest <- wilcox.test(expression ~ grade, data=rt)
    
    normalGeneMeans = mean(x[1:normalNum])
    tumorGeneMeans = mean(x[(normalNum+1):ncol(data_methy)])
    logFC = log2(tumorGeneMeans)-log2(normalGeneMeans)
    
    normalMed = median(x[1:normalNum])
    tumorMed = median(x[(normalNum+1):ncol(data_methy)])
    diffMed = tumorMed - normalMed
    if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
      return(c(normalGeneMeans, tumorGeneMeans, logFC, wilcoxTest$p.value))
    }
  })
  res_df <- data.frame(t(do.call(cbind, res_l)))
  names(res_df) <- c('normalGeneMeans', 'tumorGeneMeans', 'logFC', 'pvalue')
  #对p值进行矫正
  fdr <- p.adjust(as.numeric(as.vector(res_df[,'pvalue'])), method="fdr")
  res_df <- cbind(res_df, FDR=fdr)
  return(res_df)
}

res_methy <- methy_diff(data, 16, 96)
#输出差异甲基化的基因
save(res_methy, file = "./result_Methylation/diff_gene_methylation_all.RData")
diff_res <- res_methy[res_methy$FDR < 0.05 & abs(res_methy$logFC) > 0.5,]
write.csv(res_methy, 
          file="./result_Methylation/diff_gene_methylation_fdr_fc05.csv")

# 差异分析2
limma_diff <- function(eset, type, FC = 2, p=0.05) {
  # eset, type 格式为 数据框（表达谱），注意检查数据格式
  # eset 列为样本，行为探针； type 为两列，第一列为样本id, 第二列为类型
  # 此函数只能用于分析两种类型
  
  # 1. 构建分类信息：1.核对分表表,2.提取实验分类信息
  ## 从表达矩阵中挑选出分类表中有的样本
  type <- type[match(colnames(log),type[, 1]),] 
  type <- factor(type[, 2])
  
  # 2. 设计实验矩阵design及对比模型contrast.matrix
  design <- model.matrix(~-1 + type)
  duibi <- paste(colnames(design), collapse = ' - ')
  contrast.matrix <- makeContrasts(contrasts = duibi,levels = design)
  
  # 3. 线性模型拟合，据对比模型行差值计算，贝叶斯检验
  fit <- lmFit(eset, design)                          # g
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit1)
  
  # 4. 检验结果报表及筛选
  dif <- topTable(fit2, coef = duibi, n = nrow(fit2), lfc = log2(FC))
  if (FC != 1) {
    dif <- dif[dif[, 'P.Value'] < p, ]
  }
  attr(dif,duibi)
  return(dif)
}
test <- limma_diff(log, Type, FC=2^0.5)
write.csv(test, 
          file="./result_Methylation/diff_gene_methylation_limma_fdr_fc05.csv")
intersect(rownames(test), rownames(diff_res)) %>% length

# ------------------------------------------------------------
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

# -----------------------------------------------
# 
# -----------------------------------------------
library(MethylMix)
data(METcancer)
data(METnormal)
data(GEcancer)
head(METcancer[, 1:4])
head(METnormal)
head(GEcancer[, 1:4])
all(colnames(METcancer) == colnames(GEcancer))
colnames(METnormal)
all(rownames(METcancer) == rownames(GEcancer))
all(rownames(METcancer) == rownames(METnormal))

library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
stopCluster(cl)

MethylMixResults$MethylationDrivers
MethylMixResults$NrComponents
MethylMixResults$MixtureStates
MethylMixResults$MethylationStates[, 1:5]
MethylMixResults$Classifications[, 1:5]

plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("MGMT", MethylMixResults, 
                             METcancer, METnormal = METnormal)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("ZNF217", MethylMixResults, 
                             METcancer, METnormal = METnormal)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("MGMT", MethylMixResults, 
                             METcancer, GEcancer, METnormal)
plots$MixtureModelPlot
plots$CorrelationPlot

plots_a <- list()
for (gene in MethylMixResults$MethylationDrivers) {
  plots_a[[gene]] <- MethylMix_PlotModel(gene, MethylMixResults, 
                                         METcancer, GEcancer, METnormal)
}

opar <- par(no.readonly=T)
par(mfrow=c(2,2))
par(opar)





