rm(list = ls());gc()
WDIR <- '~/Desktop/PR2_Methylation/'
setwd(WDIR)
options(stringsAsFactors = F)
library(magrittr);library(dplyr);library(data.table);library(survival)
library(DMwR);library(impute)
library(ChAMP)

# 1. 载入phenotype、筛选出其中的鳞癌数据的tcga_id
# https://wemp.app/posts/8576bc3f-ec8c-4743-957e-81326ec74b04
# UCSC的xena浏览器下载指定ESCA的phenotype数据
pd.all <- fread('./TCGA-ESCA.GDC_phenotype.tsv.gz')
colnames(pd.all) %>% sort()
pd <- pd.all[,c("submitter_id.samples", "sample_type.samples", "disease_type")]
# 挑出squamous样本
pd$disease_type %>% table
pd <- pd[disease_type=="Squamous Cell Neoplasms"]

# 2.
# UCSC的xena浏览器下载指定ESCA的methylation beta数据
methyd <- fread('./TCGA-ESCA.methylation450.tsv.gz')
names(methyd)
methyd[,1:5]
methydf <- setDF(methyd[,2:ncol(methyd)], 
                 rownames = methyd[, `Composite Element REF`])
methydf[1:5,1:5]
# 去除有缺失值的位点（行）
#methydf <- methydf[complete.cases(methydf),]
methydf <- methydf[-manyNAs(methydf, 0.3),]
# 保留鳞癌数据
pd$submitter_id.samples %in% names(methydf) %>% sum
lie <- pd$submitter_id.samples[pd$submitter_id.samples %in% names(methyd)]
lie
methydf <- methydf[, lie]
### 整理为ChAMP对象
beta <- as.matrix(methydf)
# 处理beta信息，beta信号值矩阵里面不能有NA值
beta <- impute.knn(beta) 
sum(is.na(beta))
beta <- beta$data
beta[1:5,1:5]
beta <- beta + 0.00001
###处理pd信息
pd[1:5,]
names(pd) <- c('sampleID', 'class', 'patient')
pd$patient <- substr(pd$sampleID, 1, 12)
pd$class %>% table
pd$class <- ifelse(as.numeric(substr(pd$sampleID,14,15))<10,"Tumor","Normal")
pd <- setDF(pd)
pd <- pd[pd$sampleID %in% colnames(beta),]
all(pd$sampleID == colnames(beta))
rownames(pd) <- pd$sampleID

# 转化为ChAMP对象
myLoad <- champ.filter(beta = beta , pd = pd) #这一步已经自动完成了过滤
dim(myLoad$beta)
save(myLoad, file = './result/step1_myLoad.Rdata')
#save.image(file = 's1_squamous_data.RData')

#load('./result/step1_myLoad.Rdata')
myLoad$beta[1:4, 1:4]
myLoad$pd[1:4,]

norm_file = "./result/step2_champ_myNorm.Rdata"
if(!file.exists(norm_file)){
  myNorm <- champ.norm(beta=myLoad$beta, arraytype="450K", cores=8)
  save(myNorm, file = norm_file)
}
load(norm_file)


# 归一化过程产生了缺失值,需要将有NA的样本和它们的配对样本一起删掉
num.na <- apply(myNorm, 2, function(x)(sum(is.na(x))))
table(num.na)
##> num.na
##>      0 258616 260092 264579 
##>     61      1      1      1
#names(num.na) = colnames(myNorm)
#dt = names(num.na[num.na>0])
#dn = str_replace(dt,"-01","-11")
#keep = setdiff(colnames(myNorm),c(dt,dn))
#myNorm = myNorm[,keep]
#pd = myLoad$pd
#pd = pd[pd$sampleID %in% keep,]
identical(pd$sampleID, colnames(myNorm))


# 主成分分析
library(FactoMineR)
library(factoextra) 
dat <- t(myNorm)

group_list=pd$class
table(group_list)

dat.pca <- PCA(dat, graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE, 
             legend.title = "Groups")

# 热图
cg <- names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
ac <- data.frame(group=group_list)
rownames(ac) <- colnames(myNorm)  
pheatmap(myNorm[cg,],show_colnames =F,show_rownames = F,
         annotation_col=ac)
dev.off()

# 相关关系矩阵热图
pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F)
#剔除聚类失败的样本
#pn = c("TCGA-CV-5971-01","TCGA-CV-6953-11","TCGA-CV-6955-11")
#drop = str_sub(colnames(myNorm),1,12) %in% str_sub(pn,1,12)
#table(drop)
##> drop
##> FALSE  TRUE 
##>    52     6
#myNorm = myNorm[,!drop]
#dim(myNorm)
##> [1] 412481     52
#
#pd = pd[!(pd$patient %in% str_sub(pn,1,12)),]
#identical(pd$sampleID,colnames(myNorm))
##> [1] TRUE
#save(pd, myNorm, file = "./Rdata/step2_filtered_pd_myNorm.Rdata")


# 差异分析
group_list <- pd$class
myDMP <- champ.DMP(beta = myNorm, pheno=group_list)
head(myDMP$Tumor_to_Normal)

df_DMP <- myDMP$Tumor_to_Normal
df_DMP=df_DMP[df_DMP$gene!="",]
logFC_t <- 0.3
P.Value_t <- 10^-2
df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$logFC) > logFC_t,
                        ifelse(df_DMP$logFC > logFC_t ,'UP','DOWN'),'NOT') 
table(df_DMP$change) 
save(df_DMP, file = "./result//step3.df_DMP.Rdata")

library(dplyr)
library(ggplot2)
library(tibble)
dat  = rownames_to_column(df_DMP)
for_label <- dat %>% head(3)
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
# 标记某些个点
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = rowname),
    data = for_label,
    color="black"
  )
volcano_plot

# 热图
cg <-  rownames(df_DMP[df_DMP$change != "NOT",])
plot_matrix <- myNorm[cg,]
annotation_col <- data.frame(Sample=pd$class) 
rownames(annotation_col) <- colnames(plot_matrix)
ann_colors <- list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
pheatmap(plot_matrix,show_colnames = T,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors,show_rownames = F)

## 生存批量
library(survival)
library(survminer)
#load("./Rdata/step2_filtered_pd_myNorm.Rdata")
#load("./Rdata/step3.df_DMP.Rdata")

cg <- rownames(df_DMP[df_DMP$change != "NOT",])
myNorm_tumor <- myNorm[cg,1:100]
colnames(myNorm_tumor) <- substr(colnames(myNorm_tumor), 1, 15)

suv_dat <- data.table::fread("./TCGA-ESCA.survival.tsv.gz")
suv_dat$sample <- substr(suv_dat$sample,1,15)
suv_dat <- suv_dat[suv_dat$sample %in% colnames(myNorm_tumor),]
suv_dat <- suv_dat[substr(suv_dat$sample,14,15)=="01",] 
pd$sampleID <- substr(pd$sampleID, 1, 15)
suv_dat <- merge(pd, suv_dat, by.x = "sampleID", by.y = "sample")
myNorm_tumor <- myNorm_tumor[, suv_dat$sample]
identical(colnames(myNorm_tumor), suv_dat$sample)
#> [1] TRUE
library(survival)
logrankP <- apply(myNorm_tumor, 1, function(x){
  #x <- myNorm_tumor[1,]
  suv_dat$group <- ifelse(x>mean(x),"High","Low")
  res <- coxph(Surv(OS.time, OS)~group, data=suv_dat)
  beta <- coef(res)
  se <- sqrt(diag(vcov(res)))
  p.val <- 1 - pchisq((beta/se)^2, 1)
})
table(logrankP<0.05) #17个CpG位点

surv_gene <- names(sort(logrankP))[1:10] 
choose_matrix <- myNorm[surv_gene,]
annotation_col <- data.frame(Sample=pd$class) 
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

library(pheatmap)
pheatmap(choose_matrix,show_colnames = T,
         annotation_col = annotation_col,
         border_color=NA,
         color = colorRampPalette(colors = c("white","navy"))(50),
         annotation_colors = ann_colors)

gs=head(surv_gene,4)
exprSet = myNorm_tumor
meta = suv_dat
splots <- lapply(gs, function(g){
  meta$gene=ifelse(exprSet[g,]>median(exprSet[g,]),'high','low')
  sfit1=survfit(Surv(OS.time, OS)~gene, data=meta)
  ggsurvplot(sfit1,pval =TRUE, data = meta)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2)


#rm(list = ls())
#load(file = '/result/step3.df_DMP.Rdata')
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(stringr)

length(unique(df_DMP$gene))
s2e <- bitr(unique(df_DMP$gene), fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
df_DMP=merge(df_DMP,s2e, by.y='SYMBOL', by.x='gene')
table(!duplicated(df_DMP$))

gene_up= unique(df_DMP[df_DMP$change == 'UP','ENTREZID'])
gene_down=unique(df_DMP[df_DMP$change == 'DOWN','ENTREZID'])
gene_diff=c(gene_up,gene_down)
gene_all=unique(df_DMP$ENTREZID)

b <- bitr_kegg(gene_diff, fromType = "ENTREZID", 
               toType = 'UNIPROT', OrgDb = org.Hs.eg.db)
keytypes(org.Hs.eg.db)

kkgo_file = "./result//kkgo_file.Rdata"
if(!file.exists(kkgo_file)){
  kk <- enrichKEGG(gene = gene_diff,
                   universe = gene_all,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  enrichke
  go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all") 
  save(kk,go,file = kkgo_file)
}
load(kkgo_file)

dev.off()
dotplot(kk)

barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=45))









#load('./s1_squamous_data.RData')
# 载入生存信息
survd <- fread('./TCGA-ESCA.survival.tsv.gz')
methydf[1:5,1:5]
methydf <- as.data.frame(t(methydf), check.names=F)
methydf %>% rownames() %>% substr(1,12) %>% duplicated() %>% sum
methydf %>% rownames() %>% substr(14,17) %>% table

survdf <- setDF(survd)
survdf <- arrange(survdf, `_PATIENT`, sample, )
survdf <- survdf[!duplicated(survdf$`_PATIENT`),]
survdf[,1] %>% substr(14,17) %>% table

survdf <- survdf[survdf[,1] %in% rownames(methydf),]
methydf <- methydf[survdf[,1],]
all(rownames(methydf) == survdf[,1])
rownames(survdf) <- survdf[,1];survdf$sample <- NULL 

# 合并
baseline <- cbind(survdf, methydf)
baseline[1:4,1:5]

baseline$OS.time <- as.numeric(baseline$OS.time)
str(baseline[,1:5])

# 单因素Cox回归分析
library(survival)
surv_data <- Surv(time = baseline$OS.time, event = baseline$OS)
baseline$surv_data <- with(baseline, surv_data)

Unicox <- function(x) {
  #library(survival)
  #cox <- coxph(as.formula(paste0('surv_data~', x)), data = baseline)
  #goxsum <- summary(cox)
  goxsum <- coxph(as.formula(paste0('surv_data~', x)), 
                  data = baseline) %>% summary()
  HR <- round(goxsum$coefficients[,2], 2)
  Pvalue <- round(goxsum$coefficients[,5], 3)
  CI <- paste0(round(goxsum$conf.int[,3:4],2), collapse = '-')
  
  Unicox <- data.frame('characteristics' = x,
                       'Hazard Ratio' = HR,
                       'CI 95' = CI,
                       'P Value' = Pvalue)
  return(Unicox)
}


#load('./parData.RData')
Unicox1 <- function(x) {
  survival::coxph(as.formula(paste0('surv_data~', x)), 
        data = baseline) %>% summary()
}
Unicox2 <- function(x) {
  c('Hazard Ratio' = round(x$coefficients[,2], 2),
    'CI 95' = paste0(round(x$conf.int[,3:4],2), collapse = '-'),
    'P Value' = round(x$coefficients[,5], 3))
}

library(parallel)
cl.cores <- detectCores()
cl <- makeCluster(cl.cores)
clusterEvalQ(cl, c('survival', 'dplyr', 'Unicox1', 'Unicox2'))
clusterExport(cl, c("methydf", "baseline"))  #需要把这两个数据加载到并行计算环境里面
#system.time({ Univar <- parSapply(cl=cl, names(methydf)[1:10], Unicox, simplify=F) })
#system.time({ Univar <- parSapply(cl=cl, Univar, Unicox2, simplify=F) %>% do.call(rbind,.) })
#stopCluster(cl)
#save.image('./s2_squamous_data.RData')
#save(methydf, baseline, file = 'parData.RData')

stime <- Sys.time()
Univar <- parSapply(cl=cl, names(methydf), Unicox1, simplify=F)
Univar2 <- parSapply(cl=cl, Univar, Unicox2, simplify=F) %>% do.call(rbind,.)
stopCluster(cl)
#system.time({ Univar <- lapply(names(methydf)[1:30], Unicox) })
#Univar <- plyr::ldply(Univar, data.frame)
Univar2 <- data.frame(Univar2)
Univar2[,3] <- as.numeric(Univar2[,3])
Univar2$characteristics <- rownames(Univar2)
sig_univar <- Univar2[Univar2$P.Value < 0.001, ]
Univar2[Univar2$P.Value < 0.001,] %>% rownames() %>% length()
endtime <- Sys.time() - stime
save.image('./s3_squamous_data.RData')

#load('./s3_squamous_data.RData')
# 多因素Cox回归分析
options(expressions = 500000)
library(survival)
fml2 <- as.formula(
  paste0('surv_data~', 
         paste0(Univar2$characteristics[Univar2$P.Value < 0.001], 
                collapse = '+')))
multicox <- coxph(fml2, data = baseline)
multisum <- summary(multicox)

HR <- round(multisum$coefficients[,2], 2)
Pvalue <- round(multisum$coefficients[,5], 3)
CIL <- round(multisum$conf.int[,3],2)
CIU <- round(multisum$conf.int[,4],2)
CI <- paste0(CIL, '-' , CIU)

Multicox <- data.frame('characteristics' = rownames(multisum$coefficients),
                       'Hazard Ratio' = HR,
                       'CI 95' = CI,
                       'P Value' = Pvalue)
sig_multi <- Multicox[Multicox$P.Value < 0.05,]
# 合并Cox回归结果
cox_analysis <- merge(sig_univar, sig_multi, by = 'characteristics', all = T)
c <- c('Symbol', 'Hazard Ratio', 'CI', 'Pvalue','Hazard Ratio', 'CI', 'Pvalue')
cox_h <- rbind(c, cox_analysis)
write.table(cox_h, file = '../cox_analysis.csv', sep = ',',
            row.names = F, col.names = F)





multicox <- read.csv('Multi_sum2.csv')
statu <- read.csv('statu_time.csv')

sig_cox <- subset(multicox, Pr...z.. < 0.05)

exprs_statu <- merge(exprs, statu, by = 'gene_id')
exprs_statu <- exprs_statu[!duplicated(exprs_statu$gene_id), ]
row.names(exprs_statu) <- exprs_statu[, 1]

exprs_s <- exprs_statu[,2:6]
exprs_s <- exprs_s[, order(names(exprs_s))]
statu_i <- exprs_statu[,7:8]

names(exprs_s)
names(sig_cox)
exprs_s[1:3,1:5]

sig_cox <- sig_cox[order(sig_cox[, 1]), ]
sig_cox[,1] == names(exprs_s)


exprs_s$risk <- apply(
  sapply(1:5, function(i) { sig_cox$coef[i] * exprs_s[, i]}), 
  1, sum)

exprs_s$risktype <- 'High'
exprs_s$risktype[exprs_s$risk < median(exprs_s$risk)] <- 'Low'
statu_i <- cbind(exprs_s, statu_i)

library(survival)
Lung <- lung
surv_data <- Surv(time = statu_i$futime, event = statu_i$fustat)
surv_fit <- survfit(surv_data ~ statu_i$risktype)
sig_p <- survdiff(surv_data ~ statu_i$risktype)

plot(surv_fit, conf.int = 'none', col = c('red', 'blue'), lwd = 3, 
     mark.time = T, xlab = 'Time', ylab = 'Survival Probability')
legend(0.11, 0.2, c('High risk', 'Low risk'), 
       col = c('red', 'blue'), lwd = 2)
text(6000, 0.8, 'p = 7.17e-07', lwd = 3, col = 'red')

library(ROCR)
data(ROCR.simple)
pred <- prediction(statu_i$risk, statu_i$fustat) 
statu_i$fustat <- factor(statu_i$fustat)
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

summary(statu_i$risk)
