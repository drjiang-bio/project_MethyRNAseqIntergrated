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
#load(file = 'Rdata/step3.df_DMP.Rdata')
library(clusterProfiler)
library(org.Hs.eg.db)

length(unique(df_DMP$gene))
s2e <- bitr(unique(df_DMP$gene), fromType = "SYMBOL",
            toType = c( "ENTREZID"),
            OrgDb = org.Hs.eg.db)
df_DMP=merge(df_DMP,s2e, by.y='SYMBOL', by.x='gene')
table(!duplicated(df_DMP$ENTREZID))

gene_up= unique(df_DMP[df_DMP$change == 'UP','ENTREZID'])
gene_down=unique(df_DMP[df_DMP$change == 'DOWN','ENTREZID'])
gene_diff=c(gene_up,gene_down)
gene_all=unique(df_DMP$ENTREZID)


kkgo_file = "./result//kkgo_file.Rdata"
if(!file.exists(kkgo_file)){
  kk <- enrichKEGG(gene = gene_diff,
                   universe = gene_all,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)
  go <- enrichGO(gene_diff, OrgDb = "org.Hs.eg.db", ont="all") 
  save(kk,go,file = kkgo_file)
}
load(kkgo_file)

dotplot(kk)
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=45))
