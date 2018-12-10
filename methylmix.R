rm(list = ls());gc()
library(magrittr)
library(dplyr)
options(stringsAsFactors = F)
projectwd <- getwd() %>% print

# 注释
load('./result_RNAseq/diffAll_DESeq2.RData')
res_df <- data.frame(res_d)
DEGs <- res_df[abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, ]
# 
annot <- read.csv('./origin/esemble_symbol.csv')
annot[1:3,1:3]
index <- rownames(DEGs) %in% annot[,1];sum(index)
DEGs <- DEGs[index, ]
DEGs$symbol <- annot[match(rownames(DEGs), annot[,1]), 2]
DEGs$type <- annot[match(rownames(DEGs), annot[,1]), 3]
table(DEGs$type)

mRNA <- dplyr::filter(DEGs, type == 'protein_coding')
write.csv(mRNA, file = './result_mix/mRNA_diff.csv')
# 标准化表达数据，去重
norm_counts <- read.csv(file="./result_RNAseq/exprs_ESCC_DESeq2_normalized.csv",
                        check.names = F, row.names = 1)
norm_counts$symbol <- annot[match(rownames(norm_counts), annot[,1]), 2]
mRNA_eset <- dplyr::filter(norm_counts, symbol %in% unique(mRNA$symbol))
sum(!duplicated(mRNA_eset$symbol))

library(parallel)
eset <- data.frame(mclapply(mRNA_eset[,-93], function(x){
  tapply(x, mRNA_eset$symbol, mean)
}, mc.cores = 5), check.names = F)
write.csv(eset, file = './result_mix/DEGs_norm_exprs.csv')
cancer_eset <- eset[, -grep('11A$', names(eset))]
write.csv(cancer_eset, file = './result_mix/DEGs_cancer_norm_exprs.csv')


# exprs data
cancer_eset <- read.csv('./result_mix/DEGs_cancer_norm_exprs.csv', 
                        check.names = F, row.names = 1)
# methylation data
methydata <- data.frame(data, check.names = F)

intersect(rownames(cancer_eset), rownames(test)) %>% length()
intersect(names(methydata), names(cancer_eset)) %>%  length()

methy_cancer <- dplyr::select(methydata, one_of(names(cancer_eset)))
methy_normal <- dplyr::select(methydata, ends_with('11A'))

# intergrate data
library(MethylMix)
methy_cancer <- as.matrix(methy_cancer)
methy_normal <- as.matrix(methy_normal)
cancer_eset <- as.matrix(cancer_eset)

index <- intersect(rownames(cancer_eset), rownames(methy_cancer))
methy_cancer <- methy_cancer[rownames(methy_cancer) %in% index,]
methy_normal <- methy_normal[rownames(methy_normal) %in% index,]
cancer_eset  <- cancer_eset [rownames(cancer_eset ) %in% index,]

head(methy_cancer [, 1:4])
head(methy_normal[, 1:4])
head(cancer_eset[, 1:4])

all(colnames(methy_cancer ) == colnames(cancer_eset))
colnames(methy_normal)
all(rownames(methy_cancer) == rownames(cancer_eset))
all(rownames(methy_normal) == rownames(cancer_eset))

library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)
MethylMixResults <- MethylMix(methy_cancer, cancer_eset, methy_normal)
stopCluster(cl)
save(MethylMixResults, file = './result_interg/methylmix_result.RData')

MethylMixResults$MethylationDrivers
MethylMixResults$NrComponents
MethylMixResults$MixtureStates
MethylMixResults$MethylationStates[, 1:5]
MethylMixResults$Classifications[, 1:5]
MethylMixResults$Models[1]

plots <- MethylMix_PlotModel("ZNF69", MethylMixResults, methy_cancer)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("ZNF69", MethylMixResults, 
                             methy_cancer, METnormal = methy_normal)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("GLRB", MethylMixResults, 
                             methy_cancer, METnormal = methy_normal)
plots$MixtureModelPlot

plots <- MethylMix_PlotModel("GLRB", MethylMixResults, 
                             methy_cancer, cancer_eset, methy_normal)
plots$MixtureModelPlot
plots$CorrelationPlot

plots_a <- list()
for (gene in MethylMixResults$MethylationDrivers) {
  plots_a[[gene]] <- MethylMix_PlotModel(gene, MethylMixResults, 
                                         METcancer, GEcancer, METnormal)
}

drivers <- MethylMixResults$MethylationDrivers
write.csv(drivers, file = './result_interg/drvies.csv')

# ----------------
eset <- read.csv('./result_mix/DEGs_norm_exprs.csv', 
                 check.names = F, row.names = 1)
names(eset)
eset_dr <- eset[rownames(eset) %in% drivers,]
eset_dr_t <- data.frame(t(eset_dr), check.rows = F, check.names = F)
eset_dr_t$type <- sapply(rownames(eset_dr_t), function(x) {
  ifelse(grepl('11A$',x), 'Normal', 'Tumor')})
df <- eset_dr_t
library(ggplot2)
ggplot(df, aes(x=type, y=ABCB6, fill=type)) + geom_boxplot()
ggplot(df, aes(x=type, y=ABCB6, color=type)) + 
  geom_boxplot(outlier.colour = NA, width=.4) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill=NA)

diff_dr <- dplyr::filter(mRNA, symbol %in% drivers)
intersect(rownames(test), drivers) %>% length()

drivers


# -----------------------
# 2.Cox回归分析
# ------------------------
# 载入差异基因与生存数据
DEGs <- as.character(read.csv('./result_interg/drvies.csv', row.names = 1)[,1])
clic <- read.csv('./origin/clic/clinical.cart.2018-11-03/clic_sub4_status.csv')
# 生存数据整理
clic <- dplyr::arrange(clic, desc(submitter_id))
t1 <- clic$days_to_death == '--'
t2 <- !grepl('[0-9]',clic$days_to_last_follow_up)
t3 <- xor(t1,t2) 
clic[grep('F', t3),]

clic[t1,3] <- clic[t1,4]
clic$vital_status <- ifelse(clic$vital_status == 'dead',1,0)
clic$days_to_last_follow_up <- NULL
names(clic)

eset_dr_t[1:3,1:3]
clic[1:5,1:3]
expr[1:5,1:3]

# 合并
eset_dr_t$id <- rownames(eset_dr_t)
eset_dr_t$id <- substr(eset_dr_t$id,1,12)
duplicated(eset_dr_t$id) %>% sum
expr <- eset_dr_t[!duplicated(eset_dr_t$id),]
names(clic) <- c('id', 'fustat', 'futime')

baseline <- merge(expr, clic, by = 'id')
baseline[1:3,94:97]
baseline_type <- baseline$type
baseline$type <- NULL
baseline$futime <- as.numeric(baseline$futime)
str(baseline[,92:96])
baseline_log <- baseline
baseline_log[1:3,1:3]
names(baseline_log)[95]
for (i in 2:94) {
  baseline_log[i] <- log2(baseline_log[,i] + 1)
}
sum(sapply(baseline[,2:94], function(x) sum(x == 0)))
log2(baseline[,'WFIKKN2'])

# 单因素Cox回归分析
library(survival)
surv_data <- Surv(time = baseline_log$futime, event = baseline_log$fustat)
baseline$surv_data <- with(baseline_log, surv_data)

Unicox <- function(x) {
  FML <- as.formula(paste0('surv_data~', x))
  cox <- coxph(FML, data = baseline_log)
  goxsum <- summary(cox)
  HR <- round(goxsum$coefficients[,2], 2)
  Pvalue <- round(goxsum$coefficients[,5], 3)
  CI <- paste0(round(goxsum$conf.int[,3:4],2), collapse = '-')
  
  Unicox <- data.frame('characteristics' = x,
                       'Hazard Ratio' = HR,
                       'CI 95' = CI,
                       'P Value' = Pvalue)
  return(Unicox)
}

Univar <- lapply(drivers, Unicox)
Univar <- plyr::ldply(Univar, data.frame)

sig_univar <- Univar[Univar$P.Value < 0.05, ]
Univar$characteristics[Univar$P.Value < 0.05]

# 多因素Cox回归分析
fml2 <- as.formula(
  paste0('surv_data~', 
         paste0(Univar$characteristics[Univar$P.Value < 0.05], 
                collapse = '+')))
multicox <- coxph(fml2, data = baseline_log)
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




















