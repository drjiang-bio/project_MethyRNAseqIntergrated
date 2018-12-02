#### RNAseq 及 Methylation 联合分析

**目录结构：**

![](https://github.com/drjiang-bio/project_MethyRNAseqIntergrated/blob/master/screenshots/ProjectContents-2018-12-03.png)

**工作流程：**

* 1. 批量读入RNAseq，差异分析（DESeq2, edgeR, limma）
  2. 批量读入Methylation, 整合成symbol--beta value，差异分析(wilcoxTest)
  3. 对DEGs和wilcox result 取交集