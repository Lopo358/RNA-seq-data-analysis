# https://ucdavis-bioinformatics-training.github.io/2018-September-Bioinformatics-Prerequisites/friday/limma_biomart_vignettes.html


# Combine all count data files across samples into one data.frame object
pkg <- c("edgeR", "DESeq2", "ggplot2", "tidyverse")
BiocManager::install(pkg)

## Check that you can load the packages, i.e., that the following runs without errors
suppressPackageStartupMessages({
  invisible(lapply(pkg, library, character.only = TRUE))
})

dir = "./gtf"
temp = list.files(path = dir, "*gene_count_matrix.csv", recursive = TRUE, full.names = TRUE)
temp_files = lapply(temp, read.csv)
head(temp_files[[1]])
lapply(temp_files, dim)
countData = reduce(temp_files, full_join)
head(countData)
colnames(countData)

colData <- read.csv(file.path(dir2, "PHENO_DATA.csv"), row.names=1) # PHENO_DATA lists the sample names and the groups they are in (i.e., control, treatment1, treatment2)
colData
group = factor(colData$Group, levels = c("control", "treatment1", "treatment2"))
dge$samples$group = group

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))


genes = do.call(rbind, str_split(geneid, pattern = "\\|"))
genes = as.data.frame(genes)
colnames(genes) = c("ENTREZID", "SYMBOL")

all(rownames(colData) %in% colnames(countData))
countData1 <- countData[, rownames(colData)]
countData1[is.na(countData1),] # some genes were only found in some of the samples
countData1 = countData1 %>% mutate(across(everything(), .fns = ~replace_na(.,0))) # replace NA values with 0 counts

all(rownames(colData) == colnames(countData))
#Create a DESeqDataSet from count matrix and labels
dds <- DESeqDataSetFromMatrix(countData = countData1,
                                colData = colData, design = ~ Group)
dds
#Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Group") #DESEQ2 plotPCA

# https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
dge <- DGEList(counts=countData1)

rownames(dge$samples) == rownames(colData)
dge$samples$group = colData$Group

# remove rows with zero or very low counts
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

# 15.4 Differential expression: limma-trend
logCPM <- cpm(dge, log=TRUE, prior.count=3)

contr.matrix <- makeContrasts(
  ControlvsTreat1 = control - treatment1,
  ControlvsTreat2 = control - treatment2,
  Treat1vsTreat2 = Treat1 - Treat2
  levels = colnames(design))

fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=1:3)

de <- decideTests(fit)
vennDiagram(de, circle.col=c("turquoise", "plum4", "lightgreen"), main = "Empirical Bayes stat for DE", cex = c(1,1,1))

png(file.path(out, "mean_variance_trend_model.png"), width = 1000)
par(mfrow=c(1,2))
v <- voom(dge, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()

# gene labels
geneid = rownames(fit)
length(geneid) # 
genes = do.call(rbind, str_split(geneid, pattern = "\\|"))
dim(genes)
genes = as.data.frame(genes)
colnames(genes) = c("ENTREZID", "SYMBOL")
class(genes)
table(genes$ENTREZID == genes$SYMBOL)
genes[genes$ENTREZID == genes$SYMBOL,]$SYMBOL = ""
genes

# volcano  plot 1
library(ggrepel)

clrs = rep("black", nrow(efit))
clrs[abs(fit$coefficients[,1]) > 2 & -log10(fit$p.value[,1]) > 2] = "red"

dat1 = as.data.frame(cbind(genes = genes$SYMBOL, log2FoldChange = fit$coefficients[,1], neglog10pval = -log10(fit$p.value[,1]), clrs = clrs))
str(dat1)
dat1$log2FoldChange = as.numeric(dat1$log2FoldChange)
dat1$neglog10pval = as.numeric(dat1$neglog10pval)
dat1$gene_labs = ""
dat1$gene_labs[abs(dat1$log2FoldChange) > 2 & dat1$neglog10pval > 2] = dat1$genes[abs(dat1$log2FoldChange) > 2 & dat1$neglog10pval > 2]
head(dat1)

ggplot(dat1, aes(x = log2FoldChange, y = neglog10pval, col = clrs, label = gene_labs)) +
  geom_point() +
  geom_text_repel(color = "gray30") +
  scale_color_manual(values=c("gray", "red")) +
  theme_bw() +
  geom_hline(yintercept = 2, col = "gray") +
  geom_vline(xintercept = c(-2,2), col = "gray") +
  theme(legend.position = "none") +
  ggtitle("Control vs Treatment1") +
  xlab("log2 Fold Change") +
  ylab("-log10(p-value)")
ggsave(file.path(out, "volcanoplot_eBayes_Control-vs-Treatment1.png"))

