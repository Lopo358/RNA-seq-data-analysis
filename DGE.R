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

