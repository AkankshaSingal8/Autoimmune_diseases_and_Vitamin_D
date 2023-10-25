# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(hgu133plus2.db)
library(clusterProfiler)

# load series and platform data from GEO

gset <- getGEO("GSE24287", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXX00000000000000000000000",
               "000000XXXXXXXXXXXXXXXXXX1111111111111111111111111")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Group1","Group2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="holm", sort.by="B", number=250)

tT <- subset(tT, select=c("Gene.symbol","adj.P.Val","P.Value","t","B","logFC"))

results <- tT[abs(tT$logFC)>0.5 & tT$adj.P.Val < 0.05, ]  

de_genes <- results[,1]
gene_ids <- AnnotationDbi::select(hgu133plus2.db, keys = de_genes, columns = "SYMBOL", keytype = "SYMBOL")

#perform the enrichment analysis
go_results <- enrichGO(gene = gene_ids$SYMBOL,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "SYMBOL",
                       ont = "BP",
                       pAdjustMethod = "holm",
                       pvalueCutoff = 0.05)
#View the enriched pathways
View(go_results@result)

