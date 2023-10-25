
#Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(hgu133plus2.db)
library(clusterProfiler)

# load series and platform data from GEO

gset <- getGEO("GSE108712", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL24460", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("100101100101110010101100001010X1011010000001101010",
               "00000110111101001110100X10110011101010101011110001",
               "10010000010100000011101010101101010100100010110110",
               "10111100101010100011101111100111101010101010101010",
               "11101010101010101010010001011011010101101101010110",
               "11011011101010010101001010100010111001001110010000",
               "11010101001010101001010100101011010101001011011001",
               "01001010101010101010010X10101101011011100110101000",
               "1001010110100X001010010101001001001010101101011000",
               "10100101100101010010100001011111000111000101001010",
               "10010100101010000001010100101010000001111010110110",
               "1010010111010110101101100")
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
groups <- make.names(c("tumour","normal"))
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
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("Name","adj.P.Val","P.Value","t","B","logFC"))


results <- tT[abs(tT$logFC)>0.75 & tT$adj.P.Val < 0.05, ]  

de_genes <- results[,1]
de_genes
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
