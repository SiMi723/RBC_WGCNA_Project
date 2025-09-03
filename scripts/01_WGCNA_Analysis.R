# ||======================================================================||
# ||                      SECTION 0: SETUP WORKSPACE                      ||
# ||======================================================================||

# --- 0.1: Clear workspace and set seed for reproducibility ---
rm(list = ls())
set.seed(12345)

# --- 0.2: Install and load necessary packages ---
# This code checks if packages are installed and installs them if not.
packages <- c("WGCNA", "DESeq2", "clusterProfiler", "org.Hs.eg.db", "pheatmap", "ggplot2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  # Install CRAN packages
  install.packages(packages[!installed_packages])
  # Install Bioconductor packages
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

# --- 0.3: Set up directories ---
# Note: RStudio projects automatically set the working directory.
# This structure assumes you are running from the root of the project.
data_dir <- "data"
results_dir <- "results"
plot_dir <- file.path(results_dir, "plots")
table_dir <- file.path(results_dir, "tables")


# ||======================================================================||
# ||                SECTION 1: DATA LOADING & PRE-PROCESSING              ||
# ||======================================================================||

# --- 1.1: Load gene expression and clinical trait data ---
# The expression data should have genes as rows and samples as columns.
# The first column should be the gene names.
counts_df <- read.csv(file.path(data_dir, "expression_data.csv"), row.names = 1)

# The trait data should have samples as rows and traits as columns.
# The first column should be the sample names.
traits_df <- read.csv(file.path(data_dir, "clinical_traits.csv"), row.names = 1)

# --- 1.2: Sanity Check - Ensure sample names match ---
# It's critical that the sample names in your count data columns
# match the sample names in your trait data rows.
if (!all(colnames(counts_df) == rownames(traits_df))) {
  stop("Sample names in expression data and trait data do not match!")
}

# --- 1.3: Filter low-count genes ---
# A gene must have at least 10 reads in total across all samples to be included.
keep <- rowSums(counts_df) >= 10
counts_df <- counts_df[keep, ]


# ||======================================================================||
# ||                    SECTION 2: DATA NORMALIZATION                     ||
# ||======================================================================||

# --- 2.1: Create DESeq2 object ---
# We use DESeq2 to normalize the data, which accounts for library size differences.
# The design formula is simple since we just need the transformation.
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = traits_df,
                              design = ~ 1) # No specific design needed for VST

# --- 2.2: Apply Variance Stabilizing Transformation (VST) ---
# VST is a robust normalization method for WGCNA as it stabilizes variance
# across the range of mean expression values.
vsd <- vst(dds, blind = TRUE)
normalized_counts <- assay(vsd)

# --- 2.3: Prepare data for WGCNA ---
# WGCNA requires genes to be in columns and samples in rows, so we transpose the matrix.
datExpr <- as.data.frame(t(normalized_counts))


# ||======================================================================||
# ||          SECTION 3: WGCNA - SOFT THRESHOLDING POWER SELECTION        ||
# ||======================================================================||

# --- 3.1: Choose a set of soft-thresholding powers ---
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# --- 3.2: Call the network topology analysis function ---
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# --- 3.3: Plot the results to select the best power ---
png(file.path(plot_dir, "01_soft_thresholding_power.png"), width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red") # This line corresponds to a standard R^2 cutoff

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# --- 3.4: Set the chosen power ---
# From the plot, we choose the power that is the lowest to reach the 0.85 line.
# For example, let's say we chose 6.
softPower <- 6


# ||======================================================================||
# ||             SECTION 4: WGCNA - NETWORK & MODULE CREATION             ||
# ||======================================================================||

# --- 4.1: Construct the gene network and identify modules in one step ---
# This function is the core of WGCNA. It will take some time to run.
net <- blockwiseModules(datExpr,
                        power = softPower,
                        TOMType = "unsigned",
                        minModuleSize = 30, # Minimum 30 genes per module
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25, # A higher value merges more modules
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

# --- 4.2: Plot the module dendrogram ---
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

png(file.path(plot_dir, "02_module_dendrogram.png"), width = 12, height = 9, units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# ||======================================================================||
# ||            SECTION 5: WGCNA - MODULE-TRAIT RELATIONSHIP              ||
# ||======================================================================||

# --- 5.1: Define the number of genes and samples ---
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# --- 5.2: Recalculate MEs with color labels ---
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# --- 5.3: Correlate MEs with the clinical trait ---
# We are interested in the 'RBC_Clogging' column of our traits data.
trait <- traits_df$RBC_Clogging
moduleTraitCor <- cor(MEs, trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# --- 5.4: Create the module-trait heatmap ---
# This is the most important plot of the analysis.
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png(file.path(plot_dir, "03_module_trait_heatmap.png"), width = 7, height = 10, units = "in", res = 300)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "RBC Clogging",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Module-Trait Correlation")
dev.off()


# ||======================================================================||
# ||                SECTION 6: ENRICHMENT OF KEY MODULE                   ||
# ||======================================================================||

# --- 6.1: Identify the most significant module ---
# From the heatmap, select the module with the highest absolute correlation.
# For this example, we assume it's the "turquoise" module.
target_module <- "turquoise"

# --- 6.2: Extract the genes from this module ---
module_genes <- colnames(datExpr)[moduleColors == target_module]

# --- 6.3: Perform Gene Ontology (GO) enrichment analysis ---
# We need to convert gene symbols to Entrez IDs for clusterProfiler.
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = module_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP", # Biological Process
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05)

# --- 6.4: Plot the enrichment results ---
p <- dotplot(go_enrich, showCategory = 15) + ggtitle("GO Enrichment for Turquoise Module")
ggsave(file.path(plot_dir, "04_GO_enrichment_dotplot.png"), plot = p, width = 10, height = 8)

# --- 6.5: Save the enrichment results and gene list to tables ---
write.csv(as.data.frame(go_enrich), file.path(table_dir, "turquoise_module_go_enrichment.csv"))
write.csv(data.frame(GeneSymbol = module_genes), file.path(table_dir, "turquoise_module_genes.csv"), row.names = FALSE)

print("Analysis complete. Check the 'results' folder for outputs.")
