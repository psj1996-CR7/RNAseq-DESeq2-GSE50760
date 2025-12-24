# ============================================================
# RNA-Seq Differential Expression Analysis (GSE50760)
# ============================================================
# Author: PRASHANTH S JAVALI
# Description: Pipeline for processing raw RNA-seq count data,performing differential expression (DESeq2), and conducting pathway enrichment analysis (clusterProfiler).

# 0. INSTALL LIBRARIES 
install.packages("BiocManager")
cran_pkgs <- c("ggplot2","dplyr","pheatmap","ggrepel","readr","R.utils")
bioc_pkgs <- c( "DESeq2","clusterProfiler","org.Hs.eg.db","enrichplot", "DOSE")
install.packages(cran_pkgs)
BiocManager::install(bioc_pkgs)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(readr)
library(R.utils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tibble)

# 2. SETUP DIRECTORIES & DATA
dir.create("data", showWarnings = FALSE)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

# 3. Downloading GEO Dataset
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE50760&format=file"
geo_tar <- "data/GSE50760_RAW.tar"


# Untar the downloaded file
untar(geo_tar, exdir = "data/raw")

gz_files <- list.files("data/raw", pattern = "\\.gz$", full.names = TRUE)
if (length(gz_files) > 0) {
  sapply(gz_files, gunzip, overwrite = TRUE)
}

# 4. READ & MERGE COUNT FILES
raw_dir <- "data/raw"
files <- list.files(raw_dir, pattern = "txt$", full.names = TRUE)
if (length(files) == 0) stop("No count files found.")
read_sample <- function(f) {
  sample <- strsplit(basename(f), "_")[[1]][1]
  df <- read.table(
    f,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )

  df <- df[, 1:2]
  colnames(df) <- c("gene", sample)

  # Convert counts to numeric
  df[[sample]] <- as.numeric(df[[sample]])

  # Remove duplicated genes
  df <- df[!is.na(df[[sample]]), ]
  df <- df[!duplicated(df$gene), ]

  return(df)
}

# Merge all samples
counts <- read_sample(files[1])
for (f in files[-1]) {
  counts <- full_join(counts, read_sample(f), by = "gene")
}

# cleanup
rownames(counts) <- counts$gene
counts <- counts[, -1]
counts[is.na(counts)] <- 0
counts <- round(as.matrix(counts))

# Handeling Duplicates (Moved after counts is defined)
counts_df <- as.data.frame(counts) %>% 
  tibble::rownames_to_column(var = "gene")
counts_df_distinct <- counts_df %>% 
  distinct(gene, .keep_all = TRUE)
rownames(counts_df_distinct) <- counts_df_distinct$gene
counts <- counts_df_distinct %>% 
  dplyr::select(-gene) %>% 
  as.matrix()
counts[is.na(counts)] <- 0
counts <- round(counts)


# 5. LOAD METADATA
meta_file <- "metadata.csv"
coldata <- read.csv(meta_file, stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample
coldata$sample <- NULL
coldata$condition <- factor(coldata$condition,levels = c("normal", "CRC"))
stopifnot(all(colnames(counts) == rownames(coldata)))

# 6. DESEQ2 ANALYSIS
keep <- rowSums(counts >= 10) >= 3
dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = coldata,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "normal", "CRC"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
write.csv(
  res_df[order(res_df$padj), ],
  "results/GSE50760_DE_results.csv",
  row.names = FALSE
)

# 7. VISUALIZATION: PCA
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
ggsave("results/plots/PCA.png",
       ggplot(pca, aes(PC1, PC2, color = condition)) +
         geom_point(size = 4) + theme_classic() +
         ggtitle("PCA Plot"),
       width = 6, height = 5, dpi = 300)

## VOLCANO PLOT 
res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    status = case_when(
      log2FoldChange >  1 & padj < 0.05 ~ "Up",
      log2FoldChange < -1 & padj < 0.05 ~ "Down",
      TRUE                             ~ "NS"
    ),
    log2FC_plot = pmax(pmin(log2FoldChange, 6), -6)
  )

# Select top genes separately
top_up <- res_df %>%
  filter(status == "Up") %>%
  arrange(padj) %>%
  head(5)

top_down <- res_df %>%
  filter(status == "Down") %>%
  arrange(padj) %>%
  head(5)

label_genes <- bind_rows(top_up, top_down)

ggsave(
  "results/plots/Volcano.png",
  ggplot(res_df,
         aes(log2FC_plot, -log10(padj), color = status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 0, 1),
               linetype = c("dashed", "solid", "dashed"),
               color = "grey40") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed",
               color = "grey40") +
    geom_text_repel(
      data = label_genes,
      aes(label = gene),
      size = 4,
      max.overlaps = Inf
    ) +
    scale_x_continuous(limits = c(-6, 6)) +
    scale_color_manual(values = c(
      Down = "#E64B35",
      NS   = "grey70",
      Up   = "#4DBBD5"
    )) +
    theme_classic(base_size = 14) +
    labs(
      title = "Volcano Plot: CRC vs Normal",
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value"
    ),
  width = 7, height = 6, dpi = 300
)

top20 <- head(order(res$padj), 20)
mat <- assay(vsd)[top20, ]
mat <- mat - rowMeans(mat)

jpeg(
  filename = "results/plots/Heatmap.jpg",
  width = 10,      # inches
  height = 8,      # inches
  units = "in",
  res = 300,
  quality = 100
)

pheatmap(
  mat,
  annotation_col = coldata,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA
)

dev.off()

# TOP GENE BOXPLOT
top_gene_name <- res_df$gene[which.min(res_df$padj)]
df_gene <- plotCounts(dds, gene = top_gene_name, intgroup = "condition", returnData = TRUE)
ggsave("results/plots/TopGene_Boxplot.png",
       ggplot(df_gene, aes(condition, count, fill = condition)) +
         geom_boxplot() + geom_jitter(width = 0.2) +
         theme_classic() + ggtitle(paste("Expression of", top_gene_name)),
       width = 5, height = 6, dpi = 300)


# 8. PATHWAY ANALYSIS

# 8.1. Filter for significant genes
sig_genes_df <- res_df %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
# 8.2. Map Symbols to Entrez IDs
gene_map <- bitr(sig_genes_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
merged_df <- merge(sig_genes_df, gene_map, by.x = "gene", by.y = "SYMBOL")
gene_list <- merged_df$log2FoldChange
names(gene_list) <- merged_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)


# 9. GO ENRICHMENT ANALYSIS (BP, MF, CC)

ontologies <- c("BP", "MF", "CC")

for (ontology in ontologies) {
  # suppressMessages keeps the console clean from package-specific output
  ego <- suppressMessages(
    enrichGO(
      gene          = names(gene_list),
      OrgDb         = org.Hs.eg.db,
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      readable      = TRUE
    )
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego),
              paste0("results/GO_Enrichment_", ontology, "_Results.csv"),
              row.names = FALSE)

    p <- dotplot(ego, showCategory = 15) +
         ggtitle(paste("GO Enrichment:", ontology)) +
         theme_minimal() + # Theme adjustment for a cleaner look
         theme(plot.title = element_text(hjust = 0.5))

    ggsave(paste0("results/plots/GO_Dotplot_", ontology, ".png"),
           p, width = 8, height = 8, dpi = 300)
  }
}


# 10. KEGG PATHWAY ANALYSIS
kk <- enrichKEGG(gene = names(gene_list), organism = 'hsa', pvalueCutoff = 0.05)

if (!is.null(kk)) {
  kk_readable <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(as.data.frame(kk_readable), "results/KEGG_Results.csv")
  ggsave("results/plots/KEGG_Dotplot.png", dotplot(kk, showCategory = 15) + ggtitle("KEGG Pathways"), width = 8, height = 8, dpi = 300)
}


# 11. NETWORK PLOT (CNETPLOT)

if (!is.null(ego)) {
  p_net <- cnetplot(ego, categorySize = "pvalue", showCategory = 5, foldChange = gene_list, circular = TRUE, colorEdge = TRUE)
  ggsave("results/plots/Gene_Pathway_Network.png", p_net, width = 10, height = 8, dpi = 300)
}