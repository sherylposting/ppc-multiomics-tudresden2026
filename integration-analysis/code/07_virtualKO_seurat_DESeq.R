library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(future)
library(DESeq2)

# global variables - check these ------------------------------------------

DATA_PATH <- "C:/Users/shery/Documents/bioinfo_windows/ppc_multiomics_tudresden2026/data/GSE100619_RAW/GSM2689399_e14_matrix.dge.txt.gz"

OUTPUT_DIR <- "V:/MARIA/Sheryl/integration_results/virtual_KO"

# load seurat counts data -------------------------------------------------

# fread to read in counts matrix without eating RAM
dt <- fread(DATA_PATH)

genes <- dt[[1]]
dt[[1]] <- NULL

counts <- as(as.matrix(dt), "dgCMatrix")
rownames(counts) <- genes

# Initialize the Seurat object with the raw (non-normalized data).
seu <- CreateSeuratObject(counts = counts, project = "E14")

# 21753 features across 15243 samples within 1 assay 
# Active assay: RNA (21753 features, 0 variable features)

# preprocessing ---------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "preprocessing.pdf"), width = 10, height = 8)

# mitochondrial dna percentage, counts outliers
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

hist(seu[["percent.mt"]]$percent.mt)
quantile(seu[["percent.mt"]]$percent.mt, probs = seq(0, 1, 0.01)) # 95% is 5.79 so i guess cut off any cells with >6% mito DNA
quantile(seu[["nFeature_RNA"]]$nFeature_RNA, probs = seq(0, 1, 0.01))[c("5%", "95%")] # 560.0, 3371.9
quantile(seu[["nCount_RNA"]]$nCount_RNA, probs = seq(0, 1, 0.01))[c("5%", "95%")] # 733.1, 8365.9

seu <- subset(
  seu,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 3400 &
    percent.mt < 6
)

# normalize
seu <- NormalizeData(seu)

# find variable features
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seu), 10)
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scale data
all.genes <- rownames(seu)
seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA", "percent.mt"))

# PCA
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
VizDimLoadings(seu, dims = 1:2, reduction = "pca")
DimPlot(seu, reduction = "pca") + NoLegend()
ElbowPlot(seu, ndims=100)
# dude there's two "elbows" i guess the first 40 PCs are relevant??

dev.off()

# UMAP --------------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "virtual_KO_seurat.pdf"), width = 10, height = 8)

plan(multisession, workers = 4)

seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu, resolution = 0.5)

seu <- RunUMAP(seu, dims = 1:40)
DimPlot(seu, reduction = "umap")

# cluster markers ---------------------------------------------------------

# find our own identifying positive markers of each UMAP cluster
# markers <- FindAllMarkers(
#   seu,
#   test.use = "bimod",
#   only.pos = TRUE
# )
# 
# markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1)

# use cluster markers from paper (for tSNE)
markers_paper <- list(
  Alpha  = c("Gcg", "Arx", "Pou3f4"),
  Beta   = c("Ins2", "Ins1", "Nkx6-1"),
  BP     = c("Spp1", "Cdh1", "Nkx6-1"),
  Endo   = c("Pecam1"),
  EP     = c("Neurog3"),
  Eryth  = c("Hbb-bt", "Hemgn"),
  Hepato = c("Afp", "Alb"),
  MΦ     = c("Mrc1", "Lyz2"),
  Mast   = c("Cpa3", "Cd52", "Nkg7"),
  Mes    = c("Itm2a", "Col3a1", "Vim", "Nkx2-5", "Col1a2", "Ptn"),
  MesO   = c("Upk3b", "Lrrn4", "Aldh1a2"),
  Neuron = c("Phox2b", "Ascl1", "Hand2", "Gap43", "Phox2a", "Tlx2", "Tbx3"),
  PaSC   = c("Nes", "Vim", "Vegfa", "Ncam1", "Pdgfra", "Pdgfrb"),
  TC     = c("Rbpjl", "Ptf1a", "Nr5a2", "Cpa1")
)

# find which markers from the paper aren't in the dataset + need to be renamed
lapply(markers_paper, function(x) {
  x[!x %in% rownames(seu)]
  })

grep("ngn", rownames(seu), value = TRUE, ignore.case = TRUE)

# fixed: nkx6-1, Neurog3, nkx2-5. k they're all good now

seu <- AddModuleScore(
  seu,
  features = markers_paper,
  name = "CellTypeScore"
)

score.cols <- paste0(
  "CellTypeScore",
  seq_along(markers_paper)
)

colnames(seu@meta.data)[
  match(score.cols, colnames(seu@meta.data))
] <- paste0("score_", names(markers_paper))

scores <- seu@meta.data[
  ,
  paste0("score_", names(markers_paper))
]

seu$marker_celltype <- names(markers_paper)[
  max.col(scores, ties.method = "first")
]

DimPlot(
  seu,
  reduction = "umap",
  group.by = "marker_celltype"
)

# tSNE --------------------------------------------------------------------

seu <- RunTSNE(
  seu,
  dims = 1:40
)

DimPlot(
  seu,
  reduction = "tsne",
  group.by = "marker_celltype"
)

# ok so the UMAP looks like cluster markers should be different but tSNE (same as paper) looks mostly ok

# subset aldh1b1 cells ----------------------------------------------------

# plot all ALDHs on tSNE
FeaturePlot(
  seu,
  features = all.genes[grepl("Aldh", all.genes)],
  reduction = "tsne"
)

FeaturePlot(
  seu,
  features = "Aldh1b1",
  reduction = "tsne"
)

# subset only cells with detectable aldh1b1
expr <- FetchData(seu, vars = "Aldh1b1")

# percentage of cells in each cluster (defined by paper) expressing aldh1b1
data.frame(
  cluster = seu@meta.data$marker_celltype,
  express_pct = expr$Aldh1b1 > 0
) |>
  aggregate(express_pct ~ cluster, data = _, FUN = mean) |>
  transform(express_pct = express_pct * 100)

# cluster express_pct
# 1    Alpha    6.567797
# 2     Beta    7.006369
# 3       BP   40.461216
# 4     Endo    2.702703
# 5       EP   28.499157
# 6    Eryth    8.785942
# 7   Hepato   21.052632
# 8     Mast    4.790419
# 9      Mes    4.025045
# 10    MesO    4.278922
# 11      MΦ    2.234637
# 12  Neuron    1.785714
# 13    PaSC    4.453441
# 14      TC   63.724305

VlnPlot(
  seu,
  features = "Aldh1b1",
  group.by = "marker_celltype"
)

# plot on tSNE our dummy dataset, e-cad sorting and double positive aldh1b1 cells
df <- FetchData(
  seu,
  vars = c("Aldh1b1", "Cdh1")
)

# fetch seurat whichCells for dummy subset to replicate our dataset - e-cad / cdh1 sorting
cdh1_cells <- subset(seu, 
       cells = WhichCells(
         seu,
         expression = Cdh1 > 0
  )
)

# double positive aldh1b1 cells
subset_cells <- subset(seu, 
                     cells = WhichCells(
                       seu,
                       expression = (Cdh1 > 0 & Aldh1b1 > 0 & Ptf1a > 0 & Cpa1 > 0)
                       # subset cdh1, aldh1b1, and TC identity cells
                     )
)

subset_cells_data <- FetchData(subset_cells, vars = "Aldh1b1")

DimPlot(
  seu,
  reduction = "tsne",
  cells.highlight = list(
    "Cdh1" = rownames(cdh1_cells@meta.data),
    "Cells" = rownames(subset_cells@meta.data)
  ),
  cols.highlight = c("brown1", "dodgerblue"),
  cols = "grey80"
)

# subset top and bottom 25% quartiles of aldh1b1+cdh1 cells
idx_25 <- quantile(subset_cells_data$Aldh1b1, probs = seq(0, 1, 0.25))["25%"]

bottom_25 <- WhichCells(cdh1_cells,
           expression = Aldh1b1 < idx_25)

idx_75 <- quantile(subset_cells_data$Aldh1b1, probs = seq(0, 1, 0.25))["75%"]

top_25 <- WhichCells(cdh1_cells,
                        expression = Aldh1b1 > idx_75)

# 0%       25%       50%       75%      100% 
# 0.6980327 1.3392364 1.7396888 2.1866297 3.7143529 

top_25 <- subset(
  seu,
  cells = top_25
) # 276 cells

top_25_counts <- GetAssayData(
  top_25,
  assay = "RNA",
  layer = "counts"
)

bottom_25 <- subset(
  seu,
  cells = bottom_25
)

bottom_25_counts <- GetAssayData(
  bottom_25,
  assay = "RNA",
  layer = "counts"
) # 1246 cells

dev.off()

# DEseq -----------------------------------------------------------------

pdf(file.path(OUTPUT_DIR, "virtual_DESeq.pdf"), width = 10, height = 8)

# pseudobulk according to run1, run2, run3 for each group
sample <- unique(sub("^[^.]+.", "top.", colnames(top_25_counts)))

top_25_agg <- sapply(sample, function(s) {
  Matrix::rowSums(
    top_25_counts[, sample == s, drop = FALSE]
  )
})

sample <- unique(sub("^[^.]+.", "bot.", colnames(top_25_counts)))

bottom_25_agg <- sapply(sample, function(s) {
  Matrix::rowSums(
    bottom_25_counts[, sample == s, drop = FALSE]
  )
})

counts_mat <- cbind(top_25_agg, bottom_25_agg)

# create coldata
coldata <- data.frame("condition" = c(
  rep("top", 3),
  rep("bottom", 3)
  )
)
rownames(coldata) <- colnames(counts_mat)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

res <- results(dds,
               contrast = c("condition", "bottom", "top")
)

res <- res[order(res$padj), ]

res <- as.data.frame(res)

# compare how well the two deseq results match ----------------------------

# load our real data
real_diffgenes <- read.delim("V:/MARIA/Sheryl/individual_hits_tables/RNA_diffgenes.tsv")

x <- res[which(res$padj < 0.2),]
y <- real_diffgenes[which(real_diffgenes$padj < 0.2),]
rownames(real_diffgenes) <- real_diffgenes$Gene_Symbol

common <- intersect(rownames(x), rownames(y))

x <- x[common, ]
y <- y[common, ]

cor(
  x$stat,
  y$stat,
  method = "spearman",
  use = "complete.obs"
)
# lol -0.0768155
# lol after subsetting to TC only, -0.09715787
# ok 0.2984615 after subsetting padj < 0.05
# wait there's only 26 shared significant genes...
# 0.3492393 if subset padj < 0.20. i don't know...

save(res, file=file.path(OUTPUT_DIR, "virtual_KO_df.RData"))

# volcano plot ------------------------------------------------------------

# prepare df for plotting
res$Gene_Symbol <- rownames(res)
res$direction <- "None"
res$direction[which(res$padj < 0.05 & res$log2FoldChange < 0)] <- "Negative"
res$direction[which(res$padj < 0.05 & res$log2FoldChange > 0)] <- "Positive"

# helper function to plot
volcanoplotter <- function(DF, TITLE, XLIM = NULL, YLIM = NULL, ...) {
  p <- ggplot(data=DF, aes(x=log2FoldChange, y=-log10(padj), col=direction))+
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point() +
    geom_text_repel(
      data = subset(DF, DF$padj < 0.03), # preferentially label higher padj
      aes(label = Gene_Symbol),
      show.legend = FALSE,
      ... # extra options here just if for some reason the labels look crazy on the plot
    ) +
    labs(
      y = "-log10(padj)",
      x = "log2 fold change",
      title = TITLE,
    ) +
    scale_color_manual(values=c(
      "None" = "gray", 
      "Positive" = "brown1", 
      "Negative" = "dodgerblue"))
  
  if (!is.null(XLIM))
    p <- p + xlim(XLIM)
  
  if (!is.null(YLIM))
    p <- p + ylim(YLIM)
  
  print(p)
}

volcanoplotter(res, "Virtual KO: top vs bottom quartile of Aldh1b1 expression in E14.5 tip cells", XLIM = c(-7.5, 7.5), YLIM = c(0,200))

volcanoplotter(res, "Virtual KO: top vs bottom quartile of Aldh1b1 expression in E14.5 tip cells, zoomed", XLIM = c(-7.5, 7.5), YLIM = c(0,100))

dev.off()
