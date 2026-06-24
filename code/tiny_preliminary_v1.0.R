library(methylKit)
library(ggplot2)
library(ggrepel)

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>
# 
# # load workspace object if present and skip loading / analysis
# if (file.exists("code/methylKit_preliminary.RData")) {
#   load("code/methylKit_preliminary.RData")
# } else {
#   
# load bismark .cov file names from folder -------------------------------------

file.list <- as.list(list.files("data/EM_seq_files/tiny_cov", pattern = "\\.cov.gz$", full.names = TRUE)) # name list of all .cov files in folder

sample.names <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")

# methylKit analyses ------------------------------------------------------

# load files as methylRawList object
tiny_RawList=methRead(file.list,
               sample.id=sample.names,
               assembly="mm9", # just annotation
               treatment=c(0,0,0,0,0,0,0,1,1,1,1,1),
               context="CpG", # bismark .cov / bedgraph by default returns cpg context only
               mincov = 2, # default is 10x coverage
               pipeline = "bismarkCoverage"
)


# pull summary statistics into dataframe ----------------------------------
for(i in 1:length(tiny_RawList)){
  tiny_RawList_df_temp = data.frame(chr=tiny_RawList[[i]]$chr,
                               start=tiny_RawList[[i]]$start,
                               end=tiny_RawList[[i]]$end,
                               coverage=tiny_RawList[[i]]$coverage)
}


# filtered: discard cpg that have less than 10x coverage or 100% methylation
filtered.tiny_RawList=filterByCoverage(tiny_RawList,lo.count=3,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

# unite: retain only cpg that is covered across all samples
meth <- unite(filtered.tiny_RawList, min.per.group=1L)

# correlation matrix
corr <- as.matrix(read.table("data/EM_seq_files/corr.txt")) # already done in previous run... for some reason the methylkit function doesn't return the object so i saved it as a .txt file

# pca
pca <- PCASamples(meth, obj.return=TRUE)

nrow(meth)
dim(percMethylation(meth))

# save workspace objects for easy loading later ------------------------------------------

save.image(file = "code/methylKit_tiny_preliminary.RData")


# -------------------------------------------------------------------------

load("code/methylKit_tiny_preliminary.RData")

myDiff=calculateDiffMeth(meth)

print(myDiff[order(myDiff$qvalue),])

DE_df <- data.frame(myDiff, direction = c('negative','positive')[(myDiff$meth.diff>0)+1])
DE_df[abs(DE_df$qvalue) > 0.05,]$direction <- 'none'

ggplot(data=DE_df, aes(x=meth.diff, y=-log10(qvalue), col=direction))+
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 2) +
  scale_color_manual(values = c("cyan3", "grey", "red3"))

volcanoplotter <- function(DE, title, ylim = c(0, 80)) {
  
  # creating new df out of DE_renamed (gene IDs instead of ENSEMBL) for sig genes labeling
  DE_df <- as.data.frame(DE_renamed)
  volcano_labels <- DE_df[DE_df$padj < 0.05 & abs(DE_df$log2FoldChange) > 1,]
  volcano_labels <- subset(volcano_labels, log2FoldChange >= -10 & log2FoldChange <= 10) # so we stop labeling out-of-bounds dots
  
  # make volcano plot 
  DE_volcano <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = direction)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 2) + 
    geom_text_repel(data = volcano_labels, 
                    aes(label=rownames(volcano_labels)), size = 3, 
                    show.legend = FALSE) +
    scale_color_manual(values = c("cyan3", "grey", "red3")) +
    coord_cartesian(ylim = ylim, xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Direction', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
         title = title, subtitle = 'padj < 0.05 and absolute log2FoldChange > 1') +
    scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis
  
  return(list(
    DE_volcano = DE_volcano,
    tophits = volcano_labels[order(volcano_labels$padj),])
  )
  
}


# plots ----------------------------------------------------------------

# pca plot
# plotting df prep
plot_pca <- as.data.frame(pca$x)
plot_pca$group <- c('WT','WT','WT','WT','WT','WT','WT','KO','KO','KO','KO','KO')

ggplot(plot_pca, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(plot_pca)))
