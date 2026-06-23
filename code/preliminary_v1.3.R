# link R packages installed in temp directory on cluster
.libPaths(c("/data/horse/ws/shli842i-p_dna15_1/rpacks", .libPaths()))

library(methylKit)
library(ggplot2)
library(ggrepel)

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

# set wd to home folder on cluster (writeable)
setwd('/home/shli842i/p_dna15')

load("code/methRawList_load-methraw_v1.0.RData")
  
# filtered: discard cpg that have <=3x coverage or <=30x
filtered.methRawList=filterByCoverage(methRawList,lo.count=3,lo.perc=NULL,
                              hi.count=30, hi.perc=NULL)

# unite: merge all samples into one object
methBase <- unite(filtered.methRawList, destrand=FALSE, min.per.group = 1L)

# # correlation matrix
# corr <- as.matrix(read.table("data/EM_seq_files/corr.txt")) # already done in previous run... for some reason the methylkit function doesn't return the object so i saved it as a .txt file

# pca
pca <- PCASamples(methBase, obj.return=TRUE)

# coverage distribution data frame ----------------------------------------

methRawList.dfList <- lapply(methRawList, getData)

# omit_huge <- function(x){
#   print('anothuh one')
#   return(x$coverage[which(x$coverage<30)])
#   }
# 
# methRawList.dfList <- lapply(methRawList.dfList, omit_huge)
# 
# par(mfrow=c(3,4), mar=c(4,4,2,1), oma=c(0,0,3,0))
# 
# for(i in 1:length(methRawList.dfList)){
#   hist(methRawList.dfList[[i]], xlim=c(0,30), main=sample.names[i])
# }

# mtext("Methylation Coverage (counts)", side=3,line=1,outer=TRUE,cex=2,font=2)


# pca plot ----------------------------------------------------------------

# plotting df prep
pca_plotdf <- as.data.frame(pca$x)
pca_plotdf$group <- c('WT','WT','WT','WT','WT','WT','WT','KO','KO','KO','KO','KO')

# plot
pca_plot <- ggplot(pca_plotdf, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca_plotdf)))

ggsave(
  filename = "results/pca-plot_preliminary_v1.2.pdf",
  plot = pca_plot,
  width = 8,
  height = 6,
  units = "in"
)

save("methBase", file = "code/methBase_preliminary_v1.3.RData")