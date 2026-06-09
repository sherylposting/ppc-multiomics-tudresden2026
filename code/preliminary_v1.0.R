library(methylKit)
library(ggplot2)
library(ggrepel)

# following this methylKit tutorial:
# https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html
# input bismark .cov format:
# <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count unmethylated>

# load workspace object if present and skip loading / analysis
if (file.exists("code/methylKit_preliminary.RData")) {
  load("code/methylKit_preliminary.RData")
} else {
  
  # load bismark .cov file names from folder -------------------------------------
  
  file.list <- as.list(list.files("data/EM_seq_files", pattern = "\\.cov.gz$", full.names = TRUE)) # name list of all .cov files in folder
  
  sample.names <- list("L188015_WT1","L188016_WT2","L188017_WT3","L188018_WT4","L188019_WT5","L188020_WT6","L188021_WT7","L188022_KO1","L188023_KO2","L188024_KO3","L188025_KO4","L188026_KO5")
  
  # methylKit analyses ------------------------------------------------------
  
  # load files as tabix... supposedly easier on the laptop's memory but not sure if it's actually faster...
  myobj=methRead(file.list,
                 sample.id=sample.names,
                 assembly="mm9", # just annotation
                 treatment=c(0,0,0,0,0,0,0,1,1,1,1,1),
                 context="CpG", # bismark .cov / bedgraph by default returns cpg context only
                 mincov = 2, # default is 10x coverage
                 pipeline = "bismarkCoverage",
                 dbtype = "tabix",
                 dbdir = "methylDB" # creates tabix database directory called methylDB
  )
  
  # filtered: discard cpg that have less than 10x coverage or 100% methylation
  filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
                                  hi.count=NULL,hi.perc=99.9)
  
  # unite: retain only cpg that is covered across all samples
  meth <- unite(myobj, destrand=FALSE)
  
  # correlation matrix
  corr <- as.matrix(read.table("data/EM_seq_files/corr.txt")) # already done in previous run... for some reason the methylkit function doesn't return the object so i saved it as a .txt file
  
  # pca
  pca <- PCASamples(meth, obj.return=TRUE)
  
  # save workspace objects for easy loading later ------------------------------------------
  
  save.image(file = "code/methylKit_preliminary.RData")

}

myDiff=calculateDiffMeth(meth)

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")

# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)


# plots ----------------------------------------------------------------

# pca plot
# plotting df prep
plot_pca <- as.data.frame(pca$x)
plot_pca$group <- c('WT','WT','WT','WT','WT','WT','WT','KO','KO','KO','KO','KO')

ggplot(plot_pca, aes(PC1, PC2, color=group)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(plot_pca)))
