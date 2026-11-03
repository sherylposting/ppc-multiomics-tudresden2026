This project provides scripts for the analysis of bulk RNA-seq using DESeq2, ATAC-seq using DiffBind, and EM-seq using DMRichR.

# Overview
Bulk RNA-seq, ATAC-seq (chromatin accessibility), and EM-seq (DNA CpG methylation) data were used to infer differential activity of genes, transcription factors, and pathways, then integrated in order to identify consensus hits for follow-up investigation. Cdh1+ cells from the pancreata of day E14.5 mice were flow sorted in order to isolate pancreatic progenitors. Collection of bulk datasets was performed at the TU Dresden Genome Center, followed by quality control, trimming, read-mapping, and alignment with nf-core to the GRCm39 Mus musculus genome (June 2020). Differentially expressed genes (DEGs), differentially accessible chromatin peaks (DARs), and differentially methylated regions (DMRs) were respectively identified using DESeq2, DiffBind, and DMRichR. This project can be found under `bfx2557_RNA_seq_E14.5/`, `Sheryl/atac_seq_results/`, and `Sheryl/DMRichR_results/`.

Gene ontology (GO) and KEGG pathway enrichment testing were performed on the gene- or region-level differences using clusterProfiler and rGREAT. For RNA-seq, differential transcription factor (TF) activities between KO and WT samples were inferred from gene expression data using decoupleR / the DoRothEA database. For ATAC-seq and EM-seq data, HOMER was used to infer differentially accessible TF binding motifs. TF predictions were integrated using Robust Rank Aggregation (RRA). Within each analysis, TFs were ranked according to their adjusted p-value, with higher-ranked TFs representing stronger candidate regulators. The resulting list was ordered according to the nominal RRA scores. The code for this project can be found under `Sheryl/integration_results/R_scripts/`.

As a separate comparison of transcriptional programs associated with cell type-specific Aldh1b1 expression, publicly available scRNA-seq data from day E14 mouse pancreatic cells (GEO dataset GSE100622) was visualized and used for correlation analysis. The code for this project can be found under `Sheryl/integration_results/R_scripts/`. The raw data can be found as Seurat dge matrices under `Sheryl/data/GSE100619_RAW/`.

# Instructions for use
- The Github repo is only a copy of the code which is already present on the GAVALABS server, and is purely for the purposes of version control. All code must be run directly on GAVALABS, as you will otherwise not have access to the raw data.
- Code for the DMRichR analysis was written for use on the TU Dresden HPC, using R 4.2.1 and Bioconductor 3.16. Besides that, all other code was written for R 4.4.2 and Bioconductor 3.20 on Windows.
- Before running the code, ensure that:
  - you have installed all packages used in the analysis
  - you have edited all absolute paths which are specified at the top of all R scripts. (For example, all the absolute paths use Windows path formatting, and assume that you have mounted the GAVALABS server as V:)

### R Shiny App: PPC Multiomics Explorer
This repo features an R Shiny app, for convenient display of all the analyses (except for single-cell / Seurat). To use it, first install the Shiny package (`install.packages("shiny")`). Then open `Sheryl/integration_results/shiny_app/global.R` and press the "Run App" button which should appear at the top.

<img width="1917" height="978" alt="Screenshot 2026-11-03 203738" src="https://github.com/user-attachments/assets/021402ac-0fca-415b-ad21-8dd09d7c6cce" />

### RNA-seq analysis
- The RNA-seq analysis was previously performed on E14.5 and P10 pancreatic progenitors using the Dresden Genome Center's pipeline. However, my full analysis focuses only on E14.5.
- All code for the RNA-seq analysis is available in `bfx2557_RNA_seq_E14.5/`.
- For the purposes of this analysis, only the final table of differentially expressed genes from DESeq2 is used.

### ATAC-seq analysis
- The analysis was re-run on E14.5 pancreatic progenitors (`bfx2668_ATAC_E14.5/`) as well as P10 (`bfx2743_ATAC_P10`). It takes input from raw MACS2 narrow peak data which was previously processed.
- All code for the ATAC-seq analysis is available in `Sheryl/atac_seq_results/E14.5/R_scripts/`.
- The code is adapted from Adina's existing scripts (`Analysis_Adina/`). The main change is re-running the analysis using the mm39 reference genome, rather than mm10.

### EM-seq analysis
- The analysis was only run on day E14.5 pancreatic progenitors (`bfx2800_EMseq_E14.5/`). It was fully run on the TU Dresden HPC.
- Initially, the raw data was processed into cytosine reports files using Bismark. However, Raj performed a new analysis that attained higher coverage, which is available under `bfx2800_EMseq_E14.5/Raj_nfcore_analysis_Aug2026/CpG_reports`. It is also available on the HPC under `/projects/p_dna15/data/EM_seq_files/CpG_reports`.
- All code is available for this project under `/projects/p_dna15/`. However, as of the time when I am writing this, you will likely have to copy all the code over to your home folder in order to run it. I believe that SLURM jobs only have write permissions in the temp workspace and in your home directory, and not necessarily in the projects directory. Thus, any scripts which produce new results will need to work in your home directory.

#### A hacky installation of DMRichR-mm39
See the original DMRichR documentation: https://www.benlaufer.com/DMRichR/articles/DMRichR.html

The DMRichR package is no longer being maintained and is not compatible with the mm39 reference genome. As a result, you will need to implement some hacky solutions in order to run it with mm39...
1. Assuming you will conduct the analysis in your home directory, clone this repo to your home directory: `git clone https://github.com/sherylposting/ppc-multiomics-tudresden2026`. It will also redundantly clone some of the other scripts that I used for the full project. If desired, you can delete those other directories using `rm -r other-directories`.
2. I recommend setting up a new temporary workspace for yourself under `/data/horse/ws/` to hold the R package library and other temp files. The instructions for this can be found on the TU Dresden HPC compendium.
3. Additionally, the initial installation process for the regular DMRichR was a HUGE hassle due to outdated version incompatibilities with all its dependencies. I literally copied my library over to `/projects/p_dna15/rpacks_4.2.1_DMRichR/`, so hopefully for future purposes, you can simply copy this over to your workspace.
4. In your home directory on the HPC, clone the DMRichR-mm39 repo: `git clone https://github.com/sherylposting/DMRichR-mm39`. This is a fork of DMRichR which I edited for compatibility with mm39. However, it is still limited to running on R 4.2.1 and Bioconductor 3.16, leaving a lot of room for improvement to bring it fully up to date. For example, it is not yet able to run ontology enrichments with mm39.
5. Edit the `install-DMRichR.sbatch` script. You will have to edit LIBPATH to your new workspace (after copying rpacks_4.2.1.DMRichR/), and edit `export DMRICHR_PACKAGE_PATH` to the location of your cloned DMRichR-mm39 repo. Run this script to install it. It will then go on to run the `DMRichR/code/install-DMRichR.R` script. This is set up as a SLURM job due to the long time it will take to install all the dependencies, but theoretically you could also try doing this as an interactive session in RStudio, with `install-DMRichR.R` directly.
6. You will also have to clone the source package for the mm39 TxDb reference genome: `wget https://bioconductor.org/packages/3.17/data/annotation/src/contrib/TxDb.Mmusculus.UCSC.mm39.knownGene_3.17.0.tar.gz`. Then, go into the source and edit the DESCRIPTION file. Change the line which says: `Depends: GenomicFeatures` to `Depends: GenomicFeatures (>= 1.50.4)`. This fix is necessary so as not to disrupt the version of GenomicFeatures used by DMRichR-mm39.
7. Now, in your next R session, install it in the R console: `install.packages("your/path/to/TxDb.Mmusculus.UCSC.mm39.knownGene_3.17.0.tar.gz", repos=NULL, type="source")`.

#### Running DMRichR on HPC
If you have checked your logs files and the installation appears to have run successfully, you can then run `run_DMRichR.sbatch` to run the full DMRichR pipeline. First, to prepare, you will have to change the following global variables:
- In `run_DMRichR.sbatch`:
  - `LIBPATH` to your R packages library / "rpacks".
  - `VERSION` is an arbitrary version number which you can set to whatever you want. It will append to your working directory and create a new subfolder for the results. If you are running it for the first time you could name this something like "results" or "v1.0".
  - `WORKDIR` to your working directory. I used my home directory for this.
  - `SAMPLEINFO_PATH` to your sample_info.xlsx, which is required by DMRichR.
  - `DATADIR` to the input Bismark CPG reports. By default, you can keep this as `/projects/p_dna15/data/EM_seq_files/CpG_reports`.
  - `RESPATH` is generally used by DMRichR to seek any additional resources that may be needed for annotation purposes, such as `mm39.cpgislands.bed`, which was downloaded from the UCSC genome browser. It will run fine with the resources that already exist in `/projects/p_dna15/data/`.
- In `run_DMR.R`: you will also need to change the DMRichR-specific parameters, under the "DMRichR global parameters" section.

If it has successfully run, there will be a variety of results folders created by DMRichR in your working directory. DMRichR is prone to stopping mid-analysis due to errors in the code or OOM errors from SLURM. It saves some .RData files after each part of the analysis. Therefore, you may need to comment out chunks of the pipeline in `run_DMR.R`, uncomment the lines that load the .RData from the previous step, and re-run the job.

#### Running HOMER on HPC
The scripts assume that you have installed HOMER using conda in your home directory (following these instructions: http://homer.ucsd.edu/homer/introduction/install.html). You will then have to run `run_HOMER.sbatch`. This will first run `HOMER_v1.R` in order to prepare the hyper- and hypo-methylated BED files, then it will run the HOMER tool from command line.
