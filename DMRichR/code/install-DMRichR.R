# these should be specified in the sbatch script
LIBPATH <- Sys.getenv("LIBPATH", unset = "/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR")
DMRICHR_PACKAGE_PATH <- Sys.getenv("DMRICHR_PACKAGE_PATH", unset = "/home/shli842i/DMRichR")

CRAN_SNAPSHOT <- Sys.getenv("CRAN_SNAPSHOT", unset = "https://packagemanager.posit.co/cran/2023-04-25") # approx the last time DMRichR was updated

# link R packages installed in temp directory on cluster
.libPaths(LIBPATH)

cat("Loaded namespaces at startup:\n")
print(loadedNamespaces())

cat("Library paths:\n")
print(.libPaths())

dir.create(LIBPATH, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIBPATH, .libPaths()))

# Set CRAN before the first package installation
options(repos = c(CRAN = CRAN_SNAPSHOT),
       timeout = 3600)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
        "BiocManager",
        lib = LIBPATH,
#        repos = CRAN_SNAPSHOT,
        type = "source"
    )
}

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages(
        "remotes",
        lib = LIBPATH,
#        repos = CRAN_SNAPSHOT,
        type = "source"
    )
}

# Combine Bioconductor 3.16 with the dated CRAN snapshot
repos <- BiocManager::repositories(version = "3.16")
repos["CRAN"] <- CRAN_SNAPSHOT
options(repos = repos)

cat("Library paths:\n")
print(.libPaths())

cat("Repositories:\n")
print(getOption("repos"))

install.packages(
    c("ggplot2", "BSgenome.Mmusculus.UCSC.mm39", "org.Mm.eg.db"),
    lib = LIBPATH,
    repos = repos,
    type = "source"
)

remotes::install_local(
    DMRICHR_PACKAGE_PATH,
    lib = LIBPATH,
    dependencies = TRUE,
    upgrade = "never",
    force = TRUE,
    repos = repos,
    type = "source"
)