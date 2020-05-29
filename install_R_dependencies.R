required_pkgs <- c("plyr", "ShortRead", "ggplot2", "dada2", "optparse",
                   "Biostrings", "DECIPHER", "foreach", "doParallel",
                   "phylotools", "parallel")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

for (pkg in required_pkgs)
  if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg)
