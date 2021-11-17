# Set required packages.
req_pack <- c("iCluster",
              "NMF",
              "SNFtool",
              "OmicsPLS",
              "mixOmics",
              "MOSS",
              "bigstatsr",
              "BiocManager",
              "ggplot2",
              "viridis",
              "Rmisc",
              "MOSim")

# Check the required packages are installed.
is_installed <- req_pack %in% rownames(installed.packages())
names(is_installed) <- req_pack

# Install packages from CRAN if needed.
install.packages(req_pack[!is_installed & req_pack != "MOSim"])

# Install packages from Bioconductor if needed.
if (!is_installed[10]) BiocManager::install("MOSim")
