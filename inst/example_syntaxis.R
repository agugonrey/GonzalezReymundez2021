# Load required packages.
library("MOSim")
library("MOSS")

# Set random seed.
set.seed(340)

# Set types of omic blocks to simulate.
omics_list <- c("RNA-seq", "miRNA-seq","DNase-seq")

# Set simulation options by omic block.
omics_options <- c(

  omicSim(#Name of the omic.
    "RNA-seq",
    # Limit the number of features to simulate
    totalFeatures = 1000),

  omicSim("miRNA-seq",
          # Limit the number of features to simulate
          totalFeatures = 50,
          # Modify the percentage of regulators with effects.
          regulatorEffect = list(
            'activator' = 0.0,
            'repressor' = 0.2,
            'NE' = 0.8
          )),
  omicSim("DNase-seq",
          # Limit the number of features to simulate
          totalFeatures = 2000,
          # Modify the percentage of regulators with effects.
          regulatorEffect = list(
            'activator' = 0.0,
            'repressor' = 0.2,
            'NE' = 0.8)))

# Create multi-omic data.
multi_simulation <- mosim(omics = omics_list,
                          numberReps = 100,
                          numberGroups = 3,
                          times = 1,
                          omicsOptions = omics_options)

# Turning simulated omics into MOSS' inputs.
out_sim <- lapply(omicResults(multi_simulation), t)

# Get simulation information.
out_settings <- omicSettings(multi_simulation,
                             association = FALSE,
                             only.linked = TRUE)

# Run example.
set.seed(345)
out_moss <- moss(data.blocks = out_sim,
                 method = "pls",
                 resp.block = 1,
                 scale.arg = TRUE,
                 norm.arg = TRUE,
                 K.Y = 3,
                 nu.v = seq(1,200,by=2),
                 nu.u = seq(1,100,by=2),
                 alpha.v = 0.5,
                 alpha.u = 0.5,
                 K.X = 50,
                 use.fbm = TRUE,
                 nu.parallel = TRUE,
                 tSNE = TRUE,
                 cluster = TRUE,
                 plot =TRUE)
