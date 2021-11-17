##-----------------------------------------------------------------
## -The following lines, preceded by '##', represent a suggestion to
## run the R code below in a cluster using SLURM-

##!/usr/bin/bash --login

##SBATCH --job-name=mossAnPerf
##SBATCH --time=20:00:00
##SBATCH --nodes=1
##SBATCH --cpus-per-task=20
##SBATCH --mem=100gb
##SBATCH --array=1-10000

#cd $SLURM_SUBMIT_DIR

#R CMD BATCH --no-save --no-restore
# '--args jobID='${SLURM_ARRAY_TASK_ID}
# analytic_performance_moss.R x_perf_out_${SLURM_ARRAY_TASK_ID}

#------------------------------------------------------
# analytic_performance_moss.R:

#This script estimates MOSS' performance at
# detecting features with signal.

args <- commandArgs(TRUE)

# Evaluating arguments.
for(i in seq_along(args)){
  eval(parse(text=args[[i]]))
}

# Get job ID.
if(!exists("jobID")){jobID <- 1}

# Load required packages.
library("MOSim")
library("MOSS")

# Set simulation scenarios.
scenarios <- expand.grid(n=c(1e3,1e4),
                         p=c(1e2,1e3,1e4),
                         signal=c(0.05,0.2),
                         rep=seq_len(1e3))

# Exclude the largest so MOSim can run.
scenarios <- scenarios[!(scenarios$n == 1e4 &
                           scenarios$p == 1e4), ]

# Set scenario and random repetition.
set.seed(34*scenarios[jobID, "rep"])

# Set types of omic blocks to simulate.
omics_list <- c("RNA-seq", "miRNA-seq","DNase-seq")

# Set simulation options by omic block.
omics_options <- c(

  omicSim(#Name of the omic.
    "RNA-seq",
    # Limit the number of features to simulate
    totalFeatures = scenarios[jobID, "p"]),

  omicSim("miRNA-seq",
          # Limit the number of features to simulate
          totalFeatures = scenarios[jobID, "p"],
          # Modify the percentage of regulators with effects.
          regulatorEffect = list(
            'activator' = 0.0,
            'repressor' = scenarios[jobID, "signal"],
            'NE' = 1 - scenarios[jobID, "signal"]
          )),
  omicSim("DNase-seq",
          # Limit the number of features to simulate
          totalFeatures = scenarios[jobID, "p"],
          # Modify the percentage of regulators with effects.
          regulatorEffect = list(
            'activator' = 0.0,
            'repressor' = scenarios[jobID, "signal"],
            'NE' = 1 - scenarios[jobID, "signal"])))

# Create multi-omic data.
multi_simulation <- mosim(omics = omics_list,
                          numberReps = scenarios[jobID, "n"],
                          numberGroups = 3,
                          times = 1,
                          omicsOptions = omics_options)

# Turning simulated omics into MOSS inputs.
out_sim <- lapply(omicResults(multi_simulation), t)

# Get simulation information.
out_settings <- omicSettings(multi_simulation,
                             association = FALSE,
                             only.linked = TRUE)

# Define  background features (i.e. noise).
features_labels <- rep("Background",
                       3 * scenarios[jobID, "p"])

# Define regulatory modules (i.e. signal).
features_labels[unlist(lapply(out_sim,colnames)) %in%
                  unique(c(out_settings[[1]]$ID[out_settings[[1]]$DE],
                           out_settings[[2]]$ID,
                           out_settings[[3]]$ID))] <- "Signal"

# Run MOSS
out_moss <- moss(data.blocks = out_sim,
                 nu.v = round(seq(1,3*scenarios[jobID, "p"])),
                 alpha.v = 0.5,
                 axes.pos = c(1, 2),
                 exact.dg = TRUE,
                 use.fbm = TRUE,
                 nu.parallel = TRUE,
                 cluster = list(eps_range=c(0,1),
                                    eps_res=10,
                                    min_clus_size=2),
                 plot = TRUE,
                 lib.thresh = FALSE)

# Measure variable selection performance.
tmp <- out_moss$feat_signatures$signatures$Cluster != 0 &
  out_moss$feat_signatures$signatures$Candidate == TRUE

pred_pos <- rep(FALSE, length(features_labels))
pred_pos[unlist(lapply(out_sim, colnames)) %in%
           out_moss$feat_signatures$signatures$Feature_name[tmp]] <- TRUE
pred_neg <- !pred_pos
true_pos <- features_labels == "Signal"
true_neg <- !true_pos
TP <- sum(pred_pos & true_pos)
FP <- sum(pred_pos & true_neg)
TN <- sum(pred_neg & true_neg)
FN <- sum(pred_neg & true_pos)

perf <- c("ACC" = (TP + TN) / (TP + TN + FP + FN),
          "PRE" = 1 - (FP / (FP + TP)),
          "SPE" = TN / (TN + FP),
          "SEN" = TP / (TP + FN))

# Store results.
res <- data.frame("n"=scenarios[jobID, "n"],
                  "p"=scenarios[jobID, "p"],
                  "Signal"=scenarios[jobID, "signal"],
                  "Perf_met"=names(perf),
                  "Perf_value"=perf,
                  "Repetition"=scenarios[jobID, "rep"])
# Save results.
save(res,
     file=paste0("moss_perf_",jobID,".RData"))

# Exit R without saving.
q(save = "no")
