##-----------------------------------------------------------------
## -The following lines, preceded by '##', represent a suggestion to
## run the R code below in a cluster using SLURM-

## #!/usr/bin/bash --login
## #!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper

## #SBATCH --job-name=omic_integration
## #SBATCH --time=10:00:00
## #SBATCH --mem=30G
## #SBATCH --array=101-500

## #SBATCH --output=/mnt/research/quantgen/logs/%A_%a

## cd $SLURM_SUBMIT_DIR

## if [ -d "output_files" ]
## then
## echo "Storing output results in output_files/"
## else
##   mkdir output_files
## fi


## R CMD BATCH --no-save --no-restore '--args jobID='${SLURM_ARRAY_TASK_ID} omic_integration_comparisson.R output_files/x_${SLURM_ARRAY_TASK_ID}.out


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
if (!exists("jobID")) jobID <- 1

# Load required packages.
library("iCluster")
library("NMF")
library("SNFtool")
library("OmicsPLS")
library("mixOmics")
library("MOSS")
library("pROC")

# Use gene expression from breast tumors from TCGA.
data("breast.TCGA")

PE <- rbind(breast.TCGA$data.train$protein)
GE <- rbind(breast.TCGA$data.train$mrna)

# Get proportion of features with signal.
PROP_SIGNAL <- c(0.1,0.8)

# Create object(s) to store results.
RES1 <- NULL

# Create a response omic.
set.seed(7*jobID+jobID*3)

# Shuffle features.
signal_features_PE <- sample(seq_len(ncol(PE)))
signal_features_GE <- sample(seq_len(ncol(GE)))

# Break correlations among features and generate noise.
#i <- sample(seq_len(nrow(PE)),round(nrow(PE)*0.8))
#PE <- PE[i,]
#GE <- GE[i,]
PE <- scale(PE)
GE <- scale(GE)
PE[is.na(PE)] <- 0
GE[is.na(GE)] <- 0
PE_perm <- apply(PE,2,sample)
GE_perm <- apply(GE,2,sample)

for (prop_signal in PROP_SIGNAL) {

  # Get proportion of signal.
  prop_signal_PE <- round(prop_signal * ncol(PE))
  prop_signal_GE <- round(prop_signal * ncol(GE))

  # Set features labels.
  features_labels_PE <- rep("Backround",ncol(PE))
  features_labels_PE[signal_features_PE[seq_len(prop_signal_PE)]] <- "Signal"

  PE_perm[, signal_features_PE[seq_len(prop_signal_PE)]] <-
    PE[,signal_features_PE[seq_len(prop_signal_PE)]]

  features_labels_GE <- rep("Backround",ncol(GE))
  features_labels_GE[signal_features_GE[seq_len(prop_signal_GE)]] <- "Signal"

  GE_perm[, signal_features_GE[seq_len(prop_signal_GE)]] <-
    GE[,signal_features_GE[seq_len(prop_signal_GE)]]

  features_labels <- c(features_labels_PE,features_labels_GE)
  true <- ifelse(features_labels == "Signal", 1, 0)

  # Fit MOSS.
  out_moss <- MOSS::moss(data.blocks = list(PE_perm,GE_perm),
                         method = "pca",nu.v=100,alpha.v = 0,
                         verbose = TRUE,norm.arg = TRUE,
                         use.fbm = TRUE,
                         K.X=2)

  # Get ROCs.
  AUC <- smooth(roc(predictor = out_moss$sparse$v[,1]^2 ,
                    response = true),
                method="density",n=500)

  # Store performance results.
  res <- data.frame(method="MOSS",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)
  RES1 <- rbind(RES1,res)

  # Fit mixOmics.
  out_mixOmics <- mixOmics::pca(X = cbind(PE_perm,GE_perm),
                                ncomp = 2)

  AUC <- smooth(roc(predictor=out_mixOmics$loadings$X[,1]^2,
                    response=ifelse(true,
                                    1,0)),
                method="density",n=500)
  # Store performance results.
  res <- data.frame(method="mixOmics",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)

  RES1 <- rbind(RES1,res)

  # Fit iCluster.
  out_iCluster <- iCluster::iCluster2(datasets = list(PE_perm,GE_perm),
                                      k = 4,
                                      lambda = c(1e3,1e3))

  AUC <- smooth(roc(predictor=out_iCluster$W[,1]^2,
                    response=ifelse(true,
                                    1,0)),
                method="density",n=500)

  # Store performance results.
  res <- data.frame(method="iCluster",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)

  RES1 <- rbind(RES1,res)

  # Fit SNF.
  blocks <- lapply(list(PE_perm,GE_perm),
                   standardNormalization)

  K <-  20;		# number of neighbors,
  alpha <-  0.5;  	# hyperparameter.
  T0 <-  20;
  Dist1 <-  (dist2(blocks[[1]],blocks[[1]]))^(1/2)
  Dist2 <-  (dist2(blocks[[2]],blocks[[2]]))^(1/2)

  #Similarity graphs.
  W1 <-  affinityMatrix(Dist1, K, alpha)
  W2 <-  affinityMatrix(Dist2, K, alpha)

  # Rank features.
  out_snf <- rankFeaturesByNMI(blocks,
                               SNF(list(W1,W2), K, T0))
  AUC <- smooth(roc(predictor=unlist(out_snf[[1]]),
                    response=ifelse(true,
                                    1,0)),
                method="density",n=500)
  # Store performance results.
  res <- data.frame(method="SNF",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)

  RES1 <- rbind(RES1,res)

  # Fit NMF.
  out_nmf <- nmf(rbind(t(ifelse(PE_perm < 0, 0, PE_perm)),
                       t(ifelse(GE_perm < 0, 0, GE_perm))),
                 rank=3,
                 method='lee')

  AUC <- smooth(roc(predictor=featureScore(out_nmf)^2,
                    response=ifelse(true,
                                    1,0)),
                method="density",n=500)
  # Store performance results.
  res <- data.frame(method="NMF",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)

  RES1 <- rbind(RES1,res)

  # Fit OmicsPLS
  o2m_out <- o2m(X = GE_perm,
                 Y = PE_perm,
                 n=2,
                 nx=2,
                 ny=2)

  # Measure variable selection performance.
  V <- rbind(o2m_out$`C.`,o2m_out$`W.`)

  AUC <- smooth(roc(predictor=rowMeans(V)^2,
                    response=ifelse(true,
                                    1,0)),
                method="density",n=500)

  # Store performance results.
  res <- data.frame(method="OmicsPLS",
                    prop_signal=prop_signal,
                    SENS=AUC$sensitivities,
                    SPEC=-AUC$specificities,
                    rep=jobID)

  RES1 <- rbind(RES1,res)
}

save(RES1,
     file=paste0("results/omic_int_meth_perf_data1",jobID,".RData"))

q(save="no")
