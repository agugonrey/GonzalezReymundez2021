args <- commandArgs(TRUE)

# Evaluating arguments.
for(i in seq_along(args)){
  eval(parse(text=args[[i]]))
}

# Get job ID.
if(!exists("jobID")) jobID <- 1

# Load required packages.
library("MOSS")
library("pROC")
library("bigstatsr")

# Get candidate omics
# Download  from https://data.mendeley.com/datasets/r8p67nfjc8/1.
load("GE.RData")
load("METH.RData")
load("CNV.RData")

# Get proportion of features with signal.
PROP_SIGNAL <- c(0.1,0.8)

# Create object(s) to store results.
RES1 <- NULL

# Create a response omic.
set.seed(7*jobID+jobID*3)

# Shuffle features.
signal_features_GE <- sample(seq_len(ncol(GE)))
signal_features_METH <- sample(seq_len(ncol(METH)))
signal_features_CNV <- sample(seq_len(ncol(CNV)))

# Break correlations among features and generate noise.
i <- sample(seq_len(nrow(GE)),round(nrow(GE)*0.8))
GE_perm <- apply(GE[i,],2,sample)
METH_perm <- apply(METH[i,],2,sample)
CNV_perm <- apply(CNV[i,],2,sample)

for (prop_signal in PROP_SIGNAL) {

  # Get proportion of signal by omic.
  prop_signal_GE <- round(prop_signal * ncol(GE))
  prop_signal_METH <- round(prop_signal * ncol(METH))
  prop_signal_CNV <- round(prop_signal * ncol(CNV))

  # Set features labels.
  features_labels_GE <- rep("Backround",ncol(GE))
  features_labels_METH <- rep("Backround",ncol(METH))
  features_labels_CNV <- rep("Backround",ncol(CNV))
  features_labels_GE[signal_features_GE[seq_len(prop_signal_GE)]] <- "Signal"
  features_labels_METH[signal_features_METH[seq_len(prop_signal_METH)]] <- "Signal"
  features_labels_CNV[signal_features_CNV[seq_len(prop_signal_CNV)]] <- "Signal"

  GE_perm[, signal_features_GE[seq_len(prop_signal_GE)]] <-
    GE[i,signal_features_GE[seq_len(prop_signal_GE)]]
  METH_perm[, signal_features_METH[seq_len(prop_signal_METH)]] <-
    METH[i,signal_features_METH[seq_len(prop_signal_METH)]]
  CNV_perm[, signal_features_CNV[seq_len(prop_signal_CNV)]] <-
    CNV[i,signal_features_CNV[seq_len(prop_signal_CNV)]]

  # Fit MOSS.
  out_moss <- moss(data.blocks = list(GE_perm,METH_perm,CNV_perm),
                   method = "pca",
                   verbose = TRUE,
                   use.fbm = TRUE,
                   scale.arg = TRUE,norm.arg=TRUE,
                   K.X=50)

  # Measure performance for different degrees of sparsity.
  pred <- out_moss$dense$v[,1]^2
  features_labels <- c(features_labels_GE,features_labels_METH,features_labels_CNV)
  true <- ifelse(features_labels == "Signal", 1, 0)

  # Get ROCs.
  AUC <- smooth(roc(predictor = pred / sum(pred) ,
                    response = true),
                method="density",n=500)

  # Store performance results.
  res <- data.frame(
    prop_signal=prop_signal,
    SENS=AUC$sensitivities,
    SPEC=-AUC$specificities,
    rep=jobID)
  RES1 <- rbind(RES1,res)
}

save(RES1,
     file=paste0("results/moss_perf_real_dat_",jobID,".RData"))

q(save="no")

