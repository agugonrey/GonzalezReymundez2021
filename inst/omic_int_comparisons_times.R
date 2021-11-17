##-----------------------------------------------------------------
## -The following lines, preceded by '##', represent a suggestion to
## run the R code below in a cluster using SLURM-

##!/usr/bin/bash --login

##SBATCH --job-name=omicIntTimes
##SBATCH --time=70:30:00
##SBATCH --nodes=1
##SBATCH --cpus-per-task=20
##SBATCH --mem=100gb
##SBATCH --array=1-112

#cd $SLURM_SUBMIT_DIR

#R CMD BATCH --no-save --no-restore
# '--args jobID='${SLURM_ARRAY_TASK_ID}
# omic_int_comparisons_times.R x_times_out_${SLURM_ARRAY_TASK_ID}

#------------------------------------------------------
# omic_int_comparisons_times.R:

#This script compares MOSS against several
# omic integration packages in terms of
#computational time.

args <- commandArgs(TRUE)

# Evaluating arguments.
for(i in seq_along(args)){
  eval(parse(text=args[[i]]))
}

# Get job ID.
if(!exists("jobID")){jobID <- 1}

# Load required packages.
library("iCluster")
library("NMF")
library("SNFtool")
library("OmicsPLS")
library("mixOmics")
library("MOSS")
library("bigstatsr")

# Simulate omic blocks
sim_blocks <- simulate_data(moss=42)$sim_blocks[-4]

# Define scenarios in terms of number of samples and features.
size<-expand.grid(n=c(1e2,1e3,1e4,1e5),
                  p=c(1e3,1e4,1e5,1e6),
                  Method=c("spca","moss",
                           "moss_big","iCluster",
                           "o2m","nmf","snf"))[jobID, ]
# Get number of samples.
n <- as.numeric(size[1])

# Get number of features.
p <- as.numeric(size[2])

# Get omic integration method.
Method <- size[3]

# Create an FBM for Method = 'moss_big'.
if (Method == "moss_big") {
  sim_blocks_ext <- lapply(1:3, function(i) {
    m <- FBM(nrow =  n, ncol =  floor(p/3))
    m[1:100,1:200] <- sim_blocks[[i]][1:100,1:200]
    return(m)
  })
}

# Create regular matrices for every other method.
if (Method !="moss_big") {
  sim_blocks_ext <- lapply(1:3, function(i) {
    m <- matrix(rnorm(p*n), n, floor(p/3))
    m[1:100,1:200] <- sim_blocks[[i]][1:100,1:200]
    return(m)
  })
}

# Calculate computational time by method.

if(Method=="spca"){
  x <- do.call(what = "cbind",
               args = sim_blocks_ext )

  ptm<-proc.time()
  out_spca <- try(spca(X = x,
                       ncomp = 2,
                       keepX = c(479,479)),silent = TRUE)
  timeSPCA<-proc.time()-ptm
  if (class(out_spca) != "try-error")
    out <- data.frame(n=n,p=p,Model="SPCA",time=c(timeSPCA[3]))
  else out<-data.frame(n=n,p=p,Model="SPCA",time=NA)

  save(out, file=paste0("Size_",n,"_",p,"SPCA.RData"))
}

if(Method=="moss"){
  ptm<-proc.time()
  out_moss_NoTun <- try(moss(data.blocks = sim_blocks_ext ,
                             method = "pca",
                             scale.arg = F,norm.arg = F,
                             K.X = 2,
                             nu.v = 479,
                             alpha.v = 1),silent = TRUE)
  timeMOSSNoTun<-proc.time()-ptm

  if (class(out_moss_NoTun) != "try-error")
    out<-data.frame(n=n,p=p,Model=c("MOSS_NoTun"),
                    time=c(timeMOSSNoTun[3]))
  else out<-data.frame(n=n,p=p,Model=c("MOSS_NoTun"),time=NA)
  save(out, file=paste0("Size_",n,"_",p,"mossnotun.RData"))
}

if(Method=="moss_big"){

  ptm<-proc.time()
  out_moss_NoTun_big <-
    try(moss(data.blocks = sim_blocks_ext ,
                                 method = "pca",
                                scale.arg = F,
                                 norm.arg = F,
                                 K.X = 2,
                                 nu.v = 479,
                                 alpha.v = 1),silent = TRUE)
  timeMOSSNoTunBig<-proc.time()-ptm

  if (class(out_moss_NoTun_big) != "try-error")
    out<-data.frame(n=n,p=p,Model=c("MOSS_NoTun_big"),
                    time=c(timeMOSSNoTunBig[3]))
  else out<-data.frame(n=n,p=p,Model=c("MOSS_NoTun_big"),time=NA)
  save(out, file=paste0("Size_",n,"_",p,"mossnotunBig.RData"))
}

if(Method=="iCluster"){
  ptm<-proc.time()
  out_icluster <- try(iCluster::iCluster2(sim_blocks_ext ,
                                          k = 3,
                                          scalar = TRUE),silent = TRUE)
  timeIcluster<-proc.time()-ptm
  if (class(out_icluster) != "try-error")
    out<-data.frame(n=n,p=p,Model=c("iCluster"),
                    time=c(timeIcluster[3]))
  else out<-data.frame(n=n,p=p,Model=c("iCluster"),time=NA)
  save(out, file=paste0("Size_",n,"_",p,"iCluster.RData"))
}

if (Method=="snf") {
  out_sim <- lapply(sim_blocks_ext, standardNormalization)

  ptm<-proc.time()
  K <- 20;		# Number of neighbors.
  alpha <- 0.5;  	# Hyperparameter.
  T0 <- 20;
  Dist1 <-  (dist2(out_sim[[1]],out_sim[[1]]))^(1/2)
  Dist2 <-  (dist2(out_sim[[2]],out_sim[[2]]))^(1/2)
  Dist3 <-  (dist2(out_sim[[3]],out_sim[[3]]))^(1/2)

  #Similarity graphs.
  W1 <-  affinityMatrix(Dist1, K, alpha)
  W2 <-  affinityMatrix(Dist2, K, alpha)
  W3 <-  affinityMatrix(Dist3, K, alpha)

  # Rank features.
  out_snf <- try(SNF(list(W1,W2,W3), K, T0),
                 silent = TRUE)
  timeSnf <-proc.time()-ptm

  out<-data.frame(n=n,p=p,Model=c("SNF"),time=c(timeSnf[3]))
  save(out, file=paste0("Size_",n,"_",p,"snf.RData"))
}

if (Method=="nmf"){
  ptm<-proc.time()
  sim_blocks_ext <- lapply(sim_blocks_ext,
                           function(x) ifelse(x > 0,x,-x))
  out_nmf <- try(nmf(do.call("rbind",lapply(sim_blocks_ext,t)),
                     rank=2, method='pe',
                     alpha=0.01, beta=1),silent = TRUE)
  timeNmf<-proc.time()-ptm

  if (class(out_nmf) != "try-error")
    out<-data.frame(n=n,p=p,Model=c("NMF"),time=c(timeNmf[3]))
  else  out<-data.frame(n=n,p=p,Model=c("NMF"),time=NA)
  save(out, file=paste0("Size_",n,"_",p,"nmf.RData"))
}

if (Method=="o2m") {
  ptm<-proc.time()
  out_o2m <- try(o2m(X=scale(do.call("cbind",sim_blocks_ext[-1])),
                     Y = scale(sim_blocks_ext[[1]]),
                     n=2,
                     nx=2,
                     sparse = TRUE,
                     keepx = 479,
                     keepy = 479,
                     ny=2),silent = TRUE)
  timeO2m<-proc.time()-ptm

  if (class(out_o2m) != "try-error")
    out<-data.frame(n=n,p=p,Model=c("OmicsPLS"),time=c(timeO2m[3]))
  else out<-data.frame(n=n,p=p,Model=c("OmicsPLS"),time=NA)
  save(out, file=paste0("Size_",n,"_",p,"o2m.RData"))
}

# Exit R without saving.
q(save="no")

