library("miceadds")
library("ggplot2")
library("RColorBrewer")

wd <- getwd()
setwd("results1")
files <- list.files()[grep("omic_int_meth_perf_data1",list.files())]
RES1 <- do.call("rbind",lapply(X=files,FUN=load.Rdata,objname="RES1"))
setwd(wd)

PROP_SIGNAL <- unique(RES1$prop_signal)
METHOD <- unique(RES1$method)

D1 <- NULL
for (prop_signal in PROP_SIGNAL) {
  for (method in METHOD) {
    tmp <- RES1$prop_signal == prop_signal & RES1$method == method
    D1 <- rbind(D1,
                data.frame(
                  method=method,
                  prop_signal=prop_signal,
                  SPEC = rowMeans(do.call("cbind",
                                          tapply(RES1$SPEC[tmp],
                                                 RES1$rep[tmp],
                                                 c))),
                  SENS =rowMeans(do.call("cbind",
                                         tapply(RES1$SENS[tmp],
                                                RES1$rep[tmp],
                                                c))),
                  sdSPEC = apply(do.call("cbind",
                                         tapply(RES1$SPEC[tmp],
                                                RES1$rep[tmp],
                                                c)),1,var),
                  sdSENS = apply(do.call("cbind",
                                         tapply(RES1$SENS[tmp],
                                                RES1$rep[tmp],
                                                c)),1,var)))

  }

}

D1$sd <- (sqrt(D1$sdSPEC) + sqrt(D1$sdSENS))/10
D1$prop_signal <- factor(D1$prop_signal,
                         levels=sort(PROP_SIGNAL),
                         labels=c("Number of \nsignal features = 34 (0.1*p)",
                                  "Number of \nsignal features = 274 (0.8*p)"))

col <- c(viridis::viridis(4,option = "B"),
         viridis::viridis(4,option = "D")[-1])

meth_col <- c("iCluster"=col[6],
              "NMF"=col[2],
              "OmicPLS"=col[4],
              "SNFtool"=col[3],
              "mixOmics"=col[5],
              "MOSS"=col[1])
D1$method <- factor(D1$method,
                    levels=c("iCluster","NMF","OmicsPLS",
                             "SNF" ,"mixOmics","MOSS"),
                    labels = names(meth_col),ordered=TRUE)

p <- ggplot(data = D1,aes(x=SPEC,y=SENS,group=method,color=method)) +
  facet_grid(.~prop_signal)+
  geom_line(size=1.1)+
  geom_ribbon(aes(ymin=SENS - sd, ymax=SENS + sd,fill=method),colour = NA)+
  scale_color_manual(values=meth_col)+
  scale_fill_manual(values=alpha(meth_col,0.2))+
  scale_x_continuous(breaks = seq(-1,0,by=0.25),
                     labels = seq(1,0,by=-0.25)) +
  scale_y_continuous(breaks = seq(0,1,by=0.25))+
  labs(y="Sensitivity",x="Specificity",
       fill="Method",
       colour="Method")+
  geom_segment(aes(x = -1, xend = 0, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")+
  theme_bw()+
  theme(legend.position=c(0.95,0.1),
        legend.justification=c(0.95,0.1),
        axis.title= element_text(colour = "black", face = "bold"),
        strip.text.x = element_text(colour = "black", face = "bold") )

# Save plot.
ggsave(p,filename = "Figure1-A.png",device = "png",
       width = 5,height = 5,units = "in")

q(save="no")
