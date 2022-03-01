library("miceadds")
library("ggplot2")
library("RColorBrewer")

wd <- getwd()
setwd("results")
files <- list.files()[grep("moss_perf_real",list.files())]
RES1 <- do.call("rbind",lapply(X=files,FUN=load.Rdata,objname="RES1"))
setwd(wd)

RES1 <- RES1[RES1$prop_signal > 0.01, ]
PROP_SIGNAL <- unique(RES1$prop_signal)

D1 <- NULL
for (prop_signal in PROP_SIGNAL) {

  tmp <- RES1$prop_signal == prop_signal
  D1 <- rbind(D1,
              data.frame(
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

D1$sd <- (sqrt(D1$sdSPEC) + sqrt(D1$sdSENS))
D1$prop_signal <- factor(D1$prop_signal,
                         levels=sort(PROP_SIGNAL),
                         labels=c("6,011 (0.1p)", "48,090 (0.8p)"))

p <- ggplot(data = D1,aes(x=SPEC,y=SENS,group=prop_signal,color=prop_signal)) +

  geom_line(size=1.1)+
  geom_ribbon(aes(ymin=SENS - sd, ymax=SENS + sd,fill=prop_signal),colour = NA)+
  scale_color_manual(values=alpha(c("black","darkgrey"),0.95))+
  scale_fill_manual(values=alpha(c("black","darkgrey"),0.2))+
  scale_x_continuous(breaks = seq(-1,0,by=0.25),
                     labels = seq(1,0,by=-0.25)) +
  scale_y_continuous(breaks = seq(0,1,by=0.25))+
  labs(y="Sensitivity",x="Specificity",
       fill="Number of signal\nfeatures",
       colour="Number of signal\nfeatures")+
  geom_segment(aes(x = -1, xend = 0, y = 0, yend = 1),
               color="darkgrey", linetype="dashed")+
  theme_bw()+

  theme(legend.position=c(0.95,0.5),legend.justification=c(0.95,0.5),
        axis.title= element_text(colour = "black", face = "bold"),
        strip.text.x = element_text(colour = "black", face = "bold"))

# Save plot.
ggsave(p,filename = "Figure1-C.png",device = "png",
       width = 5,height = 5,units = "in")

q(save="no")
