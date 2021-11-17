#-----------------------------------------------------------------------
# create_Figure3.R:

#This script uses the results of analytic_performance_moss.R
# to create Figure 3.

# Load required packages.
library("ggplot2")
library("viridis")
library("Rmisc")

# Retrieve results and store them in 'D'.
files <- list.files()
D <- NULL
for (f in files) D <- rbind(D,res)

# Get simulations summaries.
D <- summarySE(D,
               measurevar="Perf_value",
               groupvars=c("Perf_met","n","p","Signal"))

# Set number of features as factors.
D$p <- factor(D$p,levels=c(100,1000,10000),
              labels=paste0("p = ",c("100","1,000","10,000")),
              ordered=T)

# Set number of samples as factors.
D$n <- factor(D$n,levels=c(1000,10000),labels=c("1,000","10,000"),
              ordered=T)

# Set signal as factors.
D$Signal <- as.factor(D$Signal)

# Same for performance variables.
D$Perf_met <- factor(D$Perf_met,
                     levels=unique(D$Perf_met),
                     labels=c("ACC","SEN","SPE","PRE"),ordered=T)

# Create plot.
p <- ggplot(D, aes(x=n,y=Perf_value,fill=Signal)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin=Perf_value-sd, ymax=Perf_value+sd,color=Signal), width=.25,position=position_dodge(.9))+
  scale_fill_manual(values = viridis(n = 6,alpha = 0.6)[c(2,5)])+
  scale_color_manual(values = viridis(n = 6,alpha = 1)[c(2,5)])+
  facet_grid(Perf_met~p,switch="y",scales="free_x")+
  scale_y_continuous("Performance value (SD)")+
  scale_x_discrete("Number of samples (n)")+
  theme_bw()+
  theme(
    text = element_text(size=15),
    axis.text.x = element_text(vjust = 1,hjust=1,angle = 0))

# Save plot.
ggsave(p,filename = "Figure3.png",device = "png",
       width = 6,height = 5,units = "in")

# Exit R without saving.

q(save = "no")

