# Load required packages.
library("ggplot2")
library("viridis")

# Load results.
files <- list.files()
data <- NULL
for (f in files) {
  load(f, verbose = TRUE)
  data <- rbind(data, out)
}

# Turn seconds into hours.
data$timehours <- as.numeric(data$time)/3600
lev<-round(unique(data$timehours),digits = 3)

# Have a plot without missing bars
lev<-lev[order(lev)]

# Round values.
data$timehours<-round(data$timehours,3)

# Set Method as a factor.
data$Model<-factor(data$Model)

# Set time as factors.
data$time1 <- actor(data$timehours,
                    levels=as.character(unique(lev)))

# Relabel methods.
data$Model<-factor(data$Model,
                   levels=c('iCluster', 'NMF', 'OmicsPLS',
                            'SNF', 'SPCA','MOSS_NoTun',
                            'MOSS_NoTun_big'),
                   labels = c('iCluster', 'NMF', 'OmicsPLS',
                              'SNFtool', 'mixOmics',
                              'MOSS\n(reg. matrices)',
                              'MOSS\n(FBM)'),
                   ordered = TRUE)

# Define size as ordered factors.
data$n <- factor(data$n,levels = unique(data$n),
                 labels = c("n=1e2","n=1e3","n=1e4","n=1e5"),
                 ordered = T)
data$p <- factor(data$p,
                 levels = c(1e+03, 1e+04, 1e+05, 1e+06),
                 labels = c("p=1e3","p=1e4","p=1e5","p=1e6"),
                 ordered = T)

# Exclude missing times for iCluster.
tmp <- data$Model == "iCluster" &
  ((data$n == "n=1e5" & data$p == "p=1e5") |
     (data$n == "n=1e4" & data$p == "p=1e6"))

# Subset results for plot.
data <- data[!tmp,]

# Use "Method" instead of "Model".
data$Method <- data$Model

# Create plot.
plot <- ggplot(data,aes(x=Model,
                        y=time1,
                        fill=Method,
                        group=Method)) +
  facet_grid(p~n,switch = "y",scales = "free_y")+
  geom_bar(stat = "identity",position = position_dodge(0.9))+
  scale_fill_manual(values= c(viridis(4,option = "B"),
                              viridis(4,option = "D")[-1])[order(runif(7))])+
  theme_bw()+
  scale_x_discrete("Omic integration method") +
  scale_y_discrete("Time (hours)")+
  theme(strip.text = element_text(face="bold"),
        legend.key.size = unit(.75, "cm"),
        axis.text.x = element_blank(),
        axis.title = element_text(face="bold"))

# Save plot.
ggsave(plot,file= "Figure1-B.png",width = 9,height = 10,units = "in")

# Exit R without saving.
q(save = "no")
