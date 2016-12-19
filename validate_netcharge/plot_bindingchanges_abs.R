setwd("E:/Documents/Github/BindingAnal/validate_netcharge")
install.packages("ggplot2")
library('ggplot2');
# read the csv file
M <- read.csv('RDE/abs_charge.csv')

# plot the netcharge data
ggplot(M,aes(factor(M[,1]), M[,2])) +
  #scale_size_area() +
  xlab('Netcharge Changes') +
  ylab('The Ratio of Cell Binding (log scale)') +
  geom_boxplot(outlier.colour=NA) + 
  geom_jitter(width=.5, height=0, size=2) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point", shape=1, size=3, show.legend = FALSE) +
  scale_x_discrete(breaks=c(-1, 0, 1), labels=c("Decreased", "Same", "Increased")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
