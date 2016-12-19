setwd("E:/Documents/Github/BindingAnal/validate_netcharge")
install.packages("ggplot2")
library('ggplot2');
# read the csv file
#M <- read.csv('RDE/MHensley.dat')
#M <- read.csv('RDE/MOthers.dat')
M <- read.csv('RDE/MTotal.dat')

#M <- read.csv('dat_netch_bindn.csv')
#Mnew <- read.csv('Mnewlist.dat')

# calculate the mean
means <- aggregate(M[,4] ~  M[,5], M, mean)

# plot the netcharge data
ggplot(M,aes(factor(M[,5]), M[,4])) +
  #scale_size_area() +
  xlab('Changes of the net charge') +
  ylab('The relative binding avidity (log scale)') +
  geom_boxplot(outlier.colour=NA) + 
  geom_jitter(width=.5, height=0, size=2) + 
  stat_summary(fun.y=mean, colour="darkred", geom="point", shape=1, size=3, show.legend = FALSE) +
  scale_x_discrete(breaks=c(-1, 0, 1), labels=c("Decrease", "Same", "Increase")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))


# anova for netcharge
fm1 <- aov(M[,4] ~ M[,5])
anova(fm1)