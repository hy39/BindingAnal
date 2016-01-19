
# read the csv file
M <- read.csv('Mlist.dat')
Mnew <- read.csv('Mnewlist.dat')

# calculate the mean
means <- aggregate(M[,4] ~  M[,5], M, mean)

# plot the netcharge data
ggplot(M,aes(factor(M$X.1.1), M$X.1.0291)) +
scale_size_area() +
xlab('Netcharge Changes') +
ylab('The effect of amino acid change on cell binding (log scale)') +
geom_boxplot() + 
geom_jitter(width=.2) + 
stat_summary(fun.y=mean, colour="darkred", geom="point", shape=1, size=3, show.legend = FALSE) +
scale_x_discrete(breaks=c(-1, 0, 1), labels=c("Negative", "Neutral", "Positive")) +
theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
