# Matlablimits

Ch <- read.csv('Chlist.dat')

#Ch1 <- Ch[which(Ch[,4]==10),]
Ch1 <- Ch
means <- aggregate(Ch1[,3] ~  Ch1[,10], Ch1, mean)
sds <- aggregate(Ch1[,3] ~  Ch1[,10], Ch1, sd)
se <- sds[,2]

age = factor(means[,1])
netcharge = means[,2]

df <- data.frame(
  age = age,
  netcharge = netcharge
)

x1 <- Ch1[,10]
y1 <- Ch1[,3] 
limits <- aes(ymin=netcharge-se, ymax=netcharge+se)


# R
#ggplot(aes(df, x=df$age, y=df$netcharge)) +
ggplot(df, aes(df, x=age, y=netcharge)) +
	#geom_bar(position=position_dodge(), stat="identity") +
geom_point(size = 5) +
geom_errorbar(limits, width=0.25) +
       # geom_errorbar(limits, position="dodge", width=0.25) +
        ylim(16, 20)



## use the following code

df <- data.frame (
  age = factor(c(1,2,3,4,5,6)),
  highcharge = c(0.40, 0.29, 0.41, 0.39, 0.62, 0.70),
  lb = c(0.34, 0.22, 0.33, 0.31, 0.53, 0.62),
  ub = c(0.48, 0.39, 0.50, 0.49, 0.71, 0.77)
)

limits <- aes(ymin=df$lb, ymax=df$ub)

ggplot(df, aes( x=age, y=highcharge)) +
geom_point(size = 5) +
geom_errorbar(limits, width=0.25) +
scale_x_discrete(breaks=c(1,2,3,4,5,6),
                      labels=c("0-10", "10-20", "20-40", "40-65", "65-80", ">100")) +
xlab("Age Groups") +
ylab("The proportion of the isolated viruses with higher netcharge") +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))

## until here

       # geom_errorbar(limits, position="dodge", width=0.25) 



plotTop = 20;
barCenters <- barplot(height = means[,2],
            ylim = c(16, plotTop))


scale_size_area() +
xlab('Netcharge Changes') +
ylab('The effect of amino acid change on cell binding (log scale)') +
geom_boxplot() + 
geom_jitter(width=.5, height=0) 
geom_point()

plotTop = 20;
barCenters <- barplot(height = means[,2],
            ylim = c(16, plotTop))