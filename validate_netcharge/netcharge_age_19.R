# Matlablimits

Ch <- read.csv('Chlist.dat')

#Ch1 <- Ch[which(Ch[,4]==10),]
Ch1 <- Ch
agegroup <- c(0,20,40,60,80,100)
agecode <- 0;

for (a in 1:length(Ch1[,1])) {
  for (i in 1:5) {
    
    if (Ch1[a,1] >= agegroup[i] & Ch1[a,1] < agegroup[i+1]) {
      agecode[a] = i;
    }
  }
} 
Ch1[,10] <- agecode;
#can I reset column 10?
means <- aggregate(Ch1[,3] ~  Ch1[,10], Ch1, mean)
sds <- aggregate(Ch1[,3] ~  Ch1[,10], Ch1, sd)
se <- sds[,2]

#aid = which(Ch1[,10] %in% 1)
ratio <- 0;
agesize <- 0;
for (i in 1:length(agegroup)-1) {
#aid <- (Ch1[,10] == i)
aid = which(Ch1[,10] %in% i)
prev <- aid[Ch1[aid,3] >= 19]
ratio[i] <- length(prev)/length(aid)
agesize[i] <- length(aid)
}

##-----^^^ finish above---------------
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


## To plot Figure2. Net charge distribution by age group. 
## Please use the following code

df <- data.frame (
  age = factor(c(1,2,3,4,5)),
  highcharge = ratio,
  lb = c(0.1018,  0.0751,    0.0484,    0.1834,    0.1874),
  ub = c(0.1858,  0.1853,    0.1732,    0.3417,    0.3216)
)

limits <- aes(ymin=df$lb*100, ymax=df$ub*100)

ggplot(df, aes( x=age, y=highcharge*100)) +
geom_point(size = 5) +
geom_errorbar(limits, width=0.25) +
scale_x_discrete(breaks=c(1,2,3,4,5),
                      labels=c("0-20", "20-40", "40-60", "60-80", ">80")) +
xlab("Age Groups") +
ylab("The Proportion of High Netcharge (%)") +
ylim(c(0, 30)) +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))

## the function ends until here

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