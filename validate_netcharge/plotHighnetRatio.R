# Matlablimits

Ch <- read.csv('Chlist.dat')

#Ch1 <- Ch[which(Ch[,4]==10),]
Ch1 <- Ch
agegroup <- c(0,20,40,60,80,100)
agecode <- 0;
netthreshold = 19

for (a in 1:length(Ch1[,1])) {
  for (i in 1:5) {
    
    if (Ch1[a,1] >= agegroup[i] & Ch1[a,1] < agegroup[i+1]) {
      agecode[a] = i;
    }
  }
} 
Ch1[,10] <- agecode;
colnames(Ch1) <- c("age", "year", "netcharge", "X4", "X5", "X6", "X7", "X7", "X8","agegroup")

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
prev <- aid[Ch1[aid,3] >= netthreshold]
ratio[i] <- length(prev)/length(aid)
agesize[i] <- length(aid)
}



##-----^^^ finish above---------------


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
ylim(c(0, 35)) +
theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))

## the function ends until here


## the 2nd function starts from here
##plot sliding window
interval <- 8
meannet <- 0
highratio <- 0
netthreshold <- 19
id <- 1
for (i in 1:100) {
  netsubset <- subset(Ch1$netcharge , Ch1$age>i-interval & Ch1$age<i+interval)
  highratio[id] <- length(which(netsubset >= netthreshold))/length(netsubset)
  meannet[id] <- mean(netsubset)
  id <- id+1
}

highratiodf <- data.frame(1:100, highratio)
colnames(highratiodf) <- c("age", "highratio")

ggplot(highratiodf, aes(x = age, y = highratio)) +
  geom_point() +
  geom_smooth(method = "loess") +
  xlab("Age") +
  ylab("The Proportion of High Netcharge (%)") +
  theme(text = element_text(size=16))

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