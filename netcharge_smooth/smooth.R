x <- read.csv('titres.csv')
x <- read.csv('netcharge.csv')
X <- x[x[,1]>0,1]
Y <- x[x[,1]>0,2]
ct1 <- gam(Y~s(X),family=Gamma(link=log))
ct1 <- gam(Y~s(X),family=ocat(theta=NULL,link="identity",R=3.5))

plot(ct1,residuals=TRUE,pch=19,shade=TRUE)
#http://www3.nd.edu/~mclark19/learn/GAMS.pdf