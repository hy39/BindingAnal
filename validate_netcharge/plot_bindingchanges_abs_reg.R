# plot the correlation between binding avidity and absolute netcharge
# to demonstrate not only netcharge change but also the absolute netcharge will affect binding

M <- read.csv('RDE/abs_charge.csv')
plot(M$abs_netcharge, M$abs_binding, xaxt = "n", xlab = "Net charge",
       ylab = "Relative binding avidity (log scale)")

axis(1, at=-2:2, labels=16:20)

reg1 <- lm(M$abs_netcharge~M$abs_binding)
abline(reg1) +
  ylim(c(0, 80)) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))
summary(reg1)


$$ codes of working out t test
high <- M$abs_binding[M$abs_netcharge>=1]
low <- M$abs_binding[M$abs_netcharge<1]
t.test(high, low)    