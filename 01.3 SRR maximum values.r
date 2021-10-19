
library(ggplot2)

alpha=0.7;beta=2.98*10^(-5)
gen.srr <- data.frame(x=seq(from=0,to=80000,by=2))
gen.srr$y <- alpha*gen.srr$x*exp(-beta*gen.srr$x)

max(gen.srr$y);
gen.srr$x[which(gen.srr$y==max(gen.srr$y))]

ggplot(gen.srr) +
  geom_line(aes(x,y))


alpha=4.72;beta=1.06*10^(-4)
seg.srr.1 <- data.frame(x=seq(from=0,to=80000,by=2))
seg.srr.1$y <- alpha*seg.srr.1$x*exp(-beta*seg.srr.1$x)

max(seg.srr.1$y);
seg.srr.1$x[which(seg.srr.1$y==max(seg.srr.1$y))]

ggplot(seg.srr.1) +
  geom_line(aes(x,y))


alpha=0.61;beta=2.90*10^(-5)
seg.srr.2 <- data.frame(x=seq(from=0,to=80000,by=2))
seg.srr.2$y <- alpha*seg.srr.2$x*exp(-beta*seg.srr.2$x)

max(seg.srr.2$y);
seg.srr.1$x[which(seg.srr.2$y==max(seg.srr.2$y))]

ggplot(seg.srr.2) +
  geom_line(aes(x,y))


alpha=0.94;beta=1.31*10^(-4)
seg.srr.3 <- data.frame(x=seq(from=0,to=80000,by=2))
seg.srr.3$y <- alpha*seg.srr.3$x*exp(-beta*seg.srr.3$x)

max(seg.srr.3$y);
seg.srr.1$x[which(seg.srr.3$y==max(seg.srr.3$y))]

ggplot(seg.srr.3) +
  geom_line(aes(x,y))