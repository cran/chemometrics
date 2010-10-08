# alr transformation with variable 2:
source("C:/DA/Literatur, pdfs/compositional data/outlier detection/drawMahal.R")
source("C:/DA/Literatur, pdfs/compositional data/outlier detection/alr.R")
source("C:/DA/Literatur, pdfs/compositional data/outlier detection/invalr.R")
library(compositions)

# normal distributed data in real space
dat <- mvrnorm(n=40, mu=c(1,2), Sigma=matrix(c(0.4225,0.6591,0.6591,1.69), nrow=2), empirical=TRUE)
logxy <- dat[,1]
logzy <- dat[,2]
# alr-backtransformed data (in simplex)
x <- exp(logxy) / (exp(logxy)+exp(logzy)+1)
z <- exp(logzy) / (exp(logxy)+exp(logzy)+1)
y <-          1 / (exp(logxy)+exp(logzy)+1)
compdata <- as.data.frame(cbind(x, y, z))
d=compdata

# transform:
d3=alr2(d[,1],d[,2],d[,3])
da=invalr2(d3[,1],d3[,2]) #row sum = 1

pdf(file="C:/DA/_Diplomarbeit_/Figures/Fig-02-01-logratio.pdf", width=9, height=4.5)
par(mfrow=c(1,2))

par(mar=c(3,2,2,2))
plot.acomp(da, pch=20)
for (j in 1:ncol(mahal$mdX)){
  plot.acomp(invalr2(mahal$mdX[,j],mahal$mdY[,j]),type="l",add=TRUE)#,cn=c("A","F","M"))
}
title("original data in restricted space")

par(mar=c(5,4,3,2))
plot(d3,xlab="x1",ylab="x2",pch=20,xlim=c(-1,3),ylim=c(-2,6))
mahal = drawMahal(d3,apply(d3,2,mean),cov(d3))
title("alr transformed data in unrestricted space")
dev.off()

