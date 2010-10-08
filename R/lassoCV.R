lassoCV <-
function(formula,data,K = 10,fraction=seq(0,1,by=0.05),trace=FALSE,plot.opt=TRUE,
	sdfact=2, legpos = "topright", ...)
{
# Cross Validation for Lasso regression
#

require(pls)
require(lars)

    mf <<- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    X <- delete.intercept(model.matrix(mt, mf))

all.folds <- cv.folds(length(y), K)
    residmat <- matrix(0, length(fraction), K)
    SEPmat <- matrix(0, length(fraction), K)
    for (i in seq(K)) {
        omit <- all.folds[[i]]
        fit <- lars(X[-omit, , drop = FALSE], y[-omit], trace = trace,
            ...)
        fit <- predict(fit, X[omit, , drop = FALSE], mode = "fraction",
            s = fraction)$fit
        if (length(omit) == 1)
            fit <- matrix(fit, nrow = 1)
        residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
        SEPmat[, i] <- apply(y[omit] - fit, 2, sd)
        if (trace)
            cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean)
    SEP <- apply(SEPmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)

if (plot.opt){

error.bars <-
function (x, upper, lower, width = 0.02, ...)
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
#    range(upper, lower)
}

ind <- which.min(cv)
sopt<-fraction[ind]
plot(fraction,cv,ylim=range(cv-sdfact*cv.error,
 cv+sdfact*cv.error),xlab="|beta|/max|beta|",ylab="MSEP",cex.lab=1.2,
  type="n")
lines(fraction,cv)
points(fraction,cv,pch=16,cex=0.5)
error.bars(fraction,cv+sdfact*cv.error,
  cv-sdfact*cv.error,width = 1/length(fraction))
abline(h=cv[ind]+sdfact*cv.error[ind],lty=4)
abline(v=sopt,lty=4)
legend(legpos,c(paste("MSEP=",round(cv[ind],2),sep=""),
  paste("SEP=",round(SEP[ind],2),sep="")))
}


list(cv=cv,cv.error=cv.error,SEP=SEP,ind=ind,sopt=sopt,fraction=fraction)
}

