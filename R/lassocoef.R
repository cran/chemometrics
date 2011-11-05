lassocoef <-
function(formula,data,sopt,plot.opt=TRUE, ...)
{
# Plot coefficients of Lasso regression
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
    X <- pls:::delete.intercept(model.matrix(mt, mf))


mod_lasso <- lars(X,y)

aa <- apply(abs(mod_lasso$beta),1,sum)
ind <- which.min(abs(aa/max(aa)-sopt))
coef <- mod_lasso$beta[ind[1],]
numb.zero <- sum(coef==0)
numb.nonzero <- sum(coef!=0)

if (plot.opt){
plot(mod_lasso,breaks=FALSE,cex=0.4,col=gray(0.6),...)
abline(v=sopt,lty=2)
title(paste(numb.zero,"coefficients are zero,",numb.nonzero,"are not zero"))
}

list(coefficients=coef,sopt=sopt,numb.zero=numb.zero,numb.nonzero=numb.nonzero,ind=ind)
}

