pcaCV <- function (X, amax, center = TRUE, scale = TRUE, repl = 50, segments = 4, 
    segment.type = c("random", "consecutive", "interleaved"), 
    length.seg, trace = FALSE, plot.opt = TRUE, ...) 
{
    if (missing(amax)) {
        amax <- min(nrow(X) - 1, ncol(X))
    }
    X <- as.matrix(scale(X, center = center, scale = scale))
    amax <- min(amax, nrow(X) - max(sapply(segments, length)) - 1)
    optcomp <- matrix(NA, nrow = segments, ncol = repl)
    MSEP <- matrix(NA, nrow = repl, ncol = amax)
    dimnames(MSEP) <- list(paste("rep", 1:repl), paste("PC", 1:amax))
    Fit <- matrix(NA, nrow = repl, ncol = amax)
    dimnames(Fit) <- list(paste("rep", 1:repl), 1:amax)
    for (i in 1:repl) {
        if (missing(length.seg)) {
            segment <- cvsegments(nrow(X), k = segments, type = segment.type)
        }
        else {
            segment <- cvsegments(nrow(X), length.seg = length.seg, 
                type = segment.type)
        }
        MSEPall <- matrix(NA, nrow = nrow(X), ncol = amax)
        for (n.seg in 1:length(segment)) {
            seg <- segment[[n.seg]]
            obsuse <- as.numeric(unlist(segment[-n.seg]))
            Xtrain <- X[obsuse, ]
            obstest <- as.numeric(unlist(segment[n.seg]))
            Xtest <- X[obstest, ]
            if (ncol(Xtrain) > nrow(Xtrain)) {
                e <- eigen(Xtrain %*% t(Xtrain))
                Ttrain <- e$vectors %*% diag(sqrt(e$values))
                Ptrain <- t(Xtrain) %*% Ttrain %*% diag(1/e$values)
            }
            else {
                Xtrain_svd <- svd(Xtrain)
                Ptrain <- Xtrain_svd$v
            }
            Ttest <- Xtest %*% Ptrain
            for (j in 1:amax) {
                MSEPall[seg, j] <- apply((Xtest - Ttest[, 1:j] %*% t(Ptrain[, 1:j]))^2,1,sum)
            }
        }
        MSEP[i,] <- apply(MSEPall,2,sum)
        Fit[i,] <- 1-MSEP[i,]/sum(X^2)
    }
    if (plot.opt) {
        boxplot(as.data.frame(Fit), ylab = "Explained variance", 
            xlab = "Number of components", ...)
    }
    list(ExplVar = Fit, MSEP = MSEP)
}

