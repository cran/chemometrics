prm_dcv   <- function(  X,                            # matrix for X
                        Y,                            # vector for y
                        ncomp=10,                     # max number of PLS components
                        repl=10,                       # repetitions in rdcv
                        segments0=4,                  # number of outer segments
                        segments=7,                   # number of inner segments
                        segment0.type="random",       # type of outer segmentation 
                        segment.type="random",        # type of inner segmentation
                        sdfact=2,                     # parsimony factor for a_opt
                        fairct=4,                     # factor for robust pls weights
                        trim=0.2,                     # trim factor for robust pls
                        opt="median",                 # calculation method for multivariate median
                        plot.opt=FALSE, ...)
{
require(chemometrics)
n = nrow(X)                                             # number of samples

optcomp <- matrix(NA, nrow=segments0, ncol=repl)      # a_opt cv [1:segments0, 1:repl] 
b       <- matrix(NA, nrow=dim(X)[2], ncol=ncomp)
pred    <- array(NA, dim=c(n, ncomp, repl))            # y_test [1:n, 1:ncomp, 1:repl]
predopt <- matrix(NA, nrow=n, ncol=repl)              # y_test for each a_opt [1:n, 1:repl] 


  # (1) --- start repetition loop --- #
  for(i in 1:repl)                                         
  {
    cat("\n", "repl: ", i, "\n")
    segment0 <- cvsegments(n, k=segments0, type=segment0.type)      # create outer segments
  
            
            # (2) --- start outer loop --- #
            for(n.seg0 in 1:length(segment0))                       
            {
            cat("\n", "seg-nr: ", n.seg0, "\n")                     
            seg0 <- segment0[[n.seg0]]                              # test set samples
            obsuse <- as.numeric(unlist(segment0[-n.seg0]))         # calibration set samples
            d1 <- list(X=X[obsuse,],Y=Y[obsuse])                    # data for calibration
  
                # (3) --- start inner loop --- #
                  res <-  prm_cv(d1$X, d1$Y,                   # robust PLS with (inner segments)-fold CV with calibration data
                          a=ncomp,fairct=fairct, 
                          opt=opt, segments=segments,
                          segment.type=segment.type, trim=trim,
                          sdfact=sdfact, plot.opt=plot.opt)
                  
                  optcomp[n.seg0,i] <- res$optcomp                  # a_opt cv for current calibration set
                # (3) --- end inner loop --- #
          
             
             # Extract robust PLS models' regression coefficients 
             b0 <- vector(length=ncomp)                              # robust intercept
             for(n.comp in 1:ncomp)                                  # for all numbers of PLS components 
             {
                rcal <- prm(d1$X, d1$Y, a=n.comp, fairct=fairct, opt=opt)   # robust PLS with entire current calibration set
    
                b[,n.comp]  <- rcal$coef                             # robust regr.coeff
                b0[n.comp] <- rcal$intercept                         # robust intercept
              }
                b0 <- matrix(rep(b0, nrow(X[-obsuse,])), ncol=ncomp, byrow=TRUE) 
    
             # test set predicted y  
             pred[-obsuse,,i]  <-  as.matrix(scale(X[-obsuse,], center=rcal$mx, scale=FALSE ))%*% b + b0                #+ ym  aus prm
             predopt[-obsuse,i] <- drop(pred[-obsuse,optcomp[n.seg0,i],i])
             }
             # (2) --- end outer loop --- #
  }
  # (1) --- end repetition loop --- #


# a_final as most frequent a_opt cv
afinaldistr <- table(optcomp)/sum(table(optcomp))
afinal <- as.numeric(names(which.max(afinaldistr)))

# SEPopt considering only residuals from a_opt cv models
resopt <- predopt - c(Y)

biasopt <- apply(resopt, 2, mean)
SEPopt <- sqrt(apply((resopt-biasopt)^2,2,sum)/(prod(dim(resopt)) - 1))

# SEPfinal considering all residuals for each number of PLS components
residcomp <- pred - c(Y)

biascomp <- apply(residcomp, 2, mean)
dimr <- dim(residcomp)
biascomp1 <- array(biascomp, c(dimr[2], dimr[1], dimr[3]))
biascomp1 <- aperm(biascomp1, c(2,1,3))
#SEPfinal <- sqrt(apply((residcomp-biascomp)^2,2,sum)/(prod(dim(resopt)) - 1))
SEPfinal <- sqrt(apply((residcomp - biascomp1)^2,2,sum)/(prod(dim(resopt)) - 1))  
# SEPtrim considering residuals at a_final only, excluding trim-% (e.g. 20%) of largest residuals
medres <- median(residcomp[, afinal,])
residcomp1 <- residcomp[, afinal,] - medres
absres <- abs(residcomp1)
reslarge <- quantile(absres, 1-trim)
SEPtrim <- sd(residcomp1[absres <= reslarge])

list(b=b, resopt=resopt, predopt=predopt, optcomp=optcomp, residcomp=residcomp, 
     pred=pred, SEPopt=SEPopt,afinal=afinal, SEPfinal=SEPfinal, SEPtrim=SEPtrim)
}
 
 
