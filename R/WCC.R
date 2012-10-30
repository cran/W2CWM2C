WCC <-
function(inputDATA, Wname, J, lmax, Hpdf, Wpdf) { 

  if(is.ts(inputDATA) != "TRUE")
  stop("The input data is not a time series, please check the ts 
  function in the R manual pages. Bye, thank you for your interest 
  in our program. \n")

  namesSMI  <- colnames(inputDATA)
  cat("colnames", namesSMI) 
  NL        <- dim(inputDATA)[1]

##Function from SOWAS (Maraun & Kurths 2004)###
 wtcolors <- function(ncolors){
     upside <- rainbow(ncolors,start=0,end=.7)
     down   <- upside
     for (i in 1:ncolors) {
      down[i] <- upside[ncolors+1-i]
     }
    down
   }
  ncolors  <- 128 
  colors   <- wtcolors(ncolors)

  #:: Param. para ourgraph wav. cros. cor.
  xx  <- seq(-lmax, lmax, 1)
  Nx  <- length(xx)
  yy  <- 1:J
  Ny  <- length(yy)

  #:: The auxiliary barcolour function  
  rangecol <- (0:(ncolors-1)) / (ncolors-1)
  byseq    <- floor(ncolors/10)
  atlab    <- seq(2, length(rangecol), by=byseq)
  Lrangcol <- rangecol[atlab]
 
   modwt.inputDATAx  <- modwt(inputDATA[,1], Wname, n.levels = J)
   modwt.inputDATAy  <- modwt(inputDATA[,2], Wname, n.levels = J) 
   bw.modwinputDATAx <- brick.wall(modwt.inputDATAx, Wname)
   bw.modwinputDATAy <- brick.wall(modwt.inputDATAy, Wname)
 
 returns.cross.cor <- NULL
   for(i in 1:J) {
    wcc <- spin.correlation(bw.modwinputDATAx[[i]],
             bw.modwinputDATAy[[i]], lmax)
    returns.cross.cor <- cbind(returns.cross.cor, wcc)
   }

   returns.cross.cor           <- ts(as.matrix(returns.cross.cor),
                                      start=-lmax, frequency=1)
   dimnames(returns.cross.cor) <- list(NULL, paste("Level", 1:J))
   lags                        <- length(-lmax:lmax)
   lower.ci <- tanh(atanh(returns.cross.cor) - qnorm(0.975) /
                sqrt(matrix(trunc(NL/2^(1:J)), nrow=lags, ncol=J,
                byrow=TRUE) - 3))
   upper.ci <- tanh(atanh(returns.cross.cor) + qnorm(0.975) /
                sqrt(matrix(trunc(NL/2^(1:J)), nrow=lags, ncol=J,
                byrow=TRUE) - 3))
 

   #: plot the WCC 
   mat.zero <- array(NaN, c(Nx, Ny))
   setlinM  <- array(NaN, c(J, 1))
   for (i in 1:J) {
    idx <- which ( upper.ci[,i] < 0 | lower.ci[,i] > 0 )
    mat.zero[idx,i] <- returns.cross.cor[idx,i]
    setlinM[i,1] <- which(returns.cross.cor[,i] == max(returns.cross.cor[,i]))
   }
   labv     <- xx[setlinM[,1]]
   X1       <- labv
   X2       <- labv
   p51      <- seq(0.5, J, 1) 
   p52      <- seq(1.5, J+1, 1) 
   rangev   <- seq(min(returns.cross.cor), max(returns.cross.cor),
                   length.out=ncolors)
   rangebar <- matrix(rangev,nrow=1,ncol=ncolors,byrow=TRUE)
   Lrangbar <- rangebar[atlab]

   VEClab <- seq(-lmax, lmax, 5)
   VECJ   <- rep(0, J) 
   j <- 1:J
   VECJ <- 2^(j-1) 
   pdf(file=paste("./wcc_", namesSMI[1],  namesSMI[2], ".pdf", 
    sep=""), width = Wpdf, height = Hpdf)
    layout(matrix(c(1,2), ncol=2, byrow=TRUE), widths=c(4,1))
    image(xx, yy, z=mat.zero, col=colors, xlab="Lag (days)",
     ylab="Scale", main=paste(namesSMI[1], "vs.", namesSMI[2]), yaxt="n")
    segments(X1, p51, X2, p52, lty=2, lwd=2)
    abline(v=VEClab, lty=3, lwd=2)
    axis(2, at=1:J, labels=VECJ)
    image(z=rangebar, axes=FALSE, col=colors, frame.plot=TRUE,
     yaxt="n", xaxt="n") 
    axis(2, at=Lrangcol, labels=round(Lrangbar, digits=2), las=2)
   dev.off()

  return(list(returns.cross.cor=returns.cross.cor))

  }

