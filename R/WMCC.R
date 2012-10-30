WMCC <-
function(inputDATA, Wname, J, lmax, Hpdf, Wpdf) {

  if(is.ts(inputDATA) != "TRUE")
   stop("The input data is not a time series, please check the ts 
    function in the R manual pages. Bye, thank you for your interest 
    in our program. \n")

  NAMES <- colnames(inputDATA)

  #: Compute dimensions 
  MNL <- dim(inputDATA)
  M  <- MNL[2] #No. columns
  N  <- MNL[1] #No. elements
  if(M >= N) stop("Be careful with the input data, there 
   are more columns (variables) than number of elements.")

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
 
  #: Comp. the MODWT 
  LIST <- as.list("NULL")
  for (l in 1:M) {
    tmp.wt  <- modwt(inputDATA[,l], Wname, J)
    tmp.bw  <- brick.wall(tmp.wt, Wname)
    LIST[l] <- list(tmp.bw)
  }

  LS                <- wave.multiple.cross.correlation(LIST, lmax)
  returns.cross.cor <- as.matrix(LS$xy.mulcor[1:J,])
  YmaxR             <- LS$YmaxR
  lags              <- length(-lmax:lmax)

  lower.ci <- tanh(atanh(returns.cross.cor) - qnorm(0.975) /
  sqrt(matrix(trunc(N/2^(1:J)), nrow=J, ncol=lags)- 3))
  upper.ci <- tanh(atanh(returns.cross.cor) + qnorm(0.975) /
  sqrt(matrix(trunc(N/2^(1:J)), nrow=J, ncol=lags)- 3))

  mat.zero <- array(NaN, c(Nx, Ny))
  setlinM  <- array(NaN, c(J, 1))
  for (i in 1:J) {
   idx <- which ( upper.ci[i,] < 0 | lower.ci[i,] > 0 )
   mat.zero[idx, i] <- returns.cross.cor[i, idx]
   setlinM[i,1] <- which(returns.cross.cor[i,] == max(returns.cross.cor[i,]))
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
  pdf(file="./WMCC_plot.pdf", width = Wpdf, height = Hpdf)
  layout(matrix(c(1,2), ncol=2, byrow=TRUE), widths=c(4,1))
  image(xx, yy, z=mat.zero, col=colors, xlab="Lag (days)",
   ylab="Scale", main=paste("WMCC"), yaxt="n")
  text(rep(lmax-7,J),1:J,labels=NAMES[YmaxR[1:J]],adj=0.25,cex=1.2,col="black")
  segments(X1, p51, X2, p52, lty=2, lwd=2)
  abline(v=VEClab, lty=3, lwd=2)
  axis(2, at=1:J, labels=VECJ)
  image(z=rangebar, axes=FALSE, col=colors, frame.plot=TRUE,
   yaxt="n", xaxt="n") 
  axis(2, at=Lrangcol, labels=round(Lrangbar, digits=2), las=2)
  dev.off()
 
 return(list(LS=LS))
 
}

