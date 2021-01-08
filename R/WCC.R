######################################################################
#:: WCC function: wavelet correlation -bivariate case from the       #
#:: R package W2CWM2C                                                #
#:: Programmed by Josue M. Polanco Martinez a.k.a jomopo             #
#:: josue.m.polanco@gmail.com                                        #
######################################################################
#:: Copyright (C) 2012, 2015, 2021 Josue M. Polanco Martinez 
#   This file is part of W2CWM2C 
#
#   W2CWM2C is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published 
#   by the Free Software Foundation, either version 3 of the License, 
#   or (at your option) any later version.
#
#   W2CWM2C is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with W2CWM2C. If not, see <http://www.gnu.org/licenses/>.
######################################################################

WCC <-
function(inputDATA, Wname, J, lmax, device="screen", filename, 
          Hfig, WFig, Hpdf, Wpdf) { 

  #: Check the input data. 
  if (is.ts(inputDATA) != "TRUE")
    stop("The input data is not a time series, please check the ts 
    function in the R manual pages. Bye, thank you for your interest 
    in our program. \n")

  namesSMI  <- colnames(inputDATA)
  cat("colnames", namesSMI) 
  NL        <- dim(inputDATA)[1]

  #:: Param.  to ourgraph wav. cros. cor.
  xx  <- seq(-lmax, lmax, 1)
  Nx  <- length(xx)
  yy  <- 1:J
  Ny  <- length(yy)

  modwt.inputDATAx  <- modwt(inputDATA[,1], Wname, n.levels = J)
  modwt.inputDATAy  <- modwt(inputDATA[,2], Wname, n.levels = J) 
  bw.modwinputDATAx <- brick.wall(modwt.inputDATAx, Wname)
  bw.modwinputDATAy <- brick.wall(modwt.inputDATAy, Wname)
 
  returns.cross.cor <- NULL
   for (i in 1:J) {
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
       setlinM[i,1]    <- which(returns.cross.cor[,i] == max(returns.cross.cor[,i]))
   }

  #:: Palette!
  Ncol     <- length(xx)*J*8
  Palette  <- diverge_hcl(Ncol, c=c(100,0), l=c(50,90), power=1.3)
  #:: colorbar! 
  rangev   <- seq(min(mat.zero, na.rm=TRUE), max(mat.zero, na.rm=TRUE), 
               length.out=J+1)
  rangebar <- matrix(rangev, nrow=1, ncol=J+1, byrow=TRUE)

  labv    <- xx[setlinM[,1]]
  X1      <- labv
  X2      <- labv
  p51     <- seq(0.5, J, 1) 
  p52     <- seq(1.5, J+1, 1) 
  VEClab  <- seq(-lmax, lmax, 5)
  VECJ    <- rep(0, J) 
  j       <- 1:J
  VECJ    <- 2^(j-1) 

  ## Devices options: png & jpg; esp & pdf! 
  if (device=="png") {
   fileout <- paste("WCC_", filename, ".png", sep="") 
   png(fileout, height=Hfig, width=WFig) 
  } 

  if (device=="jpeg" || device=="jpg") {
   fileout <- paste("WCC_", filename, ".jpg", sep="") 
   jpeg(fileout, height=Hfig, width=WFig) 
  } 

  if (device=="pdf") {
   fileout <- paste("WCC_", filename, ".pdf", sep="") 
   pdf(fileout, height=Hpdf, width=Wpdf)
  }

  if (device=="eps") {
   fileout <- paste("WCC_", filename, ".eps", sep="") 
   postscript(fileout, height=Hpdf, width=Wpdf)
  }

   layout(matrix(c(1,2), ncol=2, byrow=TRUE), widths=c(4,1))
   image(xx, yy, z=mat.zero, col=Palette, xlab="Lag (days)",
    ylab="Scale", main=paste(namesSMI[1], "vs.", namesSMI[2]), yaxt="n")
   segments(X1, p51, X2, p52, lty=2, lwd=2)
   abline(v=VEClab, lty=3, lwd=2)
   axis(2, at=1:J, labels=VECJ)
   image(z=rangebar, axes=FALSE, col=Palette, frame.plot=TRUE,
    yaxt="n", xaxt="n") 
   axis(2, at=round(seq(0,1,length.out=J+1),2), labels=round(rangebar, 
    digits=2), las=2)

  if (device != "screen") 
  dev.off()

  return(list(returns.cross.cor=returns.cross.cor))

  }
