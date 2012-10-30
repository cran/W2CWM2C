WC <-
function(inputDATA, Wname, J, Hpdf, Wpdf) { 

  #: Check the input data. 
  if(is.ts(inputDATA) != "TRUE") 
  cat("The input data is not a time series, please check the ts 
  function in the R manual pages. Bye, thank you for your interest 
  in our program. \n") 

  #: Compute the dimensions 
  MNL <- dim(inputDATA) 
  ML  <- MNL[2] #No. columns
  NL  <- MNL[1] #No. elements
  if(ML >= NL) stop("Be careful with the input data, there 
   are more columns (variables) than number of elements.") 
 
  #:: To make the combinations (without repetition)    
  Labes <- seq(1:ML) 
  combcols <- combn(1:(ML), 2) 
  combSMI  <- combn(Labes,  2)
  Ncomb    <- ncol(combcols)  

  if(ML > 7) stop("This program only tackle arrays of N X 7 (columns) 
   dimensions, if you want to use array with more columns, please 
   use the Wavelet Multiple Correlation (function WMC). Bye, thank you 
   for your interest in our program. \n") 
 
  CEX.LAB = 1 
  if(ML == 7) CEX.LAB = 0.5
  if(ML == 6) CEX.LAB = 0.6 
  if(ML == 5) CEX.LAB = 0.7 

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
  
  wavcor.modwtsDAT <- array(0, c(Ncomb, J, 3)) 
 
  AXIS <- rep(" ", Ncomb)
  AXIS <- paste(paste("C", combcols[1,1:Ncomb], sep=""), 
          paste("C", combcols[2,1:Ncomb], sep=""), sep=" ")
  to.axix <- AXIS
   

  for (k in 1:Ncomb) {
   x   			 <- combcols[,k][1]
   y  			 <- combcols[,k][2]
   modwt.inputDATAx       <- modwt(inputDATA[,x], Wname, n.levels = J)
   modwt.inputDATAy       <- modwt(inputDATA[,y], Wname, n.levels = J)
   bw.modwinputDATAx      <- brick.wall(modwt.inputDATAx,  Wname)
   bw.modwinputDATAy      <- brick.wall(modwt.inputDATAy,  Wname)
   wavcormodwtsDAT       <- wave.correlation(bw.modwinputDATAx, 
                               bw.modwinputDATAy, N=NL)
   wavcor.modwtsDAT[k,,] <- as.matrix(wavcormodwtsDAT[-(J+1),])
  } 

  #:: Checking if the zero is inside of the CI  
  for (j in 1:J) { 
   jdx <- which(wavcor.modwtsDAT[,j,2] <= 0) 
   if (length(jdx) > 0) 
    wavcor.modwtsDAT[jdx,j,1] <- NaN
  }  

  to3Dp  <- wavcor.modwtsDAT[,,1] 
  xx     <- 1:Ncomb; yy <- 1:J

  j    <- 1:J 
  VECJ <- 2^(j-1) 
  VEC1 <- seq(1.5, Ncomb, 1) 
  VEC2 <- seq(1.5, J, 1)   
  pdf(file="./WCplot.pdf", width=Wpdf, height=Hpdf)
   layout(matrix(c(1,1), ncol=2, byrow=TRUE), widths=c(4,1))
   image(xx, yy, z=to3Dp, col=colors, xlab=" ",
     ylab="Wavelet Scale", xaxt="n", yaxt="n")
   to3Dp  <- round(to3Dp, digits=2)
   for(l in 1:J) { 
    text(seq(1,Ncomb), rep(l,Ncomb), to3Dp[,l], cex=0.75)
   }
   abline(v=c(VEC1))
   abline(h=c(VEC2))
   axis(1, at=1:Ncomb, labels=to.axix, cex.axis=CEX.LAB)
   axis(2, at=1:J, labels=VECJ)
  dev.off()

  to3DpL <- apply(t(to3Dp), 2, rev)
   
  return(list(wavcor.modwtsDAT=wavcor.modwtsDAT, to3DpL=to3DpL))
}

