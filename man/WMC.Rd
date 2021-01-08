\encoding{latin1}
\name{WMC}
\alias{WMC}
\alias{wmc}
\title{Wavelet multiple correlation (multivariate case).
}
\description{
The \code{\link{WMC}} function generates a plot to the wavelet 
routine for multiple correlation \cr (\emph{wave.multiple.correlation})
from the \emph{wavemulcor} package (Fernandez-Macho 2012b). 
The \code{\link{WMC}} plot output can be displayed in the screen 
(by default) or can be saved as PNG, JPG, PDF or EPS. Furthermore, 
it also provides a way to handle multivariate time series easily as 
a list of N elements (time series).  
}
\usage{
WMC(inputDATA, Wname, J, device="screen", filename,
          Hfig, WFig, Hpdf, Wpdf)
}
\arguments{
  \item{inputDATA}{
A couple of time series as a \emph{ts} object (please, check the 
\emph{ts} manual to get more information about the \emph{ts} 
function in R).
}
  \item{Wname}{
The wavelet function or filter to use in the decomposition.
}
  \item{J}{
Specifies the depth of the decomposition.
}
  \item{device}{
The type of the output device (by default the option is ``screen'',
 and the other options are ``jpg'', ``png'', ``eps'' and ``pdf'').
}
  \item{filename}{
The output filename.
}
  \item{Hfig}{
The height of the 'jpg' or 'png' image. 

}
  \item{WFig}{
The width of the 'jpg' or 'png' image. 
}
  \item{Hpdf}{
The height of the 'eps' or 'pdf'.
}
  \item{Wpdf}{
The width of the 'eps' or 'pdf'. 
}
}
\details{
The \code{\link{WMC}} function helps to make and save easily the plot of the 
multiple correlation routine \cr (\emph{wave.multiple.correlation}) of the 
\emph{wavemulcor} package (Fernandez-Macho 2012b). The \code{\link{WMC}}
function also helps to manage easily multivariate time series to use 
the Wavelet multiple correlation routine. 
}
\value{
Output:

Output plot: \emph{screen or 'filename + .png, .jpg, .eps or .pdf'}.

Output data: The same list of elements of the funtion \emph{
wave.multiple.correlation} of the \emph{wavemulcor} package 
(Fernandez-Macho 2012b).
}
\references{
Fernandez-Macho, J. (2012a). Wavelet multiple correlation and
cross-correlation: A multiscale analysis of Euro zone stock
markets. \emph{Physica A: Statistical Mechanics and its Applications}, \cr
\strong{391}(4):1097--1104. {\doi{10.1016/j.physa.2011.11.002}}. \cr

Fernandez-Macho, J. (2012b). \emph{wavemulcor: Wavelet routine for
multiple correlation}. R package version 1.2, The Comprehensive R
Archive Network (CRAN), <URL: \url{https://cran.r-project.org/package=wavemulcor}>. \cr

\strong{Polanco-Martinez, J.} and J. Fernandez-Macho (2014). The
package 'W2CWM2C': description, features and applications. 
\emph{Computing in Science & Engineering}, \strong{16}(6):68--78. \cr 
{\doi{10.1109/MCSE.2014.96}}. \cr
}
\author{
Josue M. Polanco-Martinez (a.k.a. jomopo)..\cr
BC3 - Basque Centre for Climate Change, Bilbao, Spain.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: josue.m.polanco@gmail.com.
}
\note{
Needs \emph{wavemulcor} (to compute the \emph{wave.multiple.correlation}) 
and \emph{waveslim} packages (to compute the \emph{modwt} and the 
\emph{brick.wall}).
}
\examples{

 # This example is the wavelet multiple correlation (WMC) version of 
 # the Figure 7 in Polanco-Martinez and Fernandez-Macho (2014).

 library("wavemulcor")
 library("W2CWM2C")
 data(dataexample) 

 #:: Transform to log returns using: ln(t + deltat) - ln(t). 
 #:: The application in this example uses stock market 
 #:: indexes (it is common to use log returns instead of 
 #:: raw data). Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] #remove the dates!
 dataexample  <- dataexample[,1:5]
 lrdatex      <- apply(log(dataexample), 2, diff)
 inputDATA    <- ts(lrdatex, start=1, frequency=1)

 #Input parameters 
  Wname       <- "la8"
  J           <- 8
  compWMC     <- WMC(inputDATA, Wname, J, device="screen", NULL,
                     NULL, NULL, NULL, NULL)
}
\keyword{ Multi-variate}
\keyword{Wavelet multiple correlation}
