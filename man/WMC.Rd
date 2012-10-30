\encoding{latin1}
\name{WMC}
\alias{WMC}
\alias{wmc}
\title{Wavelet multiple correlation (multivariate case).
}
\description{
The \emph{WMC} function only generates a pdf plot (\emph{WMC_plot.pdf}) 
to the Wavelet routine for multiple correlation 
(\emph{wave.multiple.correlation}) from the \emph{wavemulcor} package 
(Fernandez-Macho 2012), but also provides a way to handle multivariate
time series easily as a list of N elements (time series).  
}
\usage{
WMC(inputDATA, Wname, J, Hpdf, Wpdf)
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
  \item{Hpdf}{
The height of the pdf file (output plot).
}
  \item{Wpdf}{
The width of the pdf file (output plot).
}
}
\details{
The \emph{WMC} function helps to make easily the plot (in pdf 
format) of the multiple correlation routine 
(\emph{wave.multiple.correlation}) of the \emph{wavemulcor} package 
(Fernandez-Macho 2012). The \emph{WMC} function also helps 
to manage easily multivariate time series to use the Wavelet 
multiple correlation routine. 
}
\value{
Output:

Output pdf file: \emph{WMC_plot.pdf}. 

Output data: The same list of elements of the funtion \emph{
wave.multiple.correlation} of the \emph{wavemulcor} package 
(Fernandez-Macho 2012).
}
\references{
Fernandez-Macho, J.. Wavelet multiple correlation and
cross-correlation: A multiscale analysis of euro zone stock
markets. \emph{Physica A: Statistical Mechanics and its Applications}, 
391(4):1097-1104, 2012. \cr

Fernandez-Macho, J. (2012). Wavemulcor Reference manual. The 
Comprehensive R Archive Network (CRAN), 
http://cran.r-project.org/web/packages/wavemulcor/index.html.
}
\author{
Josue M. Polanco-Martinez (UPV/EHU).
Dpt. of Econometrics and Statistics,
& Instituto de Economia Publica,
University of the Basque Country, Agirre Lehendakari
etorb. 83, E48015 BILBAO, Spain.
Email: bcppomaj@ehu.es, jomopo@gmail.com.
}
\note{
Needs \emph{wavemulcor} (to compute the \emph{wave.multiple.correlation}) 
and \emph{waveslim} packages (to compute the \emph{modwt} and the 
\emph{brick.wall}).
}
\examples{

 library("wavemulcor")
 library("W2CWM2C")
 data(dataexample) 

 #:: Convert. log return using: ln(t + deltat) - ln(t) 
 #:: The application in this example is with stock market 
 #:: indexes and it is common to use log returns instead 
 #:: raw data. Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] #remove the dates!
 dataexample  <- dataexample[,1:5]
 lrdatex      <- apply(log(dataexample), 2, diff)
 tslrdat      <- ts(lrdatex, start=1, frequency=1)

 #Input parameters 
  Wname       <- "la8"
  J           <- 8
  Hp          <- 6
  Wp          <- 10
  tslrdat     <- tslrdat[,1:5]
  compWMC     <- WMC(tslrdat, Wname, J, Hp, Wp)
}
\keyword{ Multi-variate}
\keyword{Wavelet multiple correlation}
