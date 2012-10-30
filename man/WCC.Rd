\encoding{latin1}
\name{WCC}
\alias{WCC}
\alias{wcc}
\title{Wavelet cross-correlation (bivariate case). 
}
\description{
The \emph{WCC} function (bivariate case) computes the wavelet
cross correlation using the \emph{spin.correlation} function 
of the \emph{waveslim} package for two time series and presents
the result as a novel plot that reduce the number of plots of 
the classical function \emph{spin.correlation}.  
}
\usage{
WCC(inputDATA, Wname, J, lmax, Hpdf, Wpdf)
}
\arguments{
  \item{inputDATA}{
A couple of time series as a \emph{ts} object (please, check 
the \emph{ts} manual to get more information about the \emph{
ts} function in R).
}
  \item{Wname}{
The wavelet function or filter to use in the decomposition.
}
  \item{J}{
Specifies the depth of the decomposition.
}
  \item{lmax}{
The maximum lag.
}
  \item{Hpdf}{
The height of the pdf file (output plot).
}
  \item{Wpdf}{
The width of the pdf file (output plot).
}
}
\details{
The \emph{WCC} function compute the Wavelet cross-correlation
between two time series and plot the results in a single pdf 
file (\emph{'wcc + var1 + var2.pdf'}) (please, see the Figure 3 of 
Polanco-Martinez and Fernandez-Macho 2012). The \emph{WCC} code 
is based on the \emph{spin.correlation} routine from Brandon 
Whitcher's \emph{waveslim} \R package Version: 1.7.1, which is 
based mainly on wavelet methodology developed in Whitcher, B., P. 
Guttorp and D.B. Percival (2000) and Gencay, Selcuk and Whitcher 
(2001).
}
\value{
Output:

Output pdf file: \emph{wcc_ + var1 + var 2 .pdf}

\emph{returns.cross.cor}: a matrix with the WCC values. 
}
\references{
Gencay, R., F. Selcuk and B. Whitcher (2001). \emph{An 
Introduction to Wavelets and Other Filtering Methods in 
Finance and Economics}, Academic Press.\cr

Polanco-Martinez, J. and J. Fernandez-Macho (2012). The
package 'W2CWM2C': description, features and applications.
To be submitted under review to \emph{Journal of 
Statistical Software}.\cr

Whitcher, B., P. Guttorp, and D.B. Percival. Wavelet analysis of
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, 105(D11):941-962, 2000. \cr

Whitcher, B. (2012). Waveslim reference manual. The Comprehensive 
R Archive Network (CRAN), \cr
http://cran.r-project.org/web/packages/waveslim/index.html\cr
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
Needs \emph{waveslim} package to calculate \emph{modwt}, 
\emph{brick.wall} and \emph{spin.correlation}. 
}
\examples{
## See the Figure 3 of Polanco-Martinez and Fernandez-Macho 
## 2012 (to be submitted under review to Journal of 
## Statistical Software). 

 library("waveslim")
 library("W2CWM2C")
 data(dataexample)  

 #:: Convert. log return using: ln(t + deltat) - ln(t) 
 #:: The application in this example is with stock market 
 #:: indexes and it is common to use log returns instead 
 #:: raw data. Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] #remove the dates!
 DAXCAC       <- dataexample[,c(3,4)] 
 lrdatex      <- apply(log(DAXCAC), 2, diff)
 tslrdat      <- ts(lrdatex, start=1, frequency=1)

 Wname   <- "la8"
 J       <- 8
 Hp      <- 6
 Wp      <- 10
 lmax    <- 30
 compWCC <- WCC(tslrdat, Wname, J, lmax, Hp, Wp)
}
\keyword{bi-variate}
\keyword{Wavelet cross correlation }
