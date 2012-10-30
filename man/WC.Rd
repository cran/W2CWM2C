\encoding{latin1}
\name{WC}
\alias{WC}
\alias{wc}
\title{Wavelet correlation (bivariate case) pairwise comparisons. 
}
\description{
The \emph{WC} function (bivariate case) computes the Wavelet 
correlation by means of the function \emph{wave.correlation} of
the \emph{waveslim} package to several time series and make a 
pairwise comparisons and plot the wavelet correlation as a single 
pdf file plot. The input data are multivariate time series and this 
function only tackle arrays with N x C (elements x columns, where the
number of columns are between 2 and 7) dimensions. 
}
\usage{
WC(inputDATA, Wname, J, Hpdf, Wpdf)
}
\arguments{
  \item{inputDATA}{
An array of multivariate time series as a \emph{ts} object 
(please, check the \emph{ts} manual to get more information 
about the \emph{ts} function in R). 
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
The \emph{WC} function compute the wavelet correlation among 
time series and plots the results in a single pdf file 
(\emph{WCplot.pdf}) showing the WC values as a table (please, 
see the Figure 1 of Polanco-Martinez and Fernandez-Macho 2012). The 
\emph{WC} code is based on the \emph{wave.correlation} routine 
from Brandon Whitcher's \emph{waveslim} \R package Version: 1.7.1,
which is based mainly on wavelet methodology developed in 
Whitcher, B., P. Guttorp and D.B. Percival (2000) and Gencay, 
Selcuk and Whitcher (2001). 
}
\value{
Output:

Output pdf file: \emph{WCplot.pdf}.

\emph{wavcor.modwtsDAT}: matrix with as many rows as levels in 
the wavelet transform object. The first column provides the 
point estimate for the wavelet correlation followed by the lower
and upper bounds from the confidence interval. \cr 

\emph{to3DpL}: A matrix (the matrix table added in the WCplot.pdf 
plot) with a J (number of wavelet scales) X C (the number of pairwise 
comparisons) dimensions. 
}
\references{
Gencay, R., F. Selcuk and B. Whitcher (2001). \emph{An 
Introduction to Wavelets and Other Filtering Methods in Finance 
and Economics}, Academic Press.\cr 

Polanco-Martinez, J. and J. Fernandez-Macho (2012). The
package 'W2CWM2C': description, features and applications.
To be submitted under review to \emph{Journal of 
Statistical Software}.\cr

Whitcher, B., P. Guttorp, and D.B. Percival. Wavelet analysis of
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, 105(D11):941-962, 2000. \cr

Whitcher, B. (2012). Waveslim reference manual. The Comprehensive 
R Archive Network (CRAN), \cr
http://cran.r-project.org/web/packages/waveslim/index.html \cr
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
\emph{brick.wall} and the \emph{wave.correlation}. 
}
\examples{
## See the Figure 1 of Polanco-Martinez and Fernandez-Macho 
## 2012 (to be submitted under review to Journal of 
## Statistical Software). 

 library("waveslim")
 library("W2CWM2C") 
 data(dataexample)   

 #:: Convert. log returns using: ln(t + deltat) - ln(t) 
 #:: The application in this example is with stock market 
 #:: indexes and it is common to use log returns instead 
 #:: raw data. Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] # remove dates!
 dataexample  <- dataexample[,1:5]
 lrdatex      <- apply(log(dataexample), 2, diff)
 tslrdat      <- ts(lrdatex, start=1, frequency=1) 

 #Input parameters 
  Wname       <- "la8"
  J           <- 8 
  Hp          <- 6
  Wp          <- 10
  tslrdat     <- tslrdat[,1:5]
  compWC      <- WC(tslrdat, Wname, J, Hp, Wp)
}
\keyword{bi-variate}
\keyword{Pairwise comparisons}
\keyword{Wavelet Correlation}
