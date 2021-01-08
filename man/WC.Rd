\encoding{latin1}
\name{WC}
\alias{WC}
\alias{wc}
\title{Wavelet correlation (bivariate case) pairwise comparisons. 
}
\description{
The \code{\link{WC}} function (bivariate case) computes the wavelet 
correlation by means of the function \emph{wave.-correlation} of
the \emph{waveslim} package to several time series, makes a 
pairwise comparisons and plot the pairwise wavelet correlations in 
descending order as a single heatmap using the \emph{colorspace} 
package. The input data are multivariate time series and \code{\link{WC}}
function only tackle arrays with N x C (elements x columns, where the
number of columns are between 2 and 7) dimensions. 
}
\usage{
WC(inputDATA, Wname, J, device="screen", filename, 
   Hfig, WFig, Hpdf, Wpdf)
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
The \code{\link{WC}} function compute the wavelet correlation among 
time series and plots the results in a single heatmap plot 
(which can be displayed in the screen or can be saved as 
PNG, JPG, PDF or EPS) showing the WC values as a table (please, 
look at Figure 3 in \strong{Polanco-Martinez} and Fernandez-Macho 2014). 
The \code{\link{WC}} code is based on the \emph{wave.correlation} routine 
from Brandon Whitcher's \emph{waveslim} \R package Version: 1.7.1,
which is based mainly on wavelet methodology developed by
Whitcher, B., P. Guttorp and D.B. Percival (2000) and Gencay, 
Selcuk and Whitcher (2001). 
}
\value{
Output:

Output plot: \emph{screen or 'filename + .png, .jpg, .eps or .pdf'}.

\emph{wavcor.modwtsDAT}: matrix with as many rows as levels in 
the wavelet transform object. The first column provides the 
point estimate for the wavelet correlation followed by the lower
and upper bounds from the confidence interval. \cr 

\emph{to3DpL}: A matrix (the matrix table added in the WC plot) 
with a J (number of wavelet scales) X C (the number of pairwise 
comparisons) dimensions, which are in descending order taking into 
account the sum of the wavelet correlation coefficients for all
(J) wavelet scales. 
}
\references{
Gencay, R., F. Selcuk and B. Whitcher (2001). \emph{An 
Introduction to Wavelets and Other Filtering Methods in Finance 
and Economics}, Academic Press.\cr 

Ihaka, R., Murrell, P., Hornik, K., Fisher, J. C. and Zeileis, A.
(2012). \emph{colorspace: Color Space Manipulation}. R package version
1.2.0, The Comprehensive R Archive Network (CRAN), <URL: \url{https://cran.r-project.org/package=colorspace}>. \cr

\strong{Polanco-Martinez, J.} and J. Fernandez-Macho (2014). The
package 'W2CWM2C': description, features and applications. 
\emph{Computing in Science & Engineering}, \strong{16}(6):68--78. \cr 
{\doi{10.1109/MCSE.2014.96}}. \cr

Whitcher, B., P. Guttorp, and D.B. Percival (2000). Wavelet analysis of
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, \strong{105}(D11):941--962. 
{\doi{10.1029/2000JD900110}}. \cr

Whitcher, B. (2012). \emph{waveslim: Basic wavelet routines for one-,
two- and three-dimensional signal processing}. R package version 1.7.1,
The Comprehensive R Archive Network (CRAN), \cr
<URL \url{https://cran.r-project.org/package=waveslim}>. \cr
}
\author{
Josue M. Polanco-Martinez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, Spain.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: josue.m.polanco@gmail.com.
}
\note{
Needs \emph{waveslim} package to calculate \emph{modwt}, 
 \emph{brick.wall} and the \emph{wave.correlation} and also 
needs the \emph{colorspace} package to plot the heatmaps.}
\examples{
## Figure 3 (Polanco-Martinez and Fernandez-Macho 2014). 

 library("colorspace")
 library("waveslim")
 library("W2CWM2C") 
 data(dataexample)   

 #:: Transforms to log returns using: ln(t + deltat) - ln(t). 
 #:: The application in this example uses stock market 
 #:: indexes (it is common to use log returns instead of
 #:: raw data). Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] # remove dates!
 dataexample  <- dataexample[,1:5]
 lrdatex      <- apply(log(dataexample), 2, diff)
 inputDATA    <- ts(lrdatex, start=1, frequency=1) 

 #Input parameters 
  Wname       <- "la8"
  J           <- 8 
  compWC      <- WC(inputDATA, Wname, J, device="screen", NULL,
                  NULL, NULL, NULL, NULL) 
}
\keyword{bi-variate}
\keyword{Pairwise comparisons}
\keyword{Wavelet Correlation}
