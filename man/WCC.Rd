\encoding{latin1}
\name{WCC}
\alias{WCC}
\alias{wcc}
\title{Wavelet cross-correlation (bivariate case). 
}
\description{
The \code{\link{WCC}} function (bivariate case) computes the wavelet
cross correlation using the \emph{spin.correlation} function 
of the \emph{waveslim} package for two time series, and presents
the result as a plot that reduce the number of plots of the 
classical function \emph{spin.correlation}. The heatmap plot is 
built using the \emph{colorspace} package and can be displayed
in the screen or can be saved as PNG, JPG, PDF or EPS.
}
\usage{
WCC(inputDATA, Wname, J, lmax, device="screen", filename,
    Hfig, WFig, Hpdf, Wpdf)
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
The \code{\link{WCC}} function compute the wavelet cross-correlation
between two time series and plot the results in a single heatmap 
plot (please, look at Figure 3 in Polanco-Martinez and Fernandez-Macho 
2014). The \emph{WCC} code is based on the \emph{spin.correlation} 
routine from Brandon Whitcher's \emph{waveslim} \R package Version: 
1.7.1, which is based mainly on wavelet methodology developed in 
Whitcher, B., P. Guttorp and D.B. Percival (2000) and Gencay, Selcuk 
and Whitcher (2001).
}
\value{
Output:

Output plot: \emph{screen or 'filename + .png, .jpg, .eps or .pdf'}. 

\emph{returns.cross.cor}: a matrix with the WCC values. 
}
\references{
Gencay, R., F. Selcuk and B. Whitcher (2001). \emph{An 
Introduction to Wavelets and Other Filtering Methods in 
Finance and Economics}, Academic Press.\cr

Ihaka, R., Murrell, P., Hornik, K., Fisher, J. C. and Zeileis, A.
(2012). \emph{colorspace: Color Space Manipulation}. R package version
1.2.0. The Comprehensive R Archive Network (CRAN), \cr
\url{http://cran.r-project.org/package=colorspace} \cr

Polanco-Martinez, J. and J. Fernandez-Macho (2014). The
package 'W2CWM2C': description, features and applications.
\emph{Computing in Science & Engineering},  16(6):68--78, 
\href{http://www.computer.org/csdl/mags/cs/preprint/06894486.pdf}{DOI: 10.1109/MCSE.2014.96}.\cr

Whitcher, B., P. Guttorp, and D.B. Percival (2000). Wavelet analysis of
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, 105(D11):941--962. \cr

Whitcher, B. (2012). \emph{waveslim: Basic wavelet routines for one-,
two- and three-dimensional signal processing}. R package version 1.7.1.
The Comprehensive R Archive Network (CRAN), \cr
\url{http://cran.r-project.org/package=waveslim} \cr
}
\author{
Josue M. Polanco-Martinez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, Spain.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: josue.m.polanco@gmail.com.
}
\note{
Needs \emph{waveslim} package to calculate \emph{modwt}, 
\emph{brick.wall} and \emph{spin.correlation} and also needs 
the \emph{colorspace} package to plot the heatmaps. 
}
\examples{
## Figure 3 (Polanco-Martinez and Fernandez-Macho 2014)

 library("colorspace")
 library("waveslim")
 library("W2CWM2C")
 data(dataexample)  

 #:: Convert. log return using: ln(t + deltat) - ln(t) 
 #:: The application in this example is with stock market 
 #:: indexes and it is common to use log returns instead of
 #:: raw data. Other kinds of pre-processing data are possible. 

 dataexample  <- dataexample[-1] #remove the dates!
 DAXCAC       <- dataexample[,c(3,4)] 
 lrdatex      <- apply(log(DAXCAC), 2, diff)
 inputDATA    <- ts(lrdatex, start=1, frequency=1)

 Wname     <- "la8"
 J         <- 8
 lmax      <- 30
 compWCC   <- WCC(inputDATA, Wname, J, lmax, device="screen", NULL,
               NULL, NULL, NULL, NULL)
}
\keyword{bi-variate}
\keyword{Wavelet cross correlation }
