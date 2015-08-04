\encoding{latin1}
\name{W2CWM2C-package}
\alias{W2CWM2C-package}
\alias{W2CWM2C}
\docType{package}
\title{W2CWM2C: a graphical tool for wavelet (cross) correlation and 
wavelet multiple (cross) correlation analysis 
}
\description{
The W2CWM2C package is a set of functions that improves the 
graphical presentations of the functions 'wave.correlation' and
'spin.correlation' (wavelet cross correlation ) (\emph{waveslim} package, 
Whitcher 2012) and the 'wave.multiple.correlation' and 
'wave.multiple.cross.correlation' (\emph{wavemulcor} package, Fernandez-Macho 
2012b). The plot outputs (heatmaps) can be displayed in the screen or 
can be saved as PNG or JPG images or as PDF or EPS formats. The W2CWM2C 
package also helps to handle the (input data) multivariate time series 
easily as a list of N elements (times series) and provides a multivariate 
data set (dataexample) to exemplify its use. A description of the package was 
published in Computing in Science & Engineering (Volume:16,  Issue: 6) on Sep. 
09, 2014, \emph{doi:10.1109/MCSE.2014.96}
}
\details{
\tabular{ll}{
Package: \tab W2CWM2C\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-08-04\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
The  W2CWM2C package contains four functions: the \code{\link{WC}} (that 
perform and plot the wavelet correlation, bivariate case), the 
\code{\link{WCC}} (that perform and plot the wavelet cross correlation, 
bivariate case), the \code{\link{WMC}} (that perform and plot the wavelet 
multiple correlation, multivariate case) and the \code{\link{WMCC}} (that 
perform and plot the wavelet multiple cross correlation, 
multivariate case). 
}
\note{
 Dependencies: \emph{waveslim}, \emph{wavemulcor} and \emph{colorspace}. 
}
\author{
Josue M. Polanco-Martinez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, Spain.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{http://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: josue.m.polanco@gmail.com \cr\cr
Acknowledgement:\cr 
We are thankful to Jan de Leeuw and Achim Zeileis for suggesting the use 
of the R package \emph{colorspace}. 
}
\references{
Fernandez-Macho, J. (2012a). Wavelet multiple correlation and 
cross-correlation: A multiscale analysis of euro zone stock 
markets. \emph{Physica A: Statistical Mechanics and its Applications}, 
391(4):1097--1104. \cr

Fernandez-Macho, J. (2012b). \emph{wavemulcor: Wavelet routine for 
multiple correlation}. R package version 1.2. The Comprehensive R 
Archive Network (CRAN), \cr
\url{http://cran.r-project.org/package=wavemulcor} \cr 

Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An 
Introduction to Wavelets and Other Filtering Methods in Finance 
and Economics}, Academic Press.\cr 

Ihaka, R., Murrell, P., Hornik, K., Fisher, J. C. and Zeileis, A. 
(2012). \emph{colorspace: Color Space Manipulation}. R package version
1.2.0, The Comprehensive R Archive Network (CRAN), \cr
\url{http://cran.r-project.org/package=colorspace} \cr

Polanco-Martinez, J. and J. Fernandez-Macho (2014). The
package 'W2CWM2C': description, features and applications. 
\emph{Computing in Science & Engineering},  16(6):68--78, 
\href{http://www.computer.org/csdl/mags/cs/preprint/06894486.pdf}{DOI: 10.1109/MCSE.2014.96}.\cr

Whitcher, B., P. Guttorp, and D.B (2000). Percival. Wavelet analysis of 
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, 105(D11):941--962.\cr

Whitcher, B. (2012). \emph{waveslim: Basic wavelet routines for one-, 
two- and three-dimensional signal processing}. R package version 1.7.1. 
The Comprehensive R Archive Network (CRAN), \cr
\url{http://cran.r-project.org/package=waveslim} \cr

Zeileis A, Hornik K, Murrell P (2009). Escaping RGBland: Selecting 
Colors for Statistical Graphics. \emph{Computational Statistics & 
Data Analysis}, 53, 3259--3270.
}
