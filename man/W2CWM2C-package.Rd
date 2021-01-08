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
'spin.correlation' (wavelet cross correlation) (\emph{waveslim} package, 
Whitcher 2012) and the 'wave.multiple.correlation' and 
'wave.multiple.cross.correlation' (\emph{wavemulcor} package, Fernandez-Macho 
2012b). The plot outputs (heatmaps) can be displayed in the screen or 
can be saved as PNG or JPG images or as PDF or EPS formats. The W2CWM2C 
package also helps to handle the (input data) multivariate time series 
easily as a list of N elements (times series) and provides a multivariate 
data set (dataexample) to exemplify its use. A description of the package was published by \strong{Polanco-Martinez} and Fernandez-Macho (2014), <doi:10.1109/MCSE.2014.96>. 
}
\details{
\tabular{ll}{
Package: \tab W2CWM2C\cr
Type: \tab Package\cr
Version: \tab 2.2\cr
Date: \tab 2021-01-07\cr
License: \tab GPL (>= 2) \cr
LazyLoad: \tab yes\cr
}
The  W2CWM2C package contains four functions: (1) \code{\link{WC}} that 
performs and plots the wavelet correlation for the bivariate case, (2) 
\code{\link{WCC}} that performs and plots the wavelet cross correlation
for the bivariate case, (3) \code{\link{WMC}} that performs and plots 
the wavelet multiple correlation for the multivariate case, and 
4) \code{\link{WMCC}} that performs and plots the wavelet multiple cross 
correlation for the multivariate case. 
}
\note{
 Dependencies: \emph{waveslim}, \emph{wavemulcor} and \emph{colorspace}. 
}
\author{
Josue M. Polanco-Martinez (a.k.a. jomopo).\cr
BC3 - Basque Centre for Climate Change, Bilbao, Spain.\cr 
Web1: \url{https://scholar.google.es/citations?user=8djLIhcAAAAJ&hl=en}.\cr
Web2: \url{https://www.researchgate.net/profile/Josue_Polanco-Martinez}.\cr
Email: josue.m.polanco@gmail.com \cr\cr
Acknowledgement:\cr 
I am thankful to Jan de Leeuw and Achim Zeileis for suggesting the use 
of the R package \emph{colorspace}. I am also thanks to Debojyoti Das, 
John Garrigan, Peterson Owusu Junior, Rim Ibrahim, and Ato Wilberforce
for reporting some bugs found in my W2CWM2C package. The author 
acknowledges to the \strong{SEPE} (Spanish Public Service of Employment) 
for its funding support.
}
\references{
Fernandez-Macho, J. (2012a). Wavelet multiple correlation and 
cross-correlation: A multiscale analysis of Euro zone stock 
markets. \emph{Physica A: Statistical Mechanics and its Applications}, \cr
\strong{391}(4):1097--1104. {\doi{10.1016/j.physa.2011.11.002}}. \cr

Fernandez-Macho, J. (2012b). \emph{wavemulcor: Wavelet routine for 
multiple correlation}. R package version 1.2, The Comprehensive R 
Archive Network (CRAN), <URL: \url{https://cran.r-project.org/package=wavemulcor}>. \cr 

Gencay, R., F. Selcuk and B. Whitcher (2001) \emph{An 
Introduction to Wavelets and Other Filtering Methods in Finance 
and Economics}, Academic Press. \cr 

Ihaka, R., Murrell, P., Hornik, K., Fisher, J. C. and Zeileis, A. 
(2012). \emph{colorspace: Color Space Manipulation}. R package version
1.2.0, The Comprehensive R Archive Network (CRAN), <URL: \url{https://cran.r-project.org/package=colorspace}>. \cr

\strong{Polanco-Martinez, J.} and J. Fernandez-Macho (2014). The
package 'W2CWM2C': description, features and applications. 
\emph{Computing in Science & Engineering}, \strong{16}(6):68--78. \cr 
{\doi{10.1109/MCSE.2014.96}}. \cr

\strong{Polanco-Martinez, J. M.} and Abadie, L. M. (2016). Analyzing 
crude oil spot price dynamics versus long term future prices: A 
wavelet analysis approach. \emph{Energies}, \strong{9}(12), 1089.  
{\doi{10.3390/en9121089}}. \cr

Whitcher, B., P. Guttorp, and D.B (2000). Percival. Wavelet analysis of 
covariance with application to atmospheric time series. \emph{
Journal of Geophysical Research - Atmospheres}, \strong{105}(D11):941--962. 
{\doi{10.1029/2000JD900110}}. \cr

Whitcher, B. (2012). \emph{waveslim: Basic wavelet routines for one-, 
two- and three-dimensional signal processing}. R package version 1.7.1. 
The Comprehensive R Archive Network (CRAN), \cr
<URL: \url{https://cran.r-project.org/package=waveslim}>. \cr

Zeileis A, Hornik K, Murrell P (2009). Escaping RGBland: Selecting 
Colors for Statistical Graphics. \emph{Computational Statistics & 
Data Analysis}, \strong{53}, 3259--3270. {\doi{10.1016/j.csda.2008.11.033}}. 
}
