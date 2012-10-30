\encoding{latin1}
\name{WMCC}
\alias{WMCC}
\alias{wmcc}
\title{Wavelet Multiple cross-correlation (multivariate case).
}
\description{
The \emph{WMCC} function (multivariate case) computes the Wavelet 
Multiple cross correlation by means of the function 
\emph{wave.multiple.cross.correlation} from the \emph{wavemulcor} 
package (Fernandez-Macho 2012) and present the result as a novel 
plot that reduce the number of plots of the classical function 
\emph{wave.multiple.cross.correlation}. The \emph{WMC} function also 
provides a way to handle multivariate time series easily as a list of 
N elements (time series). 
}
\usage{
WMCC(inputDATA, Wname, J, lmax, Hpdf, Wpdf)
}
\arguments{
  \item{inputDATA}{
An array of multivariate time series as a \emph{ts} object (please,
check the \emph{ts} manual to get more information about the \emph{
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
The \emph{WMCC} function compute the Wavelet Multiple 
Cross Correlation using the function \cr 
\emph{wave.multiple.cross.correlation} from the \emph{wavemulcor} 
package (Fernandez-Macho 2012), but the \emph{WMCC} 
function incorporates some graphical improvements (please, see the 
Figure 5 of Polanco-Martinez and Fernandez-Macho 2012), such as 
the reduction of the number of plots to present the results of the 
function \emph{wave.multiple.cross.correlation}. 
}
\value{
Output:\cr
Output pdf file: \emph{WMCC_plot.pdf}.

Output data: The same list of elements of the function 
\emph{wave.multiple.cross.correlation} of the \emph{wavemulcor} 
package (Fernandez-Macho 2012). 
}
\references{
Fernandez-Macho, J.. Wavelet multiple correlation and
cross-correlation: A multiscale analysis of euro zone stock
markets. \emph{Physica A: Statistical Mechanics and its Applications}, 
391(4):1097-1104, 2012. \cr

Fernandez-Macho, J. (2012). Wavemulcor Reference manual. The 
Comprehensive R Archive Network (CRAN), 
http://cran.r-project.org/web/packages/wavemulcor/index.html.\cr

Polanco-Martinez, J. and J. Fernandez-Macho (2012). The
package 'W2CWM2C': description, features and applications.
To be submitted under review to \emph{Journal of
Statistical Software}.
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
Needs \emph{wavemulcor} (to compute the \emph{wave.multiple.cross.correlation})
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
 lrdatex      <- apply(log(dataexample), 2, diff)
 tslrdat      <- ts(lrdatex, start=1, frequency=1)

 Wname   <- "la8"
 J       <- 8
 Hp      <- 6
 Wp      <- 10
 lmax    <- 30
 compWCC <- WMCC(tslrdat, Wname, J, lmax, Hp, Wp)
}
\keyword{Multi-variate }
\keyword{Wavelet Multiple cross-correlation}
