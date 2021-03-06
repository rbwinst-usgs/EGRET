% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plotResidTime.R
\name{plotResidTime}
\alias{plotResidTime}
\title{Plot of the residuals from WRTDS (in log concentration units) versus time}
\usage{
plotResidTime(eList, stdResid = FALSE, printTitle = TRUE, hLine = TRUE,
  tinyPlot = FALSE, col = "black", lwd = 1, cex = 0.8, cex.axis = 1.1,
  cex.main = 1.1, customPar = FALSE, ...)
}
\arguments{
\item{eList}{named list with at least the Sample and INFO dataframes}

\item{stdResid}{logical variable, if TRUE it uses the standardized residual, if FALSE it uses the actual, default is FALSE}

\item{printTitle}{logical variable if TRUE title is printed, if FALSE not printed (this is best for a multi-plot figure)}

\item{hLine}{inserts horizontal line at zero}

\item{tinyPlot}{logical variable, if TRUE plot is designed to be plotted small, as a part of a multipart figure, default is FALSE}

\item{col}{color of points on plot, see ?par 'Color Specification'}

\item{lwd}{number line width}

\item{cex}{numerical value giving the amount by which plotting symbols should be magnified}

\item{cex.axis}{magnification to be used for axis annotation relative to the current setting of cex}

\item{cex.main}{magnification to be used for main titles relative to the current setting of cex}

\item{customPar}{logical defaults to FALSE. If TRUE, par() should be set by user before calling this function
(for example, adjusting margins with par(mar=c(5,5,5,5))). If customPar FALSE, EGRET chooses the best margins depending on tinyPlot.}

\item{\dots}{arbitrary graphical parameters that will be passed to genericEGRETDotPlot function (see ?par for options)}
}
\description{
This function produces a plot of the residuals from WRTDS, expressed in natural log concentration units
versus time.
It also provides an alternative for viewing the standardized residuals, where the each residual is divided by its estimated standard error.

 Although there are a lot of optional arguments to this function, most are set to a logical default.

Data come from named list, which contains a Sample dataframe with the sample data,
and an INFO dataframe with metadata.
}
\examples{
eList <- Choptank_eList
# Water year:
plotResidTime(eList)
# Graphs consisting of Jun-Aug
eList <- setPA(eList, paStart=6,paLong=3)
plotResidTime(eList)
}
\seealso{
\code{\link{selectDays}}, \code{\link{genericEGRETDotPlot}}
}
\keyword{graphics}
\keyword{statistics}
\keyword{water-quality}

