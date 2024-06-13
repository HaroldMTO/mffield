\name{FcTime-class}
\Rdversion{1.1}
\docType{class}
\alias{FcTime-class}

\title{Class \code{"FcTime"}}
\description{
A class defining forecast (Fc) time, that is a pair of one base date and forecast steps \
or the opposite, several base dates and one forecast step.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FcTime", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{base}:}{Object of class \code{"POSIXct"} ~~ }
    \item{\code{step}:}{Object of class \code{"integer"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "FcTime" in the signature.
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("FcTime")
}
\keyword{classes}
