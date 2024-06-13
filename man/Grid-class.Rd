\name{Grid-class}
\Rdversion{1.1}
\docType{class}
\alias{Grid-class}
\alias{[,Grid-method}

\title{
A class defining a linear geographical grid
}
\description{
A linear grid is defined by a pair of vectors of lat and long for each gridpoint.
In the original grid, called the computational sphere, gridpoints are aligned on \
circles of latitude by chunks of varying length.
Term "linear" in this linear grid means that gridpoints are concatenated for all the \
circles of latitude one after the other.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Grid", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{lat}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{long}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "Grid")}: ... }
	 }
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
showClass("Grid")
}
\keyword{classes}
