\name{rangeLong}
\alias{rangeLong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Give the range of longitudes accounting for possible wrap around 180
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rangeLong(long)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long}{
%%     ~~Describe \code{long} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (long) 
{
    dlong1 = range(long\%\%360)
    dlong2 = range((long + 180)\%\%360 - 180)
    if (diff(dlong1) > 180 && diff(dlong2) > 180) 
        stop("longitudes span more than 180dg (case not supported)")
    if (diff(dlong1) <= diff(dlong2)) {
        return(dlong1)
    }
    else {
        return(dlong2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
