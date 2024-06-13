\name{csLat}
\alias{csLat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the grid latitudes, or latitudes on the computational sphere
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csLat(grid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{
%%     ~~Describe \code{grid} here~~
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
function (grid) 
{
    C2 = grid@stretch^2
    if (C2 == 1) 
        return(grid@theta)
    mu = sin(grid@theta * pi/180)
    sinla = ((C2 + 1) * mu + C2 - 1)/((C2 - 1) * mu + C2 + 1)
    theta = asin(sinla) * 180/pi
    stopifnot(all(theta > grid@theta))
    theta
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
