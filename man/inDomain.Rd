\name{inDomain}
\alias{inDomain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find gridpoints inside a lat/long domain
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inDomain(grid, domain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{
%%     ~~Describe \code{grid} here~~
}
  \item{domain}{
%%     ~~Describe \code{domain} here~~
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
function (grid, domain) 
{
    ilat = domain@ylim[1] <= grid@lat & grid@lat <= domain@ylim[2]
    xlim = (domain@xlim + 180)\%\%360 - 180
    if (diff(xlim) > 0) {
        ilon = xlim[1] <= grid@long & grid@long <= xlim[2]
    }
    else {
        ilon = grid@long <= xlim[2] | grid@long >= xlim[1]
    }
    ilat & ilon
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
