\name{degrade}
\alias{degrade}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Degrade (or "downgrade") the resolution of a grid
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degrade(grid, ilat, nlat, nlon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid}{
an object of class Grid
}
  \item{ilat}{
an integer vector of latitude indices to select from grid
}
  \item{nlat}{
an integer value giving the pace of latitudes to keep from grid
}
  \item{nlon}{
an integer value giving the pace of longitudes to keep from grid
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class Grid with a selection of lats and longs
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
new("nonstandardGenericFunction", .Data = function (grid, ilat, 
    nlat, nlon) 
{
    standardGeneric("degrade")
  }, generic = structure("degrade", package = "mffield"), package = "mffield", 
    group = list(grid = "Grid"), valueClass = character(0), signature = c("grid", 
    "ilat", "nlat", "nlon"), default = NULL, skeleton = (function (grid, 
        ilat, nlat, nlon) 
    stop("invalid call in method dispatch to 'degrade' (no default method)", 
        domain = NA))(grid, ilat, nlat, nlon))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
