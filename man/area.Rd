\name{area}
\alias{area}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Roughly calculate the area covered by a domain
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
area(dom1, dom2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dom1}{
%%     ~~Describe \code{dom1} here~~
}
  \item{dom2}{
%%     ~~Describe \code{dom2} here~~
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
function (dom1, dom2) 
{
    diff(dom1@xlim) * diff(dom1@ylim)/(diff(dom2@xlim) * diff(dom2@ylim))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
