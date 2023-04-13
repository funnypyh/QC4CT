\name{accrual.QC}
\alias{accrual.QC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identification of accrual distribution
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It provides the distribution of the time between acrruals, considering four reference distributions such as Gamma distribution, exponential distribution, four-parameter Beta distribution, and Uniform distribution. The smaller statistic of Kolmogorov-Smirnov test indicates a shorter distance between empirical and true reference distributions. 
}
\usage{
accrual.QC(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
The accrual time
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
\item{gammatest}{
Statistic and p-value of Kolmogorov-Smirnov test for Gamma distribution.
}
\item{exptest}{
Statistic and p-value of Kolmogorov-Smirnov test for Exponential distribution.
}
\item{betatest}{
Statistic and p-value of Kolmogorov-Smirnov test for four-parameter Beta distribution (whose densitiy is described in https://en.wikipedia.org/wiki/Beta_distribution).
}
\item{uniformtest}{
Statistic and p-value of Kolmogorov-Smirnov test for Uniform distribution, which is done when the test for four-parameter Beta distribution does not work due to parameter restriction. 
}
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
\dontrun{
library(ExtDist) 

y <- runif(100, min=0, max=1)
z <- cumsum(y)
res1 <- accrual.QC(z)
res1

y <- rgamma(100, shape=2, scale=0.5)
z <- cumsum(y)
res2 <- accrual.QC(z)
res2
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
