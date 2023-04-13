\name{survival.QC2}
\alias{survival.QC2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identification of survival distribution (Option 2)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It provides the distribution of survival time, considering two reference distributions such as Exponential distribution and Weibull distribution. This works when the observed survival time and accrual time are given. 
}
\usage{
survival.QC2(surv_time, arr_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv_time}{
%%     ~~Describe \code{surv_time} here~~
Survival time
}
  \item{arr_time}{
%%     ~~Describe \code{arr_time} here~~
Accrual time
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
\item{exptest}{
The weighted Kaplan-Meier statistics of Pepe and Fleming (1989) performing distance test for Exponential distribution.
}
\item{weibulltest}{
The weighted Kaplan-Meier statistics of Pepe and Fleming (1989) performing distance test for Weibull distribution.
}
}
\references{
%% ~put references to the literature/web site here ~
Pepe, M. S. and Fleming, T. R. (1989). Weighted Kaplan-Meier statistics: a class of distance tests for censored survival data. Biometrics, 497-507.

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
library(cubature) 
library(survival) 

data(data2)
res <- survival.QC2(surv_time=data2[,1], arr_time=data2[,2])
res
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
