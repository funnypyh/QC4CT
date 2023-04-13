\name{survival.QC1}
\alias{survival.QC1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Identification of survival distribution (Option 1)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It provides the distribution of survival time, considering two reference distributions such as Exponential distribution and Weibull distribution. This works when the observed survival time and censoring indicator are given. 
}
\usage{
survival.QC1(surv_time, ind, study_time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{surv_time}{
%%     ~~Describe \code{surv_time} here~~
Survival time
}
  \item{ind}{
%%     ~~Describe \code{ind} here~~
Censoring indicator taking 1 for the observed survival time and 0 for the censored survival time
}
  \item{study_time}{
%%     ~~Describe \code{study_time} here~~
Period of clinical trial
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
library(cubature) 
library(survival) 

data(data1)
study_time <- max(data1[,1])+1
res <- survival.QC1(surv_time=data1[,1], ind=data1[,2], study_time)
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
