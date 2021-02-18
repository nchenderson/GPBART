\name{FindDistMat}
\alias{FindDistMat}
\title{
BART-based Distances
}
\description{
This function computes all pairwise distances between covariate vectors
using the BART-based distance measure.
}
\usage{
FindDistMat(X)
}
\arguments{
  \item{X}{
The n x p design matrix. Each row is assumed to hold the covariate vector from one observation.
}
\item{tau.max}{
The maximum number of potential splitting points for each variable.
}
\item{weights}{
An optional vector of weights to use when computing the pairwise distances. If provided,
the length of \code{weights} should be equal to the number of columns of \code{X}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A square matrix whose number of rows equals the number of rows in the input matrix \code{X} 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("iris")
X <- as.matrix(iris[,1:4])
Dmat <- FindDistMat(X)
head(Dmat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }% use one of  RShowDoc("KEYWORDS")
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
