\name{irtreliability}
\alias{irtreliability}
\title{
Marginal and Test Reliability Coefficients with Item Response Theory
}
\description{
A function to estimate marginal and test reliability from estimated item response theory models.
}
\usage{
irtreliability(input, model, cats, relcoef = "trc", nquad = 49, SE = TRUE)
}

\arguments{
  \item{input}{
An object of class SingleGroupClass from package mirt.
}

  \item{model}{
A character vector indicating the item response theory model used, options are "GPCM" and "3-PL".
}

  \item{cats}{
A numeric vector indicating the number of possible categories for each item.
}

  \item{relcoef}{
A character vector indicting which reliability coefficients to calculate, options are "mrc" for the marginal reliability coefficient and "trc" for the test reliability coefficient.
}

  \item{nquad}{
The number of Gauss-Hermite quadrature points to be used.
}

  \item{SE}{
A logical vector denoting whether the standard errors for the reliability coefficient estimates should be calculated.
}

}
\value{
An S4 object of class 'relout' which includes the following slots
	\item{est}{The estimated coefficient.}
	\item{cov}{The estimated variance.}
	\item{pder}{The partial derivatives of the coefficient with respect to the item parameters.}
	\item{type}{The type of coefficient.}
}
\references{
Andersson, B. and Xin, T. (2018). Large Sample Confidence Intervals for Item Response Theory Reliability Coefficients. \emph{Educational and Psychological Measurement}, 78, 32-45.\cr
Cheng, Y., Yuan, K.-H. and Liu, C. (2012). Comparison of reliability measures under factor analysis and item response theory. \emph{Educational and Psychological Measurement}, 72, 52-67.\cr
Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L. and Reckase, M. D. (1984). Technical guidelines for assessing computerized adaptive tests. \emph{Journal of Educational Measurement}, 21, 347-360.\cr
Kim, S. (2012). A note on the reliability coefficients for item response model-based ability estimates. \emph{Psychometrika}, 77, 153-162.\cr
Kim, S. and Feldt, L. S. (2010). The estimation of the IRT reliability coefficient and its lower and upper bounds, with comparisons to CTT reliability statistics. \emph{Asia Pacific Education Review}, 11, 179-188.}
\author{
\email{bjoern.h.andersson@gmail.com}
}

\examples{
#Generate 2-PL data
set.seed(14)
akX <- runif(15, 0.5, 2)
bkX <- rnorm(15)
data2pl <- matrix(0, nrow = 1000, ncol = 15)

for(i in 1:1000){
	ability <- rnorm(1)
	data2pl[i,1:15] <- (1 / (1 + exp(-akX *(ability - bkX)))) > runif(15)
}

#Estimate the 2-PL IRT model with package mirt
library(mirt)
sim2pl <- mirt(data.frame(data2pl), 1, "gpcm", SE = TRUE)
mrc2pl <- irtreliability(sim2pl, "GPCM", rep(2, 15), relcoef = "mrc")
trc2pl <- irtreliability(sim2pl, "GPCM", rep(2, 15))
}