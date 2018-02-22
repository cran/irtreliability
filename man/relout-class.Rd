\name{relout-class}
\Rdversion{1.1}
\docType{class}
\alias{relout-class}

\title{Class "relout"}
\description{
Estimates of reliability coefficients and coefficients related to the reliability.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("keout", ...)}.
}
\section{Slots}{
  \describe{
	\item{est}{The estimated coefficient.}
	\item{cov}{The estimated variance.}
	\item{pder}{The partial derivatives of the coefficient with respect to the item parameters.}
	\item{type}{The type of coefficient.}
    }
}

\author{
\email{bjoern.h.andersson@gmail.com}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("relout")
}
\keyword{classes}
