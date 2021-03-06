\name{nigHessian}
\alias{nigHessian}

\title{
Calculate Two-Sided Hessian for the Normal Inverse Gaussian Distribution
}
\description{
Calculates the Hessian of a function, either exactly or approximately. Used to
obtaining the information matrix for maximum likelihood estimation.
}
\usage{
nigHessian(x, param, hessianMethod = "tsHessian",
              whichParam = 1:5, ...)
}
\arguments{
  \item{x}{Data vector.}
  \item{param}{The maximum likelihood estimates parameter vector of the
    normal inverse Gaussian distribution. The normal inverse Gaussian 
    distribution has the same sets of parameterizations as the hyperbolic 
    distribution.There are five different sets of
    parameterazations can be used in this function, the first four sets
    are listed in \code{hyperbChangePars} and the last set is the log
    scale of the first set of the parameterization, i.e.,
    \code{mu,log(delta),Pi,log(zeta)}.  }
  \item{hessianMethod}{Only the approximate method (\code{"tsHessian"})
    has actually been implemented so far.}  
  \item{whichParam}{Numeric. A number between 1 to 5 indicating which
    set of the parameterization is the specified value in argument
    \code{param} belong to.} 
  \item{\dots}{Values of other parameters of the function \code{fun} if
    required.} 
}
\details{  
  The approximate Hessian is obtained via a call to \code{tsHessian}
  from the package \code{DistributionUtils}.  \code{summary.nigFit}
  calls the function \code{nigHessian} to calculate the Hessian matrix
  when the argument \code{hessian = TRUE}.  
}
\value{
  \code{nigHessian} gives the approximate or exact Hessian matrix for
  the data vector \code{x} and the estimated parameter vector
  \code{param}.
}
\author{
  David Scott \email{d.scott@auckland.ac.nz},
  Christine Yang Dong \email{c.dong@auckland.ac.nz}
}

\examples{
### Calculate the exact Hessian using nigHessian:
param <- c(2, 2, 2, 1)
dataVector <- rnig(500, param = param)
fit <- nigFit(dataVector, method = "BFGS")
coef=coef(fit)
nigHessian(x=dataVector, param=coef, hessianMethod = "tsHessian",
           whichParam = 2)
              
### Or calculate the exact Hessian using summary.nigFit method:
### summary(fit, hessian = TRUE)

## Calculate the approximate Hessian:
summary(fit, hessian = TRUE, hessianMethod = "tsHessian")
}



