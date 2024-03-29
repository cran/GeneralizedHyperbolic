\name{hyperbFit}
\alias{hyperbFit}
\alias{print.hyperbFit}
\alias{plot.hyperbFit}
\alias{coef.hyperbFit}
\alias{vcov.hyperbFit}

\title{Fit the Hyperbolic Distribution to Data}
\description{
  Fits a hyperbolic distribution to data.
  Displays the histogram, log-histogram (both with fitted densities),
  Q-Q plot and P-P plot for the fit which has the maximum likelihood.
}

\usage{
  hyperbFit(x, freq = NULL, paramStart = NULL,
            startMethod = c("Nelder-Mead","BFGS"),
            startValues = c("BN","US","FN","SL","MoM"),
            criterion = "MLE",
            method = c("Nelder-Mead","BFGS","nlm",
                       "L-BFGS-B","nlminb","constrOptim"),
            plots = FALSE, printOut = FALSE,
            controlBFGS = list(maxit = 200),
            controlNM = list(maxit = 1000), maxitNLM = 1500,
            controlLBFGSB = list(maxit = 200),
            controlNLMINB = list(),
            controlCO = list(), ...)

  \method{print}{hyperbFit}(x,
        digits = max(3, getOption("digits") - 3), ...)

  \method{plot}{hyperbFit}(x, which = 1:4,
     plotTitles = paste(c("Histogram of ","Log-Histogram of ",
                        "Q-Q Plot of ","P-P Plot of "), x$obsName,
                        sep = ""),
     ask = prod(par("mfcol")) < length(which) & dev.interactive(), ...)

  \method{coef}{hyperbFit}(object, ...)

  \method{vcov}{hyperbFit}(object, ...)
}

\arguments{
  \item{x}{Data vector for \code{hyperbFit}. Object of class
    \code{"hyperbFit"} for \code{print.hyperbFit} and \code{plot.hyperbFit}.}
  \item{freq}{A vector of weights with length equal to \code{length(x)}.}
  \item{paramStart}{A user specified starting parameter vector
    \code{param} taking the form \code{c(mu, delta, alpha, beta)}.}
  \item{startMethod}{Method used by \code{hyperbFitStart} in calls to
    \code{\link{optim}}.}
  \item{startValues}{Code giving the method of determining starting
    values for finding the maximum likelihood estimate of \code{param}.}
  \item{criterion}{Currently only \code{"MLE"} is implemented.}
  \item{method}{Different optimisation methods to consider.
    See \bold{Details}.}
  \item{plots}{Logical. If \code{FALSE} suppresses printing of the
    histogram, log-histogram, Q-Q plot and P-P plot.}
  \item{printOut}{Logical. If \code{FALSE} suppresses printing of
    results of fitting.}
  \item{controlBFGS}{A list of control parameters for \code{optim} when using
    the \code{"BFGS"} optimisation.}
  \item{controlNM}{A list of control parameters for \code{optim}
    when using the \code{"Nelder-Mead"} optimisation.}
  \item{maxitNLM}{A positive integer specifying the maximum number of
    iterations when using the \code{"nlm"} optimisation.}
  \item{controlLBFGSB}{A list of control parameters for \code{optim} when using
    the \code{"L-BFGS-B"} optimisation.}
  \item{controlNLMINB}{A list of control parameters for \code{nlminb}
    when using the \code{"nlminb"} optimisation.}
  \item{controlCO}{A list of control parameters for \code{constrOptim}
    when using the \code{"constrOptim"} optimisation.}
  \item{digits}{Desired number of digits when the object is printed.}
   \item{which}{If a subset of the plots is required, specify a subset of
    the numbers \code{1:4}.}
  \item{plotTitles}{Titles to appear above the plots.}
  \item{ask}{Logical. If \code{TRUE}, the user is \emph{ask}ed before
    each plot, see \code{\link{par}(ask = .)}.}
  \item{\dots}{Passes arguments to \code{par}, \code{hist},
    \code{logHist}, \code{qqhyperb} and \code{pphyperb}.}
  \item{object}{Object of class \code{"hyperbFit"} for \code{coef.hyperbFit}
    and for \code{vcov.hyperbFit}.}
}
\details{
  \code{startMethod} can be either \code{"BFGS"} or
	\code{"Nelder-Mead"}.

  \code{startValues} can be one of the following:
  \describe{
    \item{\code{"US"}}{User-supplied.}
    \item{\code{"BN"}}{Based on Barndorff-Nielsen (1977).}
    \item{\code{"FN"}}{A fitted normal distribution.}
    \item{\code{"SL"}}{Based on a fitted skew-Laplace distribution.}
    \item{\code{"MoM"}}{Method of moments.}
  }
  For the details concerning the use of \code{paramStart},
  \code{startMethod}, and \code{startValues}, see
  \code{\link{hyperbFitStart}}.

  The six optimisation methods currently available are:
  \describe{
  \item{\code{"BFGS"}}{Uses the quasi-Newton method \code{"BFGS"} as
    documented in \code{\link{optim}}.}
  \item{\code{"Nelder-Mead"}}{Uses an implementation of the Nelder and
    Mead method as documented in \code{\link{optim}}.}
  \item{\code{"nlm"}}{Uses the \code{\link{nlm}} function in R.}
  \item{\code{"L-BFGS-B"}}{Uses the quasi-Newton method with box
    constraints \code{"L-BFGS-B"} as documented in \code{\link{optim}}.}
  \item{\code{"nlminb"}}{Uses the \code{\link{nlminb}} function in R.}
  \item{\code{"constrOptim"}}{Uses the \code{\link{constrOptim}}
    function in R.}
  }

  For details of how to pass control information for optimisation using
  \code{optim}, \code{nlm}, \code{nlminb} and \code{constrOptim}, see
  \code{\link{optim}}, \code{\link{nlm}}, \code{\link{nlminb}} and
  \code{\link{constrOptim}}.

  When \code{method = "nlm"} is used, warnings may be produced. These do
  not appear to be a problem.
}
\value{
  \code{hyperbFit} returns a list with components:
    \item{param}{A vector giving the maximum likelihood estimate of
      param, as \code{c(mu, delta, alpha, beta)}.}
    \item{maxLik}{The value of the maximised log-likelihood.}
    \item{method}{Optimisation method used.}
    \item{conv}{Convergence code. See the relevant documentation (either
      \code{\link{optim}} or \code{\link{nlm}}) for details on
      convergence.}
    \item{iter}{Number of iterations of optimisation routine.}
    \item{obs}{The data used to fit the hyperbolic distribution.}
    \item{obsName}{A character string with the actual \code{x} argument
      name.}
    \item{paramStart}{Starting value of param returned by call to
      \code{\link{hyperbFitStart}}.}
    \item{svName}{Descriptive name for the method finding start values.}
    \item{startValues}{Acronym for the method of finding start values.}
    \item{breaks}{The cell boundaries found by a call to
      \code{\link{hist}}.}
    \item{midpoints}{The cell midpoints found by a call to
      \code{\link{hist}}.}
    \item{empDens}{The estimated density found by a call to
      \code{\link{hist}}.}
  }

\references{
  Barndorff-Nielsen, O. (1977)
  Exponentially decreasing distributions for the logarithm of particle size,
  \emph{Proc. Roy. Soc. Lond.} A\bold{353}, 401--419.

  Fieller, N. J., Flenley, E. C. and Olbricht, W. (1992)
  Statistics of particle size data.
  \emph{Appl. Statist.} \bold{41}, 127--146.
}
\author{David Scott \email{d.scott@auckland.ac.nz},
  Ai-Wei Lee, Jennifer Tso, Richard Trendall, Thomas Tran,
  Christine Yang Dong
}
\seealso{
  \code{\link{optim}}, \code{\link{nlm}}, \code{\link{nlminb}},
  \code{\link{constrOptim}}, \code{\link{par}},
  \code{\link{hist}}, \code{\link[DistributionUtils]{logHist}} (pkg \CRANpkg{DistributionUtils}),
  \code{\link{qqhyperb}}, \code{\link{pphyperb}}, \code{\link{dskewlap}}
  and \code{\link{hyperbFitStart}}.
}

\examples{
param <- c(2, 2, 2, 1)
dataVector <- rhyperb(500, param = param)
## See how well hyperbFit works
hyperbFit(dataVector)
hyperbFit(dataVector, plots = TRUE)
fit <- hyperbFit(dataVector)
par(mfrow = c(1, 2))
plot(fit, which = c(1, 3))

## Use nlm instead of default
hyperbFit(dataVector, method = "nlm")

}

\keyword{distribution}
