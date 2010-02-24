\name{GeneralizedHyperbolic}
\alias{dghyp}
\alias{pghyp}
\alias{qghyp}
\alias{rghyp}
\alias{ddghyp}
\alias{ghypBreaks}

\title{Generalized Hyperbolic Distribution}

\description{
  Density function, distribution function, quantiles and
  random number generation for the generalized hyperbolic distribution
  with parameter vector \code{param}. Utility routines are included for
  the derivative of the density function and to find suitable break
  points for use in determining the distribution function.
}

\usage{
dghyp(x, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
      param = c(mu, delta, alpha, beta, lambda))
pghyp(q, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
      param = c(mu, delta, alpha, beta, lambda),
      small = 10^(-6), tiny = 10^(-10),
      deriv = 0.3, subdivisions = 100, accuracy = FALSE, ...)
qghyp(p, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
      param = c(mu, delta, alpha, beta, lambda),
      small = 10^(-6), tiny = 10^(-10),
      deriv = 0.3, nInterpol = 100, subdivisions = 100, ...)
rghyp(n, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
      param = c(mu, delta, alpha, beta, lambda))
ddghyp(x, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
       param = c(mu, delta, alpha, beta, lambda))
ghypBreaks(mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
           param = c(mu, delta, alpha, beta, lambda),
           small = 10^(-6), tiny = 10^(-10), deriv = 0.3, ...)
}

\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations to be generated.}
  \item{mu}{\eqn{\mu}{mu} is the location parameter. By default this is
    set to 0.}
  \item{delta}{\eqn{\delta}{delta} is the scale parameter of the distribution.
    A default value of 1 has been set.}
  \item{alpha}{\eqn{\alpha}{alpha} is the tail parameter, with a default value of 1.}
  \item{beta}{\eqn{\beta}{beta} is the skewness parameter, by default this is 0.}
  \item{lambda}{\eqn{\lambda}{lambda} is the shape parameter and dictates the
    shape that the distribution shall take. Default value is 1.}
  \item{param}{Parameter vector taking the form
    \code{c(mu, delta, alpha, beta, lambda)}.}
  \item{small}{Size of a small difference between the distribution
    function and zero or one. See \bold{Details}.}
  \item{tiny}{Size of a tiny difference between the distribution
    function and zero or one. See \bold{Details}.}
  \item{deriv}{Value between 0 and 1. Determines the point where the
    derivative becomes substantial, compared to its maximum value. See
    \bold{Details}.}
  \item{accuracy}{Uses accuracy calculated by \code{\link{integrate}}
    to try and determine the accuracy of the distribution function
    calculation.}
  \item{subdivisions}{The maximum number of subdivisions used to
    integrate the density returning the distribution function.}
  \item{nInterpol}{The number of points used in qghyp for cubic spline
    interpolation (see \code{splinefun}) of the distribution function.}
  \item{\dots}{Passes arguments to \code{uniroot}.  See \bold{Details}.}
}

\details{
  The generalized hyperbolic distribution has density

  \deqn{f(x)=c(\lambda,\alpha,\beta,\delta)\times%
    \frac{K_{\lambda-1/2}(\alpha\sqrt{\delta^2+(x-\mu)^2})}%
    {(\sqrt{\delta^2+(x-\mu)^2}/\alpha)^{1/2-\lambda}}%
    e^{\beta(x-\mu)}}{%
    f(x)=c(lambda,alpha,beta,delta)%
    (K_(lambda-1/2)(alpha sqrt(delta^2+(x-mu)^2)))/%
    ((sqrt(delta^2+(x-mu)^2)/alpha)^(1/2-lambda))%
    exp(beta(x-mu))}

  where \eqn{K_\nu()}{K_nu()} is the modified Bessel function of the
  third kind with order \eqn{\nu}{nu}, and
  \deqn{c(\lambda,\alpha,\beta,\delta)=%
    \frac{(\sqrt{\alpha^2-\beta^2}/\delta)^\lambda}%
    {\sqrt{2\pi}K_\lambda(\delta\sqrt{\alpha^2-\beta^2})}}{%
    c(lambda,alpha,beta,delta)=%
    (sqrt(alpha^2-beta^2)/delta)^lambda/%
    (sqrt(2\pi)K_lambda(delta sqrt(alpha^2-beta^2)))}

  Use \code{ghypChangePars} to convert from the
  \eqn{(\zeta, \rho)}{(zeta, rho)},
  \eqn{(\xi, \chi)}{xi, chi)}, or
  \eqn{(\bar\alpha, \bar\beta)}{(alpha bar, beta bar)} parameterisations
  to the \eqn{(\alpha, \beta)}{(alpha, beta)} parameterisation used
  above.

  \code{pghyp} breaks the real line into eight regions in order to
  determine the integral of \code{dghyp}. The break points determining
  the regions are found by \code{ghypBreaks}, based on the values of
  \code{small}, \code{tiny}, and \code{deriv}. In the extreme tails of
  the distribution where the probability is \code{tiny} according to
  \code{ghypCalcRange}, the probability is taken to be zero.  In the
  inner part of the distribution, the range is divided in 6 regions, 3
  above the mode, and 3 below. On each side of the mode, there are two
  break points giving the required three regions. The outer break point
  is where the probability in the tail has the value given by the
  variable \code{small}. The inner break point is where the derivative
  of the density function is \code{deriv} times the maximum value of the
  derivative on that side of the mode. In each of the 6 inner regions
  the numerical integration routine
  \code{\link[DistributionUtils]{safeIntegrate}} (which is a wrapper for
  \code{\link{integrate}}) is used to integrate the density
  \code{dghyp}.

  \code{qghyp} use the breakup of the real line into the same 8
  regions as \code{pghyp}. For quantiles which fall in the 2 extreme
  regions, the quantile is returned as \code{-Inf} or \code{Inf} as
  appropriate. In the 6 inner regions \code{splinefun} is used to fit
  values of the distribution function generated by \code{pghyp}. The
  quantiles are then found using the \code{uniroot} function.

  \code{pghyp} and \code{qghyp} may generally be expected to be
  accurate to 5 decimal places.

  The generalized hyperbolic distribution is discussed in Bibby and
  \enc{S�renson}{Sorenson} (2003). It can be represented as a particular
  mixture of the normal distribution where the mixing distribution is the
  generalized inverse Gaussian.  \code{rghyp} uses this representation
  to generate observations from the generalized hyperbolic
  distribution. Generalized inverse Gaussian observations are obtained
  via the algorithm of Dagpunar (1989) which is implemented in
  \code{rgig}.
}

\value{
  \code{dghyp} gives the density, \code{pghyp} gives the distribution
  function, \code{qghyp} gives the quantile function and \code{rghyp}
  generates random variates. An estimate of the accuracy of the
  approximation to the distribution function may be found by setting
  \code{accuracy = TRUE} in the call to \code{pghyp} which then returns
  a list with components \code{value} and \code{error}.

  \code{ddghyp} gives the derivative of \code{dghyp}.

  \code{ghypBreaks} returns a list with components:
   \item{xTiny}{Value such that probability to the left is less than
     \code{tiny}.}
   \item{xSmall}{Value such that probability to the left is less than
     \code{small}.}
   \item{lowBreak}{Point to the left of the mode such that the
     derivative of the density is \code{deriv} times its maximum value
     on that side of the mode.}
   \item{highBreak}{Point to the right of the mode such that the
     derivative of the density is \code{deriv} times its maximum value
     on that side of the mode.}
   \item{xLarge}{Value such that probability to the right is less than
     \code{small}.}
   \item{xHuge}{Value such that probability to the right is less than
     \code{tiny}.}
   \item{modeDist}{The mode of the given generalized hyperbolic distribution.}
 }

\references{
  Barndorff-Nielsen, O. and \enc{Bl�sild}{Blaesild}, P (1983).
  Hyperbolic distributions.
  In \emph{Encyclopedia of Statistical Sciences},
  eds., Johnson, N. L., Kotz, S. and Read, C. B., Vol. 3,
  pp. 700--707. New York: Wiley.

  Bibby, B. M. and \enc{S�renson}{Sorenson}, M. (2003).
  Hyperbolic processes in finance.
  In \emph{Handbook of Heavy Tailed Distributions in Finance}, ed.,
  Rachev, S. T. pp. 212--248. Elsevier Science B.~V.

  Dagpunar, J.S. (1989).
  An easily implemented generalised inverse Gaussian generator
  \emph{Commun. Statist. -Simula.}, \bold{18}, 703--710.

  Prause, K. (1999) \emph{The generalized hyperbolic models: Estimation,
  financial derivatives and risk measurement.} PhD Thesis, Mathematics
  Faculty, University of Freiburg.
}

\author{David Scott \email{d.scott@auckland.ac.nz},
 Richard Trendall}

\seealso{
  \code{\link{dhyperb}} for the hyperbolic distribution,
  \code{\link{dgig}} for the generalized inverse Gaussian distribution
  \code{\link[DistributionUtils]{safeIntegrate}},
  \code{\link{integrate}} for its shortfalls, \code{\link{splinefun}},
  \code{\link{uniroot}} and \code{\link{ghypChangePars}} for
  changing parameters to the \eqn{(\alpha,\beta)}{(alpha,beta)}
  parameterisation
}

\examples{
param <- c(0, 1, 3, 1, 1 / 2)
ghypRange <- ghypCalcRange(param = param, tol = 10^(-3))
par(mfrow = c(1, 2))
curve(dghyp(x, param = param), from = ghypRange[1], to = ghypRange[2],
      n = 1000)
title("Density of the\n Generalized Hyperbolic Distribution")
curve(pghyp(x, param = param), from = ghypRange[1], to = ghypRange[2],
      n = 1000)
title("Distribution Function of the\n Generalized Hyperbolic Distribution")
dataVector <- rghyp(500, param = param)
curve(dghyp(x, param = param), range(dataVector)[1], range(dataVector)[2],
      n = 500)
hist(dataVector, freq = FALSE, add = TRUE)
title("Density and Histogram of the\n Generalized Hyperbolic Distribution")
logHist(dataVector, main =
        "Log-Density and Log-Histogramof the\n Generalized Hyperbolic Distribution")
curve(log(dghyp(x, param = param)), add = TRUE,
      range(dataVector)[1], range(dataVector)[2], n = 500)
par(mfrow = c(2, 1))
curve(dghyp(x, param = param), from = ghypRange[1], to = ghypRange[2],
      n = 1000)
title("Density of the\n Generalized Hyperbolic Distribution")
curve(ddghyp(x, param = param), from = ghypRange[1], to = ghypRange[2],
      n = 1000)
title("Derivative of the Density of the\n Generalized Hyperbolic Distribution")
par(mfrow = c(1, 1))
ghypRange <- ghypCalcRange(param = param, tol = 10^(-6))
curve(dghyp(x, param = param), from = ghypRange[1], to = ghypRange[2],
      n = 1000)
bks <- ghypBreaks(param = param)
abline(v = bks)
}

\keyword{distribution}