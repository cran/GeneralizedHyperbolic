pkgname <- "GeneralizedHyperbolic"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('GeneralizedHyperbolic')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("SandP500")
### * SandP500

flush(stderr()); flush(stdout())

### Name: SandP500
### Title: S&P 500
### Aliases: SandP500
### Keywords: datasets

### ** Examples

data(SandP500)
### Consider proportional changes in the index
change <- SandP500[-length(SandP500)] / SandP500[-1]
hist(change)
### Fit hyperbolic distribution to changes
hyperbFit(change)




cleanEx()
nameEx("dghyp")
### * dghyp

flush(stderr()); flush(stdout())

### Name: GeneralizedHyperbolicDistribution
### Title: Generalized Hyperbolic Distribution
### Aliases: dghyp pghyp qghyp rghyp ddghyp
### Keywords: distribution

### ** Examples

param <- c(0, 1, 3, 1, 1/2)
ghypRange <- ghypCalcRange(param = param, tol = 10^(-3))
par(mfrow = c(1, 2))

### curves of density and distribution
curve(dghyp(x, param = param), ghypRange[1], ghypRange[2], n = 1000)
title("Density of the \n Generalized  Hyperbolic Distribution")
curve(pghyp(x, param = param), ghypRange[1], ghypRange[2], n = 500)
title("Distribution Function of the \n Generalized Hyperbolic Distribution")

### curves of density and log density
par(mfrow = c(1, 2))
data <- rghyp(1000, param = param)
curve(dghyp(x, param = param), range(data)[1], range(data)[2],
      n = 1000, col = 2)
hist(data, freq = FALSE, add = TRUE)
title("Density and Histogram of the\n Generalized Hyperbolic Distribution")
logHist(data, main = "Log-Density and Log-Histogram of\n the Generalized
      Hyperbolic Distribution")
curve(log(dghyp(x, param = param)),
      range(data)[1], range(data)[2],
      n = 500, add = TRUE, col = 2)

### plots of density and derivative
par(mfrow = c(2, 1))
curve(dghyp(x, param = param), ghypRange[1], ghypRange[2], n = 1000)
title("Density of the\n Generalized  Hyperbolic Distribution")
curve(ddghyp(x, param = param), ghypRange[1], ghypRange[2], n = 1000)
title("Derivative of the Density of the\n Generalized Hyperbolic Distribution")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dgig")
### * dgig

flush(stderr()); flush(stdout())

### Name: Generalized Inverse Gaussian
### Title: Generalized Inverse Gaussian Distribution
### Aliases: dgig pgig qgig rgig rgig1 ddgig
### Keywords: distribution

### ** Examples

param <- c(2, 3, 1)
gigRange <- gigCalcRange(param = param, tol = 10^(-3))
par(mfrow = c(1, 2))
curve(dgig(x, param = param), from = gigRange[1], to = gigRange[2],
      n = 1000)
title("Density of the\n Generalized Inverse Gaussian")
curve(pgig(x, param = param), from = gigRange[1], to = gigRange[2],
      n = 1000)
title("Distribution Function of the\n Generalized Inverse Gaussian")
dataVector <- rgig(500, param = param)
curve(dgig(x, param = param), range(dataVector)[1], range(dataVector)[2],
      n = 500)
hist(dataVector, freq = FALSE, add = TRUE)
title("Density and Histogram\n of the Generalized Inverse Gaussian")
logHist(dataVector, main =
        "Log-Density and Log-Histogram\n of the Generalized Inverse Gaussian")
curve(log(dgig(x, param = param)), add = TRUE,
      range(dataVector)[1], range(dataVector)[2], n = 500)
par(mfrow = c(2, 1))
curve(dgig(x, param = param), from = gigRange[1], to = gigRange[2],
      n = 1000)
title("Density of the\n Generalized Inverse Gaussian")
curve(ddgig(x, param = param), from = gigRange[1], to = gigRange[2],
      n = 1000)
title("Derivative of the Density\n of the Generalized Inverse Gaussian")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dhyperb")
### * dhyperb

flush(stderr()); flush(stdout())

### Name: Hyperbolic
### Title: Hyperbolic Distribution
### Aliases: dhyperb phyperb qhyperb rhyperb ddhyperb
### Keywords: distribution

### ** Examples

param <- c(0, 2, 1, 0)
hyperbRange <- hyperbCalcRange(param = param, tol = 10^(-3))
par(mfrow = c(1, 2))
curve(dhyperb(x, param = param), from = hyperbRange[1], to = hyperbRange[2],
      n = 1000)
title("Density of the\n Hyperbolic Distribution")
curve(phyperb(x, param = param), from = hyperbRange[1], to = hyperbRange[2],
      n = 1000)
title("Distribution Function of the\n Hyperbolic Distribution")
dataVector <- rhyperb(500, param = param)
curve(dhyperb(x, param = param), range(dataVector)[1], range(dataVector)[2],
      n = 500)
hist(dataVector, freq = FALSE, add =TRUE)
title("Density and Histogram\n of the Hyperbolic Distribution")
logHist(dataVector, main =
        "Log-Density and Log-Histogram\n of the Hyperbolic Distribution")
curve(log(dhyperb(x, param = param)), add = TRUE,
      range(dataVector)[1], range(dataVector)[2], n = 500)
par(mfrow = c(2, 1))
curve(dhyperb(x, param = param), from = hyperbRange[1], to = hyperbRange[2],
      n = 1000)
title("Density of the\n Hyperbolic Distribution")
curve(ddhyperb(x, param = param), from = hyperbRange[1], to = hyperbRange[2],
      n = 1000)
title("Derivative of the Density\n of the Hyperbolic Distribution")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dnig")
### * dnig

flush(stderr()); flush(stdout())

### Name: NIG
### Title: Normal Inverse Gaussian Distribution
### Aliases: dnig pnig qnig rnig ddnig
### Keywords: distribution

### ** Examples

param <- c(0, 2, 1, 0)
nigRange <- nigCalcRange(param = param, tol = 10^(-3))
par(mfrow = c(1, 2))
curve(dnig(x, param = param), from = nigRange[1], to = nigRange[2],
      n = 1000)
title("Density of the\n Normal Inverse Gaussian Distribution")
curve(pnig(x, param = param), from = nigRange[1], to = nigRange[2],
      n = 1000)
title("Distribution Function of the\n Normal Inverse Gaussian Distribution")
dataVector <- rnig(500, param = param)
curve(dnig(x, param = param), range(dataVector)[1], range(dataVector)[2],
      n = 500)
hist(dataVector, freq = FALSE, add =TRUE)
title("Density and Histogram\n of the Normal Inverse Gaussian Distribution")
logHist(dataVector, main =
        "Log-Density and Log-Histogram\n of the Normal Inverse Gaussian Distribution")
curve(log(dnig(x, param = param)), add = TRUE,
      range(dataVector)[1], range(dataVector)[2], n = 500)
par(mfrow = c(2, 1))
curve(dnig(x, param = param), from = nigRange[1], to = nigRange[2],
      n = 1000)
title("Density of the\n Normal Inverse Gaussian Distribution")
curve(ddnig(x, param = param), from = nigRange[1], to = nigRange[2],
      n = 1000)
title("Derivative of the Density\n of the Normal Inverse Gaussian Distribution")



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("dskewlap")
### * dskewlap

flush(stderr()); flush(stdout())

### Name: SkewLaplace
### Title: Skew-Laplace Distribution
### Aliases: dskewlap pskewlap qskewlap rskewlap
### Keywords: distribution

### ** Examples

param <- c(1, 1, 2)
par(mfrow = c(1, 2))
curve(dskewlap(x, param = param), from = -5, to = 8, n = 1000)
title("Density of the\n Skew-Laplace Distribution")
curve(pskewlap(x, param = param), from = -5, to = 8, n = 1000)
title("Distribution Function of the\n Skew-Laplace Distribution")
dataVector <- rskewlap(500, param = param)
curve(dskewlap(x, param = param), range(dataVector)[1], range(dataVector)[2],
      n = 500)
hist(dataVector, freq = FALSE, add = TRUE)
title("Density and Histogram\n of the Skew-Laplace Distribution")
logHist(dataVector, main =
        "Log-Density and Log-Histogram\n of the Skew-Laplace Distribution")
curve(log(dskewlap(x, param = param)), add = TRUE,
      range(dataVector)[1], range(dataVector)[2], n = 500)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ghypCalcRange")
### * ghypCalcRange

flush(stderr()); flush(stdout())

### Name: ghypCalcRange
### Title: Range of a Generalized Hyperbolic Distribution
### Aliases: ghypCalcRange
### Keywords: distribution

### ** Examples

param <- c(0, 1, 5, 3, 1)
maxDens <- dghyp(ghypMode(param = param), param = param)
ghypRange <- ghypCalcRange(param = param, tol = 10^(-3) * maxDens)
ghypRange
curve(dghyp(x, param = param), ghypRange[1], ghypRange[2])
## Not run: ghypCalcRange(param = param, tol = 10^(-3), density = FALSE)



cleanEx()
nameEx("ghypChangePars")
### * ghypChangePars

flush(stderr()); flush(stdout())

### Name: ghypChangePars
### Title: Change Parameterizations of the Generalized Hyperbolic
###   Distribution
### Aliases: ghypChangePars
### Keywords: distribution

### ** Examples

param1 <- c(0, 3, 2, 1, 2)               # Parameterization 1
param2 <- ghypChangePars(1, 2, param1)   # Convert to parameterization 2
param2                                   # Parameterization 2
ghypChangePars(2, 1, param2)             # Back to parameterization 1



cleanEx()
nameEx("ghypCheckPars")
### * ghypCheckPars

flush(stderr()); flush(stdout())

### Name: ghypCheckPars
### Title: Check Parameters of the Generalized Hyperbolic Distribution
### Aliases: ghypCheckPars
### Keywords: distribution

### ** Examples

ghypCheckPars(c(0, 2.5, -0.5, 1, 0))      # error
ghypCheckPars(c(0, 2.5, 0.5, 0, 0))       # normal
ghypCheckPars(c(0, 1, 1, -1, 0))          # error
ghypCheckPars(c(2, 0, 1, 0.5, 0))         # error
ghypCheckPars(c(0, 5, 2, 1.5, 0))         # normal
ghypCheckPars(c(0, -2.5, -0.5, 1, 1))     # error
ghypCheckPars(c(0, -1, 0.5, 1, 1))        # error
ghypCheckPars(c(0, 0, -0.5, -1, 1))       # error
ghypCheckPars(c(2, 0, 0.5, 0, -1))        # error
ghypCheckPars(c(2, 0, 1, 0.5, 1))         # skew laplace
ghypCheckPars(c(0, 1, 1, 1, -1))          # skew hyperbolic



cleanEx()
nameEx("ghypMeanVarMode")
### * ghypMeanVarMode

flush(stderr()); flush(stdout())

### Name: Specific Generalized Hyperbolic Moments and Mode
### Title: Moments and Mode of the Generalized Hyperbolic Distribution
### Aliases: ghypMean ghypVar ghypSkew ghypKurt ghypMode
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1, 2)
ghypMean(param = param)
ghypVar(param = param)
ghypSkew(param = param)
ghypKurt(param = param)
ghypMode(param = param)
maxDens <- dghyp(ghypMode(param = param), param = param)
ghypRange <- ghypCalcRange(param = param, tol = 10^(-3) * maxDens)
curve(dghyp(x, param = param), ghypRange[1], ghypRange[2])
abline(v = ghypMode(param = param), col = "blue")
abline(v = ghypMean(param = param), col = "red")



cleanEx()
nameEx("ghypMom")
### * ghypMom

flush(stderr()); flush(stdout())

### Name: ghypMom
### Title: Calculate Moments of the Generalized Hyperbolic Distribution
### Aliases: ghypMom
### Keywords: distribution

### ** Examples

param <- c(1, 2, 2, 1, 2)
mu <- param[1]
### mu moments
m1 <- ghypMean(param = param)
m1 - mu
ghypMom(1, param = param, momType = "mu")
momIntegrated("ghyp", order = 1, param = param, about = mu)
ghypMom(2, param = param, momType = "mu")
momIntegrated("ghyp", order = 2, param = param, about = mu)
ghypMom(10, param = param, momType = "mu")
momIntegrated("ghyp", order = 10, param = param, about = mu)

### raw moments
ghypMean(param = param)
ghypMom(1, param = param, momType = "raw")
momIntegrated("ghyp", order = 1, param = param, about = 0)
ghypMom(2, param = param, momType = "raw")
momIntegrated("ghyp", order = 2, param = param, about = 0)
ghypMom(10, param = param, momType = "raw")
momIntegrated("ghyp", order = 10, param = param, about = 0)

### central moments
ghypMom(1, param = param, momType = "central")
momIntegrated("ghyp", order = 1, param = param, about = m1)
ghypVar(param = param)
ghypMom(2, param = param, momType = "central")
momIntegrated("ghyp", order = 2, param = param, about = m1)
ghypMom(10, param = param, momType = "central")
momIntegrated("ghyp", order = 10, param = param, about = m1)



cleanEx()
nameEx("ghypParam")
### * ghypParam

flush(stderr()); flush(stdout())

### Name: ghypParam
### Title: Parameter Sets for the Generalized Hyperbolic Distribution
### Aliases: ghypParam ghypSmallShape ghypLargeShape ghypSmallParam
###   ghypLargeParam

### ** Examples

data(ghypParam)
plotShapeTriangle()
xis <- rep(c(0.1,0.3,0.5,0.7,0.9), 1:5)
chis <- c(0,-0.25,0.25,-0.45,0,0.45,-0.65,-0.3,0.3,0.65,
          -0.85,-0.4,0,0.4,0.85)
points(chis, xis, pch = 20, col = "red")


## Testing the accuracy of ghypMean
for (i in 1:nrow(ghypSmallParam)) {
  param <- ghypSmallParam[i, ]
  x <- rghyp(1000, param = param)
  sampleMean <- mean(x)
  funMean <- ghypMean(param = param)
  difference <- abs(sampleMean - funMean)
  print(difference)
}




cleanEx()
nameEx("gigCalcRange")
### * gigCalcRange

flush(stderr()); flush(stdout())

### Name: gigCalcRange
### Title: Range of a Generalized Inverse Gaussian Distribution
### Aliases: gigCalcRange
### Keywords: distribution

### ** Examples

param <- c(2.5, 0.5, 5)
maxDens <- dgig(gigMode(param = param), param = param)
gigRange <- gigCalcRange(param = param, tol = 10^(-3) * maxDens)
gigRange
curve(dgig(x, param = param), gigRange[1], gigRange[2])
## Not run: gigCalcRange(param = param, tol = 10^(-3), density = FALSE)



cleanEx()
nameEx("gigChangePars")
### * gigChangePars

flush(stderr()); flush(stdout())

### Name: gigChangePars
### Title: Change Parameterizations of the Generalized Inverse Gaussian
###   Distribution
### Aliases: gigChangePars
### Keywords: distribution

### ** Examples

param1 <- c(2.5, 0.5, 5)                # Parameterisation 1
param2 <- gigChangePars(1, 2, param1)   # Convert to parameterization 2
param2                                  # Parameterization 2
gigChangePars(2, 1, as.numeric(param2)) # Convert back to parameterization 1



cleanEx()
nameEx("gigCheckPars")
### * gigCheckPars

flush(stderr()); flush(stdout())

### Name: gigCheckPars
### Title: Check Parameters of the Generalized Inverse Gaussian
###   Distribution
### Aliases: gigCheckPars
### Keywords: distribution

### ** Examples

gigCheckPars(c(5, 2.5, -0.5))      # normal
gigCheckPars(c(-5, 2.5, 0.5))      # error
gigCheckPars(c(5, -2.5, 0.5))      # error
gigCheckPars(c(-5, -2.5, 0.5))     # error
gigCheckPars(c(0, 2.5, 0.5))       # gamma
gigCheckPars(c(0, 2.5, -0.5))      # error
gigCheckPars(c(0, 0, 0.5))         # error
gigCheckPars(c(0, 0, -0.5))        # error
gigCheckPars(c(5, 0, 0.5))         # error
gigCheckPars(c(5, 0, -0.5))        # invgamma



cleanEx()
nameEx("gigFit")
### * gigFit

flush(stderr()); flush(stdout())

### Name: gigFit
### Title: Fit the Generalized Inverse Gausssian Distribution to Data
### Aliases: gigFit print.gigFit plot.gigFit coef.gigFit vcov.gigFit
### Keywords: distribution

### ** Examples

param <- c(1, 1, 1)
dataVector <- rgig(500, param = param)
## See how well gigFit works
gigFit(dataVector)
##gigFit(dataVector, plots = TRUE)

## See how well gigFit works in the limiting cases
## Gamma case
dataVector2 <- rgamma(500, shape = 1, rate = 1)
gigFit(dataVector2)

## Inverse gamma
require(actuar)
dataVector3 <- rinvgamma(500, shape = 1, rate = 1)
gigFit(dataVector3)

## Use nlm instead of default
gigFit(dataVector, method = "nlm")




cleanEx()
nameEx("gigFitStart")
### * gigFitStart

flush(stderr()); flush(stdout())

### Name: gigFitStart
### Title: Find Starting Values for Fitting a Generalized Inverse Gaussian
###   Distribution
### Aliases: gigFitStart gigFitStartMoM gigFitStartLM
### Keywords: distribution

### ** Examples

param <- c(1, 1, 1)
dataVector <- rgig(500, param = param)
gigFitStart(dataVector)



cleanEx()
nameEx("gigHessian")
### * gigHessian

flush(stderr()); flush(stdout())

### Name: gigHessian
### Title: Calculate Two-Sided Hessian for the Generalized Inverse Gaussian
###   Distribution
### Aliases: gigHessian

### ** Examples

### Calculate the approximate Hessian using gigHessian:
param <- c(1,1,1)
dataVector <- rgig(500, param = param)
fit <- gigFit(dataVector)
coef <- coef(fit)
gigHessian(x = dataVector, param = coef, hessianMethod = "tsHessian",
              whichParam = 1)

### Or calculate the approximate Hessian using summary.gigFit method:
summary(fit, hessian = TRUE)



cleanEx()
nameEx("gigMeanVarMode")
### * gigMeanVarMode

flush(stderr()); flush(stdout())

### Name: Specific Generalized Inverse Gaussian Moments and Mode
### Title: Moments and Mode of the Generalized Inverse Gaussian
###   Distribution
### Aliases: gigMean gigVar gigSkew gigKurt gigMode
### Keywords: distribution

### ** Examples

param <- c(5, 2.5, -0.5)
gigMean(param = param)
gigVar(param = param)
gigSkew(param = param)
gigKurt(param = param)
gigMode(param = param)



cleanEx()
nameEx("gigMom")
### * gigMom

flush(stderr()); flush(stdout())

### Name: gigMom
### Title: Calculate Moments of the Generalized Inverse Gaussian
###   Distribution
### Aliases: gigMom gigRawMom gammaRawMom
### Keywords: distribution

### ** Examples

### Raw moments of the generalized inverse Gaussian distribution
param <- c(5, 2.5, -0.5)
gigRawMom(1, param = param)
momIntegrated("gig", order = 1, param = param, about = 0)
gigRawMom(2, param = param)
momIntegrated("gig", order = 2, param = param, about = 0)
gigRawMom(10, param = param)
momIntegrated("gig", order = 10, param = param, about = 0)
gigRawMom(2.5, param = param)

### Moments of the generalized inverse Gaussian distribution
param <- c(5, 2.5, -0.5)
(m1 <- gigRawMom(1, param = param))
gigMom(1, param = param)
gigMom(2, param = param, about = m1)
(m2 <- momIntegrated("gig", order = 2, param = param, about = m1))
gigMom(1, param = param, about = m1)
gigMom(3, param = param, about = m1)
momIntegrated("gig", order = 3, param = param, about = m1)

### Raw moments of the gamma distribution
shape <- 2
rate <- 3
param <- c(shape, rate)
gammaRawMom(1, shape, rate)
momIntegrated("gamma", order = 1, shape = shape, rate = rate, about = 0)
gammaRawMom(2, shape, rate)
momIntegrated("gamma", order = 2, shape = shape, rate = rate, about = 0)
gammaRawMom(10, shape, rate)
momIntegrated("gamma", order = 10, shape = shape, rate = rate, about = 0)

### Moments of the inverse gamma distribution
param <- c(5, 0, -0.5)
gigRawMom(2, param = param)             # Inf
gigRawMom(-2, param = param)
momIntegrated("invgamma", order = -2, shape = -param[3],
              rate = param[1]/2, about = 0)

### An example where the moment is infinite: inverse gamma
param <- c(5, 0, -0.5)
gigMom(1, param = param)
gigMom(2, param = param)



cleanEx()
nameEx("gigParam")
### * gigParam

flush(stderr()); flush(stdout())

### Name: gigParam
### Title: Parameter Sets for the Generalized Inverse Gaussian Distribution
### Aliases: gigParam gigSmallParam gigLargeParam

### ** Examples

data(gigParam)
## Check values of chi and psi
plot(gigLargeParam[, 1], gigLargeParam[, 2])
### Check all three parameters
pairs(gigLargeParam,
  labels = c(expression(chi),expression(psi),expression(lambda)))

## Testing the accuracy of gigMean
for (i in 1:nrow(gigSmallParam)) {
  param <- gigSmallParam[i, ]
  x <- rgig(1000, param = param)
  sampleMean <- mean(x)
  funMean <- gigMean(param = param)
  difference <- abs(sampleMean - funMean)
  print(difference)
}




cleanEx()
nameEx("hyperbCalcRange")
### * hyperbCalcRange

flush(stderr()); flush(stdout())

### Name: hyperbCalcRange
### Title: Range of a Hyperbolic Distribution
### Aliases: hyperbCalcRange
### Keywords: distribution

### ** Examples

par(mfrow = c(1, 2))
param <- c(0, 1, 3, 1)
hyperbRange <- hyperbCalcRange(param = param, tol = 10^(-3))
hyperbRange
curve(phyperb(x, param = param), hyperbRange[1], hyperbRange[2])
maxDens <- dhyperb(hyperbMode(param = param), param = param)
hyperbRange <- hyperbCalcRange(param = param, tol = 10^(-3) * maxDens, density = TRUE)
hyperbRange
curve(dhyperb(x, param = param), hyperbRange[1], hyperbRange[2])



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("hyperbChangePars")
### * hyperbChangePars

flush(stderr()); flush(stdout())

### Name: hyperbChangePars
### Title: Change Parameterizations of the Hyperbolic Distribution
### Aliases: hyperbChangePars
### Keywords: distribution

### ** Examples

param1 <- c(2, 1, 3, 1)                    # Parameterization 1
param2 <- hyperbChangePars(1, 2, param1)   # Convert to parameterization 2
param2                                     # Parameterization 2
hyperbChangePars(2, 1, param2)             # Back to parameterization 1



cleanEx()
nameEx("hyperbCvMTest")
### * hyperbCvMTest

flush(stderr()); flush(stdout())

### Name: hyperbCvMTest
### Title: Cramer-von~Mises Test of a Hyperbolic Distribution
### Aliases: hyperbCvMTest hyperbCvMTestPValue print.hyperbCvMTest
### Keywords: htest print

### ** Examples

param <- c(2, 2, 2, 1.5)
dataVector <- rhyperb(500, param = param)
fittedparam <- hyperbFit(dataVector)$param
hyperbCvMTest(dataVector, param = fittedparam)
dataVector <- rnorm(1000)
fittedparam <- hyperbFit(dataVector, startValues = "FN")$param
hyperbCvMTest(dataVector, param = fittedparam)



cleanEx()
nameEx("hyperbFit")
### * hyperbFit

flush(stderr()); flush(stdout())

### Name: hyperbFit
### Title: Fit the Hyperbolic Distribution to Data
### Aliases: hyperbFit print.hyperbFit plot.hyperbFit coef.hyperbFit
###   vcov.hyperbFit
### Keywords: distribution

### ** Examples

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




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("hyperbFitStart")
### * hyperbFitStart

flush(stderr()); flush(stdout())

### Name: hyperbFitStart
### Title: Find Starting Values for Fitting a Hyperbolic Distribution
### Aliases: hyperbFitStart hyperbFitStartMoM
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1)
dataVector <- rhyperb(500, param = param)
hyperbFitStart(dataVector, startValues = "FN")
hyperbFitStartMoM(dataVector)
hyperbFitStart(dataVector, startValues = "MoM")



cleanEx()
nameEx("hyperbHessian")
### * hyperbHessian

flush(stderr()); flush(stdout())

### Name: hyperbHessian
### Title: Calculate Two-Sided Hessian for the Hyperbolic Distribution
### Aliases: hyperbHessian sumX

### ** Examples

### Calculate the exact Hessian using hyperbHessian:
param <- c(2, 2, 2, 1)
dataVector <- rhyperb(500, param = param)
fit <- hyperbFit(dataVector, method = "BFGS")
coef <- coef(fit)
hyperbHessian(x = dataVector, param = coef, hessianMethod = "exact",
              whichParam = 2)
              
### Or calculate the exact Hessian using summary.hyperbFit method:
summary(fit, hessian = TRUE)


## Calculate the approximate Hessian:
summary(fit, hessian = TRUE, hessianMethod = "tsHessian")



cleanEx()
nameEx("hyperbMeanVarMode")
### * hyperbMeanVarMode

flush(stderr()); flush(stdout())

### Name: Specific Hyperbolic Distribution Moments and Mode
### Title: Moments and Mode of the Hyperbolic Distribution
### Aliases: hyperbMean hyperbVar hyperbSkew hyperbKurt hyperbMode
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1)
hyperbMean(param = param)
hyperbVar(param = param)
hyperbSkew(param = param)
hyperbKurt(param = param)
hyperbMode(param = param)



cleanEx()
nameEx("hyperbParam")
### * hyperbParam

flush(stderr()); flush(stdout())

### Name: hyperbParam
### Title: Parameter Sets for the Hyperbolic Distribution
### Aliases: hyperbParam hyperbSmallShape hyperbLargeShape hyperbSmallParam
###   hyperbLargeParam

### ** Examples

data(hyperbParam)
plotShapeTriangle()
xis <- rep(c(0.1,0.3,0.5,0.7,0.9), 1:5)
chis <- c(0,-0.25,0.25,-0.45,0,0.45,-0.65,-0.3,0.3,0.65,
          -0.85,-0.4,0,0.4,0.85)
points(chis, xis, pch = 20, col = "red")


## Testing the accuracy of hyperbMean
for (i in 1:nrow(hyperbSmallParam)) {
  param <- hyperbSmallParam[i, ]
  x <- rhyperb(1000, param = param)
  sampleMean <- mean(x)
  funMean <- hyperbMean(param = param)
  difference <- abs(sampleMean - funMean)
  print(difference)
}




cleanEx()
nameEx("mamquam")
### * mamquam

flush(stderr()); flush(stdout())

### Name: mamquam
### Title: Size of Gravels from Mamquam River
### Aliases: mamquam
### Keywords: datasets

### ** Examples

data(mamquam)
str(mamquam)
### Construct data from frequency summary, taking all observations
### at midpoints of intervals
psi <- rep(mamquam$midpoints, mamquam$counts)
barplot(table(psi))
### Fit the hyperbolic distribution
hyperbFit(psi)

### Actually hyperbFit can deal with frequency data
hyperbFit(mamquam$midpoints, freq = mamquam$counts)



cleanEx()
nameEx("momRecursion")
### * momRecursion

flush(stderr()); flush(stdout())

### Name: momRecursion
### Title: Computes the moment coefficients recursively for generalized
###   hyperbolic and related distributions
### Aliases: momRecursion
### Keywords: distribution

### ** Examples

  momRecursion(order = 12)

  #print out the matrix
  momRecursion(order = 12, "true")



cleanEx()
nameEx("nervePulse")
### * nervePulse

flush(stderr()); flush(stdout())

### Name: nervePulse
### Title: Intervals Between Pulses Along a Nerve Fibre
### Aliases: nervePulse
### Keywords: datasets

### ** Examples

data(nervePulse)
str(nervePulse)

### Fit the generalized inverse Gaussian distribution
gigFit(nervePulse)




cleanEx()
nameEx("nigCalcRange")
### * nigCalcRange

flush(stderr()); flush(stdout())

### Name: nigCalcRange
### Title: Range of a normal inverse Gaussian Distribution
### Aliases: nigCalcRange
### Keywords: distribution

### ** Examples

par(mfrow = c(1, 2))
param <- c(0, 1, 3, 1)
nigRange <- nigCalcRange(param = param, tol = 10^(-3))
nigRange
curve(pnig(x, param = param), nigRange[1], nigRange[2])
maxDens <- dnig(nigMode(param = param), param = param)
nigRange <- nigCalcRange(param = param, tol = 10^(-3) * maxDens, density = TRUE)
nigRange
curve(dnig(x, param = param), nigRange[1], nigRange[2])



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("nigFit")
### * nigFit

flush(stderr()); flush(stdout())

### Name: nigFit
### Title: Fit the normal inverse Gaussian Distribution to Data
### Aliases: nigFit print.nigFit plot.nigFit coef.nigFit vcov.nigFit
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1)
dataVector <- rnig(500, param = param)
## See how well nigFit works
nigFit(dataVector)
nigFit(dataVector, plots = TRUE)
fit <- nigFit(dataVector)
par(mfrow = c(1, 2))
plot(fit, which = c(1, 3))

## Use nlm instead of default
nigFit(dataVector, method = "nlm")




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("nigFitStart")
### * nigFitStart

flush(stderr()); flush(stdout())

### Name: nigFitStart
### Title: Find Starting Values for Fitting a normal inverse Gaussian
###   Distribution
### Aliases: nigFitStart nigFitStartMoM
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1)
dataVector <- rnig(500, param = param)
nigFitStart(dataVector, startValues = "FN")
nigFitStartMoM(dataVector)
nigFitStart(dataVector, startValues = "MoM")



cleanEx()
nameEx("nigHessian")
### * nigHessian

flush(stderr()); flush(stdout())

### Name: nigHessian
### Title: Calculate Two-Sided Hessian for the Normal Inverse Gaussian
###   Distribution
### Aliases: nigHessian

### ** Examples

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



cleanEx()
nameEx("nigMeanVarMode")
### * nigMeanVarMode

flush(stderr()); flush(stdout())

### Name: Specific Normal Inverse Gaussian Distribution Moments and Mode
### Title: Moments and Mode of the Normal Inverse Gaussian Distribution
### Aliases: nigMean nigVar nigSkew nigKurt nigMode
### Keywords: distribution

### ** Examples

param <- c(2, 2, 2, 1)
nigMean(param = param)
nigVar(param = param)
nigSkew(param = param)
nigKurt(param = param)
nigMode(param = param)



cleanEx()
nameEx("nigParam")
### * nigParam

flush(stderr()); flush(stdout())

### Name: nigParam
### Title: Parameter Sets for the Normal Inverse Gaussian Distribution
### Aliases: nigParam nigSmallShape nigLargeShape nigSmallParam
###   nigLargeParam

### ** Examples

data(nigParam)
plotShapeTriangle()
xis <- rep(c(0.1,0.3,0.5,0.7,0.9), 1:5)
chis <- c(0,-0.25,0.25,-0.45,0,0.45,-0.65,-0.3,0.3,0.65,
          -0.85,-0.4,0,0.4,0.85)
points(chis, xis, pch = 20, col = "red")


## Testing the accuracy of nigMean
for (i in 1:nrow(nigSmallParam)) {
  param <- nigSmallParam[i, ]
  x <- rnig(1000, param = param)
  sampleMean <- mean(x)
  funMean <- nigMean(param = param)
  difference <- abs(sampleMean - funMean)
  print(difference)
}




cleanEx()
nameEx("plotShapeTriangle")
### * plotShapeTriangle

flush(stderr()); flush(stdout())

### Name: plotShapeTriangle
### Title: Plot the Shape Triangle
### Aliases: plotShapeTriangle
### Keywords: distribution

### ** Examples

plotShapeTriangle()



cleanEx()
nameEx("qqghyp")
### * qqghyp

flush(stderr()); flush(stdout())

### Name: GeneralizedHyperbolicPlots
### Title: Generalized Hyperbolic Quantile-Quantile and Percent-Percent
###   Plots
### Aliases: qqghyp ppghyp
### Keywords: hplot distribution

### ** Examples

par(mfrow = c(1, 2))
y <- rghyp(200, param = c(2, 2, 2, 1, 2))
qqghyp(y, param = c(2, 2, 2, 1, 2), line = FALSE)
abline(0, 1, col = 2)
ppghyp(y, param = c(2, 2, 2, 1, 2))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("qqgig")
### * qqgig

flush(stderr()); flush(stdout())

### Name: GIGPlots
### Title: Generalized Inverse Gaussian Quantile-Quantile and
###   Percent-Percent Plots
### Aliases: qqgig ppgig
### Keywords: hplot distribution

### ** Examples

par(mfrow = c(1, 2))
y <- rgig(1000, param = c(2, 3, 1))
qqgig(y, param = c(2, 3, 1), line = FALSE)
abline(0, 1, col = 2)
ppgig(y, param = c(2, 3, 1))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("qqhyperb")
### * qqhyperb

flush(stderr()); flush(stdout())

### Name: HyperbPlots
### Title: Hyperbolic Quantile-Quantile and Percent-Percent Plots
### Aliases: qqhyperb pphyperb
### Keywords: hplot distribution

### ** Examples

par(mfrow = c(1, 2))
param <- c(2, 2, 2, 1.5)
y <- rhyperb(200, param = param)
qqhyperb(y, param = param, line = FALSE)
abline(0, 1, col = 2)
pphyperb(y, param = param)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("qqnig")
### * qqnig

flush(stderr()); flush(stdout())

### Name: nigPlots
### Title: Normal inverse Gaussian Quantile-Quantile and Percent-Percent
###   Plots
### Aliases: qqnig ppnig
### Keywords: hplot distribution

### ** Examples

par(mfrow = c(1, 2))
param <- c(2, 2, 2, 1.5)
y <- rnig(200, param = param)
qqnig(y, param = param, line = FALSE)
abline(0, 1, col = 2)
ppnig(y, param = param)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("qqskewlap")
### * qqskewlap

flush(stderr()); flush(stdout())

### Name: SkewLaplacePlots
### Title: Skew-Laplace Quantile-Quantile and Percent-Percent Plots
### Aliases: qqskewlap ppskewlap
### Keywords: hplot distribution

### ** Examples

par(mfrow = c(1, 2))
y <- rskewlap(1000, param = c(2, 0.5, 1))
qqskewlap(y, param = c(2, 0.5, 1), line = FALSE)
abline(0, 1, col = 2)
ppskewlap(y, param = c(2, 0.5, 1))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("resistors")
### * resistors

flush(stderr()); flush(stdout())

### Name: resistors
### Title: Resistance of One-half-ohm Resistors
### Aliases: resistors
### Keywords: datasets

### ** Examples

data(resistors)
str(resistors)
### Construct data from frequency summary, taking all observations
### at midpoints of intervals
resistances <- rep(resistors$midpoints, resistors$counts)
hist(resistances)
logHist(resistances)
## Fit the hyperbolic distribution
hyperbFit(resistances) 

## Actually fit.hyperb can deal with frequency data
hyperbFit(resistors$midpoints, freq = resistors$counts)




cleanEx()
nameEx("summary.gigFit")
### * summary.gigFit

flush(stderr()); flush(stdout())

### Name: summary.gigFit
### Title: Summarizing Normal Inverse Gaussian Distribution Fit
### Aliases: summary.gigFit print.summary.gigFit
### Keywords: distribution

### ** Examples

### Continuing the  gigFit(.) example:
param <- c(1,1,1)
dataVector <- rgig(500, param = param)
fit <- gigFit(dataVector)
print(fit)
summary(fit, hessian = TRUE, hessianMethod = "tsHessian")



cleanEx()
nameEx("summary.hyperbFit")
### * summary.hyperbFit

flush(stderr()); flush(stdout())

### Name: summary.hyperbFit
### Title: Summarizing Hyperbolic Distribution Fit
### Aliases: summary.hyperbFit print.summary.hyperbFit
### Keywords: distribution

### ** Examples

### Continuing the  hyperbFit(.) example:
param <- c(2, 2, 2, 1)
dataVector <- rhyperb(500, param = param)
fit <- hyperbFit(dataVector, method = "BFGS")
print(fit)
summary(fit, hessian = TRUE)



cleanEx()
nameEx("summary.nigFit")
### * summary.nigFit

flush(stderr()); flush(stdout())

### Name: summary.nigFit
### Title: Summarizing Normal Inverse Gaussian Distribution Fit
### Aliases: summary.nigFit print.summary.nigFit
### Keywords: distribution

### ** Examples

### Continuing the  nigFit(.) example:
param <- c(2, 2, 2, 1)
dataVector <- rnig(500, param = param)
fit <- nigFit(dataVector, method = "BFGS")
print(fit)
summary(fit, hessian = TRUE, hessianMethod = "tsHessian")



cleanEx()
nameEx("traffic")
### * traffic

flush(stderr()); flush(stdout())

### Name: traffic
### Title: Intervals Between Vehicles on a Road
### Aliases: traffic
### Keywords: datasets

### ** Examples

data(traffic)
str(traffic)

### Fit the generalized inverse Gaussian distribution
gigFit(traffic)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
