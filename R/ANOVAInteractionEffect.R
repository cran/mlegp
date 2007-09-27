`ANOVAInteractionEffect` <-
function(gp, effectNums, mu0 = NULL, lower = apply(gp$X,2,min), upper = apply(gp$X,2,max), maxpts = NULL) {

  if (length(effectNums) != 2) {
          stop("effectsNum must be of length 2 for ANOVAInteractionEffect")
          return (NULL)
  }

  if (is.null(mu0)) {
          mu0 = OverallMean(gp, lower = lower, upper = upper)
  }

 meanRegNoX1 = expectedMeanRegExcept(gp, effectNums[1], lower = lower, upper = upper)
 meanRegNoX2 = expectedMeanRegExcept(gp, effectNums[2], lower = lower, upper = upper)
 meanRegNoX1X2 = expectedMeanRegExcept(gp, effectNums, lower = lower, upper = upper)

  corMatrixNoX1 = createIntegratedCorMatrixExcept(gp$X, gp$beta, gp$a, effectNums[1], lower=lower,upper=upper)
  corMatrixNoX2 = createIntegratedCorMatrixExcept(gp$X, gp$beta, gp$a, effectNums[2], lower=lower,upper=upper)
  corMatrixNoX1X2 = createIntegratedCorMatrixExcept(gp$X, gp$beta, gp$a, effectNums, lower=lower,upper=upper)

  W = gp$invVarMatrix%*%(gp$Z - gp$mu)

  integratedValue=adapt(2, lower[effectNums],upper[effectNums], maxpts = maxpts, functn=IntegrateInteractionEffectSquared,
                        gp = gp, 
			meanRegNoX1 = meanRegNoX1, meanRegNoX2 = meanRegNoX2, meanRegNoX1X2 = meanRegNoX1X2,
                        corMatrixNoX1 = corMatrixNoX1, corMatrixNoX2 = corMatrixNoX2,
                        corMatrixNoX1X2 = corMatrixNoX1X2, W = W, mu0 = mu0, effectNums = effectNums)
   return (integratedValue$value * 1 / prod(upper[effectNums]-lower[effectNums]))
}

