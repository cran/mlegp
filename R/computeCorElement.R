`computeCorElement` <-
function(beta, a, theta, newTheta, holdIndex, 
	lower = rep(0,length(theta)), upper = rep(1, length(theta))) {
        ans = 1
        for (i in 1:length(theta)) {
                if (i %in% holdIndex) {
                        ans = ans * exp(-beta[i]*abs(theta[i]-newTheta[i])**a[i])
                }
                else {
                        integratedValue = integrate(Integrate1DCor,lower[i], upper[i],xknown = theta[i], 
				beta = beta[i], a.power = a[i])
                        ans = ans * integratedValue$value * 1 / (upper[i]-lower[i])
                }
        }
        return (ans)
}

