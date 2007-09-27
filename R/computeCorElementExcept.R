`computeCorElementExcept` <-
function(beta, a, theta, except, lower = rep(0,length(theta)), upper = rep(1,length(theta))) {
        #file = "integrate.out"
        ans = 1
        for (i in 1:length(theta)) {
                if (i %in% except) {
                        #do nothing!
                }
                else {
                        integratedValue = integrate(Integrate1DCor,lower[i],upper[i],xknown = theta[i], 
				beta = beta[i], a.power = a[i])
                        ans = ans * integratedValue$value * 1 / (upper[i]-lower[i])
                }
        }
        return (ans)
}

