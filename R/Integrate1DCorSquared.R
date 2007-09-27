`Integrate1DCorSquared` <-
function(x, xKnown, betaj, aj.power, corMatrix, W) {
        origCorMatrix = corMatrix
        y = x
        for (i in 1:length(x)) {
                corMatrix = origCorMatrix * exp(-betaj*abs(x[i]-xKnown)**aj.power)
                y[i] = (corMatrix%*%W)**2
        }
        return (y)
}

