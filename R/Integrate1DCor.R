`Integrate1DCor` <-
function(x, xknown, beta, a.power) {
        return (exp(-beta*abs(x-xknown)**a.power))
}

