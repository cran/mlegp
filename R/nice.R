`nice` <-
function(x, digits = 4) {
        for (i in 2:dim(x)[2]) {
                x[,i] = round(x[,i],digits)
        }
        return (x)
}

