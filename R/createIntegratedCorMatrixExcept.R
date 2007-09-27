`createIntegratedCorMatrixExcept` <-
function(X, beta, a, except, lower = rep(0,dim(X)[2]), upper = rep(1,dim(X)[2])) {

        corNewObsMatrix = matrix(0, ncol = dim(X)[1])
        for (p in 1:dim(X)[1]) {
                corNewObsMatrix[p] = computeCorElementExcept(beta, a, as.double(X[p,]), except, lower, upper)
        }

        return (corNewObsMatrix)
}

