#' Negative Log likelihood function for dirichlet data
#' @param x data for loglikelihood
#' @param alp estimated values for the dirichlet parameters
llk <- function(x, alp) {
    nrow(x) * (log(gamma(sum(alp))) - sum(log(gamma(alp))) +
                   sum((alp - 1) * colMeans(log(x))))
}
#' stabilization
#' @param x data for stabilization.
#' @param eps epsion for stabilization. It is very small number to extremely small data points
stbz <- function(x, eps = 1e-10) {
  (x + eps) / (1 + 2 * eps)
}

#' Get MME for dirichlet
#' @param x data for estimating MME.
Diri_MME <- function(x) {  # calculate MME using data
    m1 <- colMeans(x)
    m2 <- colMeans(x ^ 2)
    return(matrix((m1 * (m1 - m2)) / (m2 - m1 ^ 2) , nrow = 1))
}

#' Get MLE for dirichlet
#' @param x data for estimating MLEce
#' @param eps epsilon for iterative algorithm.
#' @param mxit maximum iteration for MLE. Default is 1e5
#' @import sirt
Diri_MLE <- function(x, eps = 1e-10, mxit = 1e5) {  # calculate ML estimator using data
    return(matrix(dirichlet.mle(x, eps = eps, maxit = mxit)$alpha ,nrow = 1))  # outputs estimated parameters in a row vector
}

#' Get MLEce for dirichlet
#' @param x data for estimating MLEce
Diri_CE <- function(x) {  # calculate closed-form estimator using data
    m <- ncol(x)
    MME <- Diri_MME(x)
    c1 <- trigamma(c(sum(MME), MME))
    denom <- 1 / c1[1] - sum(1 / c1[-1])
    c2 <- digamma(c(sum(MME), MME))
    iD <- diag(1 / c1[-1])
    iH <-
        iD %*% (diag(1, m) + matrix(1, nrow = m, ncol = m) %*% iD / denom)
    grad <- colMeans(log(x)) - (c2[-1] - c2[1])
    ans <- matrix(MME + iH %*% grad, nrow = 1)
    otp = pmax(ans, 1 / nrow(x))
    estls = list(estimation = otp)
    return(estls)
}

#' beta tilde: root-n consistent estimator
#' @param x data for estimating MLEce
Diri_CE_bt <- function(x) {
    m <- ncol(x)
    tldb <- log(Diri_MME(x))

    # constants
    c1 <- trigamma(c(sum(exp(tldb)), exp(tldb)))
    c2 <- digamma(c(sum(exp(tldb)), exp(tldb)))

    D1ll <-
        exp(tldb) * (colSums(log(x)) - nrow(x) * (c2[-1] - c2[1]))  # Gradient vector
    D2ll <- function(i, j) {
        # Hessian function
        if (i == j) {
            # diagonal elements
            D1ll[j] + nrow(x) * exp(tldb)[j] ^ 2 * (c1[1] - c1[-1][j])
        }
        else {
            # o.w.
            nrow(x) * exp(tldb)[i] * exp(tldb)[j] * c1[1]
        }
    }
    D2ll <- Vectorize(D2ll)
    iH <- solve(outer(1:m, 1:m, D2ll))   # inverse Hessian
    ans <- tldb - D1ll %*% iH  # estimation
    return(exp(ans))  # back-transformation
}

#' SD2 statistics for GCVM gof test
#' @param Data data for gof test
#' @param EST estimates for gof test
SD2fun.diri <- function(Data, EST){
  N = dim(Data)[1]
  S = dim(Data)[2]

  #Rosenblatt Transformation
  TransData = matrix(NA, N,S)

  #marginal distribution of #X1~Beta(alpha_1,sum(alpha))
  TransData[,1] = pbeta(Data[,1], EST[1], sum(EST))
  #conditional dist:1/(1-x1) * X2|X1~Beta(alpha_2,sum(alpha) - alpha_1)
  for (i in 2:S) {
    TransData[,i] = pbeta(Data[,i], EST[i], sum(EST) - sum(EST[1:i-1])) * (1 - rowSums(as.matrix(Data[,1:i-1])))
  }

  SD2par1 = (4/3)^S
  SD2par2.func = function(x) {1 + 2*x - 2*x^2}
  caled = apply(TransData, 2, SD2par2.func)
  SD2par2 <- 2*sum(apply(caled, 1, prod) )/N

  grd = expand.grid(1:N, 1:N)
  xx1 = TransData[grd[,1], ]
  xx2 = TransData[grd[,2], ]

  ## init
  SD2par3 = 1
  for (i in 1:S){
    SD2par3 = SD2par3 * (1-abs(xx1[,i]-xx2[,i]))
  }


  SD2 <- sqrt(SD2par1-SD2par2 + (2^S)*sum(SD2par3)/(N^2) )

  return(Statis=SD2)
}
