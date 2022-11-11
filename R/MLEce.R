#' Calculating a value of MLEce according to a distribution
#'
#' Provide a function that numerically computes the closed-form estimator which is asymptotically efficient and
#' can be computed faster than the maximum likelihood estimator.
#'
#' @details
#' The closed-form estimation procedure is based on root n-consistent estimators and a Fisher scoring step or a Newton step on the loglikelihood function.
#' The estimator is obtained by solving the linear equation By E.L. Lehmann.
#' This estimator follows the multivariate normal distribution with mean vector of 0 and variance matrix of inverse of Fisher Information matrix and
#' has the properties of a multivariate normal distribution.
#'
#' @param data a numeric vector or matrix.
#' @param distrib a character string "name" naming a distribution for which the corresponding density function dname
#' and the corresponding distribution function pname must be classically defined.
#' @param boots a number of iteration for calculating CI-parametric intervals.
#' @param CI.alpha a significance level of confidence intervals. default is 0.05.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "MLEce".
#' @param digits a numeric number of significant digits.
#'
#' @return an object of class “MLEce”. It is a list with the following components:
#' \item{estimation}{the parameter estimates.}
#' \item{distribution }{a character string of a distribution assuming that data set comes from.}
#' \item{stat_summary}{a numeric vector with min, 1st quantile, median, 3rd quantile, and max.}
#' \item{CI}{a matrix with confidence intervals of Estimates obtained by CI-parametric bootstrapping.}
#' \item{n}{a numeric value of data length.}
#' \item{data}{the data set.}
#' @examples
#' datt = rBiGam(100, c(1,4,5))
#' res = MLEce(datt, "BiGam", boots = 50)
#' @import stats
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom graphics title contour par points
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics legend rug
#' @export
MLEce = function (data, distrib, boots = 1000, CI.alpha = 0.05, ...) {
    if (!is.character(distrib)) {
        stop("distr must be a character string naming a distribution")
    }
    else {
        distname <- distrib
    }

    distlist <- c("BiWei", "BiGam", "Dirichlet")
    if (! distname %in% distlist){
        stop(paste("The ", distname, " function must be defined"))
    }

    if (!(is.numeric(as.matrix(data)) & length(data) > 1)) {
        stop("data must be a numeric vector of length greater than 1")
    }

    if (distname %in% c("BiWei", "BiGam", "Dirichlet")){
        variate <- "Bi"
    }
    else variate <- "Uni"


    obs <- data

    if (variate == "Uni"){
        n <- length(obs)
        stat_summary <- quantile(obs)
        meanx <- mean(obs); varx <- var(obs)
    }
    if (variate == "Bi"){
        n <- dim(obs)[1]
        stat_summary <- t(apply(obs, 2, quantile))
    }


    res <- MLEce_est(data = obs,distname = distname, boots = boots, CI.alpha = CI.alpha)
    resls <- append(res,
                    list(distribution = distname,
                         stat_summary = stat_summary,
                         n = n, data = obs,
                         CI.alpha = CI.alpha )
                    )

    class(resls) <- c("MLEce")
    return(resls)
}

#' Estimate MLEce
#' @param data a numeric vector or matrix.
#' @param distname a character string "name" naming a distribution for which the corresponding density function dname
#' and the corresponding distribution function pname must be classically defined.
#' @param boots a number of iteration for calculating CI-parametric intervals.
#' @param CI.alpha a significance level of confidence intervals. default is 0.05.
MLEce_est = function(data, distname, boots, CI.alpha){
    distls <- c("BiWei", "BiGam", "Dirichlet")

    if (!is.element(distname, distls)) {
        stop("unsupported distribution")

    }
    else {
        initb <- FALSE
    }

    if (distname == "BiWei") {
        est_ls = BiWei_CE(data)
        ce_est = est_ls$estimation
        hess = BiWei_info(ce_est, data, type = "hessian")$hessian
        n = dim(data)[1]

        asym.var = -solve(hess)
        par.se = sqrt(diag(asym.var)/n)

        ce_pb = matrix(0,nrow=boots,ncol=5)
        for (i in 1:boots){
            rg_dat = rBiWei(n, ce_est)
            ce_pb[i,] = BiWei_CE(rg_dat)$estimation
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)

        res = append(est_ls,
                     list(CI = CI, se = par.se))
    }

    if (distname == "BiGam") {
        mme_gam = BiGam_MME(data)
        estls = BiGam_CE(mme_gam, data)
        ce_est = estls$estimation
        hess = estls$hess
        n = dim(data)[1]

        testStat <- c()
        for(i in 1:1000){
          rg_dat = rBiGam(n, ce_est)
          testStat[i] <- SD2fun.gam(rg_dat, ce_est)
        }
        SD2value <- SD2fun.gam(data, ce_est)
        p.GOF <- mean(SD2value>=testStat)

        asym.var = -solve(hess)
        par.se = sqrt(diag(asym.var)/n)

        ce_pb = matrix(0,nrow=boots,ncol=3)
        for (i in 1:boots){
            rg_dat = rBiGam(n, ce_est)
            ce_pb[i,] = BiGam_CE(BiGam_MME(rg_dat),rg_dat)$estimation
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)
        res = list(estimation = ce_est, p.GOF = p.GOF, CI = CI, se = par.se)

    }

    if (distname == "Dirichlet"){
        data = stbz(data)
        ce_est = Diri_CE_bt(data)
        n = dim(data)[1]

        testStat <- c()
        for(i in 1:1000){
          rg_dat = stbz(rdirichlet(n, ce_est))
          testStat[i] = SD2fun.diri(rg_dat, ce_est)
        }
        SD2value = SD2fun.diri(data, ce_est)
        p.GOF = mean(SD2value>=testStat)

        ce_pb = matrix(0,nrow=boots, ncol = length(ce_est))
        for (i in 1:boots){
            rg_dat = stbz(rdirichlet(n, ce_est))
            ce_pb[i,] = Diri_CE_bt(rg_dat)
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)

        res = list(estimation = ce_est, p.GOF = p.GOF, CI = CI)
    }
    return(res)
}

#' Statistical inference of MLE
#'
#' Provide value of maximum likelihood estimator, a result of GOF test, and information on the CI of MLE.
#' @param data a numeric vector or matrix.
#' @param distname a character string of a distribution assuming that data set comes from (Currently, only Biweibull distribution can be input)
#' @param inits a initial vector for MLE.
#' @param boots a number of iteration for parametric bootstrapping compute confidence intervals.
#' @param CI.alpha a significance level of confidence intervals. default is 0.05.
#'
#' @return \code{MLE_est} returns a list include maximum likelihood estimators and confidence interval of estimated parameters.
#' @import stats sirt
#' @export
MLE_est = function(data, distname, inits, boots = 1000, CI.alpha = 0.05){
    distls <- c("BiWei", "BiGam", "Dirichlet")

    if (!is.element(distname, distls)) {
        stop("unsupported distribution")
    }
    else {
        initb <- FALSE
    }

    if (distname == "BiWei") {
        inits = RNCE_est(data, distname)
        ml_est = BiWei_ML(data, inits = inits)
        n = dim(data)[1]
        hess = BiWei_info(ml_est, data, type = "hessian")$hessian

        asym.var = -solve(hess)
        par.se = sqrt(diag(asym.var)/n)

        teststat_rosen=0
        for(i in 1:999){
            rg_dat = rBiWei(n, ml_est)
            tmprosen = Rosen(rg_dat,ml_est)
            teststat_rosen[i] = sum(CD2(tmprosen[[1]]),CD2(tmprosen[[2]]))
        }
        orirosen = Rosen(rg_dat,ml_est)
        p.GOF <- (sum(teststat_rosen>sum(CD2(orirosen[[1]]),CD2(orirosen[[2]])))+1)/1000

        ce_pb = matrix(0,nrow=boots,ncol=5)
        for (i in 1:boots){
            rg_dat = rBiWei(n,ml_est)
            ce_pb[i,] = BiWei_ML(rg_dat, inits = inits)
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)



        res = list(estimation.ml = ml_est, p.GOF.ml = p.GOF, CI.ml = CI, se.ml = par.se)
    }

    if (distname == "BiGam") {
        mme_gam = BiGam_MME(data)
        ml_est = BiGam_MLE(mme_gam, data)
        n = dim(data)[1]

        testStat <- c()
        for(i in 1:1000){
          rg_dat = rBiGam(n, ml_est)
          testStat[i] <- SD2fun.gam(rg_dat, ml_est)
        }
        SD2value <- SD2fun.gam(data, ml_est)
        p.GOF.ml <- mean(SD2value>=testStat)

        res = list(estimation.ml = ml_est, p.GOF.ml = p.GOF)
    }

    if (distname == "Dirichlet"){
        ml_est = Diri_MLE(data, eps = 1e-10)
        n = dim(data)[1]

        testStat <- c()
        for(i in 1:1000){
          rg_dat = rdirichlet(n, ml_est)
          testStat[i] = SD2fun.diri(rg_dat, ml_est)
        }
        SD2value = SD2fun.diri(data, ml_est)
        p.GOF.ml = mean(SD2value>=testStat)

        res = list(estimation.ml = ml_est, p.GOF.ml = p.GOF)

    }

    return(res)
}



#' Statistical inference of root-n consistent estimator
#'
#' Provide value of root-n consistent estimator
#' @param data a numeric vector or matrix.
#' @param distname a character string of a distribution assuming that data set
#'
#' @return a numeric vector of estimated parameters.
#' @import stats
#' @export
RNCE_est = function(data, distname){
    distls <- c("BiWei", "BiGam", "Dirichlet")

    if (!is.element(distname, distls)) {
        stop("unsupported distribution")
    }
    else {
        initb <- FALSE
    }

    if (distname == "BiWei") {
        marginal_cme = c(cor_method(data[,1]),cor_method(data[,2]))
        delta_mle = pnorm(nleqslv::nleqslv(0,function(delta){delta_score_probit(c(marginal_cme,delta),data)})$x)

        res = c(marginal_cme,delta_mle)
        names(res) <- c('alpha1', 'beta1', 'alpha2', 'beta2', 'delta')
    }

    if (distname == "BiGam") {
        res = BiGam_MME(data)
    }

    if (distname == "Dirichlet") {
        m1 <- colMeans(data)
        m2 <- colMeans(data ^ 2)
        res = matrix((m1 * (m1 - m2)) / (m2 - m1 ^ 2) ,nrow = 1)
    }

    return(res)

}

#' @rdname MLEce
#' @method print MLEce
#' @export
print.MLEce <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nEstimate:\n")
    if(x$distribution == "BiWei"){
        est = round(x$estimation, digits = digits)
        cat("alpha1 = ", est[1], ", ", "beta1 = ", est[2])
        cat("\n")
        cat("alpha2 = ", est[3], ", ", "beta2 = ", est[4])
        cat("\n")
        cat("delta = ", est[5])
    }

    if (x$distribution == "BiGam"){
        est = round(x$estimation, digits = digits)
        cat("shape1 = ", est[1], ", ", "shape2 = ", est[2], ", ", "scale = ", est[3])
    }

    if (x$distribution == "Dirichlet"){
        est = round(x$estimation, digits = digits)
        cat("alpha1 = ", est[1], ", ", "alpha2 = ", est[2], ", ", "alpha3 = ", est[3])
    }
}

#' Extracting estimates.
#' \code{coef} extracts estimated parameters.
#' @param object an object of class "MLEce" made by the function \code{MLEce}.
#' @param digits a numeric number of significant digits.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "MLEce".
#' @return estimated parameters are extracted from the "MLEce" class.
#' @method coef MLEce
#' @export
coef.MLEce <- function(object, digits = max(3, getOption("digits") - 3), ...){
  est = round(object$estimation, digits = digits)
  result = list(coef.out = est)
  class(result) = c('coef.MLEce')
  return(result)
}

#' @rdname coef.MLEce
#' @method print coef.MLEce
#' @export
print.coef.MLEce <- function(x, digits = max(3, getOption("digits") - 3), ...){
  print(x$coef.out)
}

#' test goodness of fit
#' @examples
#' datt = rBiGam(100, c(1,4,5))
#' res = MLEce(datt, "BiGam", boots = 50)
#' gof(res)
#' @param x an object of class "MLEce" made by the function \code{MLEce}.
#' @param digits a numeric number of significant digits.
#' @param ... not used, but exists because of the compatibility.
#' @details
#' The null hypothesis of the GCVM test is that "data follows the Bivariate Weibull distribution".
#' @return \code{gof} returns the p-value of the GCVM test.
#' @export
gof <- function(x, digits = max(3, getOption("digits") - 3), ...){
    if (!inherits(x,"MLEce")) {
        message("Error: It it not class 'MLEce'.")
    }

    if (x$distribution == "BiWei"){
      ce_est = x$estimation
      n = dim(x$data)[1]
      teststat_rosen=0
      for(i in 1:999){
        rg_dat = rBiWei(n, ce_est)
        tmprosen = Rosen(rg_dat,ce_est)
        teststat_rosen[i] = sum(CD2(tmprosen[[1]]),CD2(tmprosen[[2]]))
      }
      orirosen = Rosen(x$data,ce_est)
      p.GOF <- (sum(teststat_rosen>sum(CD2(orirosen[[1]]),CD2(orirosen[[2]])))+1)/1000
    }
    if (x$distribution == "BiGam"){
      ce_est = x$estimation
      n = dim(x$data)[1]

      ## GOF test
      testStat <- c()
      for(i in 1:1000){
        rg_dat = rBiGam(n, ce_est)
        testStat[i] <- SD2fun.gam(rg_dat, ce_est)
      }
      SD2value <- SD2fun.gam(x$data, ce_est)
      p.GOF <- mean(SD2value>=testStat)

    }
    if (x$distribution == "Dirichlet"){
      data = stbz(x$data)
      ce_est = x$estimation
      n = dim(data)[1]

      ## GOF test
      testStat <- c()
      for(i in 1:1000){
        rg_dat = stbz(rdirichlet(n, ce_est))
        testStat[i] = SD2fun.diri(rg_dat, ce_est)
      }
      SD2value = SD2fun.diri(data, ce_est)
      p.GOF = mean(SD2value>=testStat)
    }
  result = list(p.GOF = p.GOF, distribution = x$distribution)
  class(result) <- "gof"
  return(result)
}

#' @rdname gof
#' @method print gof
#' @export
print.gof <- function(x, digits = max(3, getOption("digits") - 3), ...){
  if (x$distribution == "BiWei"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-centred discrepancy (CD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the ", x$distribution, " distribution\n")
    cat(paste0('\np.value : ', x$p.GOF, '\n'))
  }
  if (x$distribution == "BiGam"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-symmetric discrepancy (SD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the ", x$distribution, " distribution\n")
    cat(paste0('\np.value : ', x$p.GOF, '\n'))
  }

  if (x$distribution == "Dirichlet"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-symmetric discrepancy (SD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the ", x$distribution, " distribution\n")
    cat(paste0('\np.value : ', x$p.GOF, '\n'))
  }
}




#' Summarizing MLEce function
#'
#' @description \code{summary} method for a class "MLEce".
#'
#' @param object an object of class "MLEce" made by the function \code{MLEce}.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "MLEce".
#' @param digits a numeric number of significant digits.
#' @method summary MLEce
#'
#' @return \code{summary} describes information about MLEce. (quantile statistics, correlation, estimates)
#' @examples
#' datt = rBiGam(100, c(1,4,5))
#' res = MLEce(datt, "BiGam", boots = 50)
#' summary(res)
#'
#' @export
summary.MLEce <- function(object,...){
    distribution <- object$distribution
    if(distribution == "BiWei"){
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)[1,2]
        ce_est <- object$estimation
        CI <- object$CI
        alpha <- object$CI.alpha
        se <- object$se
        un_lab <- paste0(alpha/2, '%'); up_lab <- paste0(1-alpha/2, '%')

        est_mat <- matrix(c(ce_est, se, CI[1,] , CI[2,] ), nrow = 5)
        colnames(est_mat) <- c("MLEce", "Std. Error", un_lab, up_lab ); rownames(est_mat) <- names(ce_est)

        result <- list("stat_summary"=stat_summary, "correlation" = rho,
                       "distribution"=distribution, "p.GOF" = object$p.GOF, "estimation_matrix" = est_mat,
                       "mle" = object$mle)
    }

    if(distribution == "BiGam"){
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)[1,2]
        ce_est <- object$estimation
        CI <- object$CI
        alpha <- object$CI.alpha
        se <- object$se
        un_lab <- paste0(alpha/2, '%'); up_lab <- paste0(1-alpha/2, '%')

        est_mat <- matrix(c(ce_est, se, CI[1,] , CI[2,] ), nrow = 3)
        colnames(est_mat) <- c("MLEce","Std. Error", un_lab, up_lab ); rownames(est_mat) <- names(ce_est)

        result <- list("stat_summary"= stat_summary, "correlation" = rho,
                       "distribution"=distribution, "estimation_matrix"=est_mat)

    }

    if(distribution == "Dirichlet"){
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)
        ce_est <- object$estimation

        CI <- object$CI
        alpha <- object$CI.alpha
        un_lab <- paste0(alpha/2, '%'); up_lab <- paste0(1-alpha/2, '%')

        est_mat <- matrix(c(ce_est, CI[1,] , CI[2,] ), nrow = length(ce_est))
        colnames(est_mat) <- c("MLEce", un_lab, up_lab )
        if (is.null(names(ce_est))) {
            rownames(est_mat) <- paste0("par",seq(1, length(ce_est)))
        } else rownames(est_mat) <- names(ce_est)



        result <- list("stat_summary"= stat_summary, "correlation" = rho,
                       "distribution"=distribution, "estimation_matrix"=est_mat)

    }

    class(result) <- "summary.MLEce"
    result
}
#' @rdname summary.MLEce
#' @method print summary.MLEce
#' @export
print.summary.MLEce <- function(x, digits = max(3, getOption("digits") - 3), ...){
    distribution <- x$distribution
    if(x$distribution == "BiWei"){
        cat("\nDistribution: Bivariate Weibull\n")
        cat("\n")
        cat("Summary:\n")
        cat("Quantile\n")
        print(round(x$stat_summary, digits = digits))
        cat("\nCorrelation\n")
        print(round(x$correlation, digits = digits))
        cat("\nEstimation\n")
        print(round(x$estimation_matrix, digits = digits))
        #if(x$mle == TRUE){
        #    cat("\n")
        #    print(round(x$estimation_matrix.ml, digits = digits))
        #}
        #if(x$mle == TRUE){
        #    cat("\nGCVM test of MLE\n")
        #    print(x$p.GOF.ml)
        #}
    }
    if(x$distribution == "BiGam"){
        cat("\nDistribution: Bivariate Gamma\n")
        cat("\n")
        cat("Summary:\n")
        cat("Quantile\n")
        print(round(x$stat_summary, digits = digits))
        cat("\nCorrelation\n")
        print(round(x$correlation, digits = digits))
        cat("\nEstimation\n")
        print(round(x$estimation_matrix, digits = digits))
    }

    if(x$distribution == "Dirichlet"){
        cat("\nDistribution:", distribution, "\n")
        cat("\n")
        cat("Summary:\n")
        cat("Quantile\n")
        print(round(x$stat_summary, digits = digits))
        cat("\nCorrelation\n")
        print(round(x$correlation, digits = digits))
        cat("\nEstimation\n")
        print(round(x$estimation_matrix, digits = digits))
    }
}

#' Providing some plots for MLEce
#' @description \code{plot} method for a class "MLEce".
#'
#' @param x an object of class "MLEce" made by the function \code{MLEce}.
#' @param which if a subset of the plots is required, specify a subset of 1:4
#' @param ask logical; if TRUE, the user is asked before each plot.
#' @param ... not used, but exists because of the compatibility.
#'
#' @method plot MLEce
#' @details
#' The first figure is a boxplot for given data.
#' The second figure is a contour line drawn by the probability density function of the estimated parameter based on MLEce.
#' the x-axis is the first column of data and the y-axis is the second column of data.
#' The third figure is a marginally fitted probability density plot for the first column of input data.
#' It provides a fitted line for each of CME, MLE and MLEce.
#' The fourth figure is a marginally fitted probability density plot for the second column of input data.
#' It can also provide a fitted line for each of CME, MLE and MLEce.
#' @return returns plots for MLEce which describe "details".
#' @import ggplot2 reshape
#' @importFrom LaplacesDemon ddirichlet
#' @importFrom graphics axis hist lines
#' @examples
#' datt = rBiGam(100, c(1,4,5))
#' res = MLEce(datt, "BiGam", boots = 50)
#' plot(res, c(1))
#' @export
plot.MLEce <- function(x, which=c(1,2,3,4),
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
    show <- rep(FALSE, 4)

    dat <- as.data.frame(x$data)
    if(dim(dat)[2] > 2){
      which = c(1)
      }
    show[which] <- TRUE
    dots <- list(...)
    nmdots <- names(dots)

    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if(x$distribution == "BiWei"){
        cme <- x$cme
        ce_est <- x$estimation

        if (is.null(colnames(dat))){xy.name = c("x", "y")}
        else xy.name = colnames(dat)

        if(show[1]==TRUE){
          meltData <- melt(dat)
          colnames(meltData) = c('variable', 'value')
          p <- ggplot(meltData, aes_string('variable', 'value')) + geom_boxplot() + facet_wrap('variable', scales="free") + labs(title = 'Boxplots for Data (Bivariate Weibull)',) + xlab('columns')
          print(p)
        }

        if(show[2]==TRUE){
            gr.x = seq(from=min(dat[,1]) * 0.9,to=max(dat[,1])*1.1,length.out = 250)
            gr.y = seq(from=min(dat[,2]) * 0.9,to=max(dat[,2])*1.1,length.out = 250)
            pars = expand.grid(gr.x,gr.y)

            bidensity_xy = t(apply(matrix(1:nrow(pars),ncol=1),1, function(k){
                exp(dBiWei(x$estimation,matrix(c(pars[k,1],pars[k,2]),ncol=2)))} ) )

            gr.z = matrix(bidensity_xy,ncol=length(gr.y),byrow=T)
            contour(gr.x,gr.y,t(gr.z),levels = c(0.001,seq(0.005,0.05,by=0.005)), main = 'Contour plot',
                    xlab = xy.name[1], ylab = xy.name[2])
            points(x$data[,1:2],pch=20,col="black",cex=0.5)
        }
        if(show[3]==TRUE){
            main.title = paste0("Marginal plot of ", xy.name[1])
            gr.x = seq(from=min(dat[,1]) * 0.9,to=max(dat[,1])*1.1,length.out = 250)
            plot(gr.x, dweibull(gr.x,shape=cme[1],scale=cme[2]),type="l",lty=3,col=2,xlab=xy.name[1],ylab="", main = main.title)
            points(gr.x,dweibull(gr.x,shape=ce_est[1],scale=ce_est[2]),type="l",lty=2,xlab="", ylab = "", col=3)
            rug(x$data[,1])
            legend("topleft", c("CME",expression(MLE[CE])), col=c(2,3), lwd=2,lty=c(3,2), cex=1.1)

        }

        if(show[4]==TRUE){
            main.title = paste0("Marginal plot of ", xy.name[2])
            gr.x = seq(from=min(dat[,2]) * 0.9,to=max(dat[,2])*1.1,length.out = 250)
            plot(gr.x, dweibull(gr.x,shape=cme[3],scale=cme[4]),type="l",lty=3,col=2,xlab=xy.name[2],ylab="", main = main.title)
            points(gr.x,dweibull(gr.x,shape=ce_est[3],scale=ce_est[4]),type="l",lty=2,xlab="", ylab = "", col=3)
            rug(x$data[,2])
            legend("topleft", c("CME",expression(MLE[CE])), col=c(2,3), lwd=2,lty=c(3,2), cex=1.1)
        }
    }
    if(x$distribution == "BiGam"){
        ce_est <- x$estimation

        mme = BiGam_MME(dat)
        est.mle = BiGam_MLE_log(mme,dat,tol2=mme)

        if (is.null(colnames(dat))){xy.name = c("x", "y")}
        else xy.name = colnames(dat)

        if(show[1]==TRUE){
          meltData <- melt(dat)
          colnames(meltData) = c('variable', 'value')
          p <- ggplot(meltData, aes_string('variable', 'value')) + geom_boxplot() + facet_wrap('variable', scales="free") + labs(title = 'Boxplots for Data (Bivariate Gamma)',) + xlab('columns')
          print(p)
        }

        if(show[2]==TRUE){
          ##Contour plots.
          gr.x=seq(from=min(dat[,1]*0.9),to=max(dat[,1]*1.1),length.out = 250)
          gr.y=seq(from=min(dat[,2]*0.9),to=max(dat[,2]*1.1),length.out = 250)
          z = outer(gr.x,gr.y,Vectorize(function(gr.x,gr.y){dBiGam(ce_est,gr.x,gr.y,log=F)}))
          contour(gr.x,gr.y,z)
          title(main="Estimated contour plot (MLEce)")
          title(xlab=xy.name[1],ylab=xy.name[2],cex.lab=1.2)
          points(dat[,1],dat[,2], pch=20)
        }

        if(show[3]==TRUE){
          hist(dat[,1], # histogram
               col = 'skyblue',
               border="black",
               prob = TRUE, # show densities instead of frequencies
               main = "Histogram & Density curve (1st column)")

          lines(density(dat[,1]), # density plot
                lwd = 2, # thickness of line
                col = "blue")
        }
        if(show[4]==TRUE){
          hist(dat[,2], # histogram
               col = 'skyblue',
               border="black",
               prob = TRUE, # show densities instead of frequencies
               main = "Histogram & Density curve (2nd column)")

          lines(density(dat[,2]), # density plot
                lwd = 2, # thickness of line
                col = "blue")
        }

    }

    if(x$distribution == "Dirichlet"){

        ce_est <- x$estimation

        if (is.null(colnames(dat))){xy.name = c("x", "y")}
        else xy.name = colnames(dat)

        if(show[1]==TRUE){
          meltData <- melt(dat)
          colnames(meltData) = c('variable', 'value')
          p <- ggplot(meltData, aes_string('variable', 'value')) + geom_boxplot() + facet_wrap('variable', scales="free") + labs(title = 'Boxplots for Data (Dirichlet)',) + xlab('columns')
          print(p)
        }

        if(show[2]==TRUE){
          gr.x=seq(from= 0.05,to=0.95,length.out = 250)
          gr.y= 1-gr.x
          z = outer(gr.x,gr.y,Vectorize(function(gr.x, gr.y){ddirichlet(cbind(gr.x, gr.y), ce_est, log=F)}))
          contour(z)
          title(main="Estimated contour plot (MLEce)")
          title(xlab=xy.name[1],ylab=xy.name[2],cex.lab=1.2)
        }

        if(show[3]==TRUE){
          hist(dat[,1], # histogram
               col = 'skyblue',
               border="black",
               prob = TRUE, # show densities instead of frequencies
               main = "Histogram & Density curve (1st column)")

          lines(density(dat[,1]), # density plot
                lwd = 2, # thickness of line
                col = "blue")
        }
        if(show[4]==TRUE){
          hist(dat[,2], # histogram
               col = 'skyblue',
               border="black",
               prob = TRUE, # show densities instead of frequencies
               main = "Histogram & Density curve (2nd column)")

          lines(density(dat[,2]), # density plot
                lwd = 2, # thickness of line
                col = "blue")
        }



    }
    invisible()
}

#' compute MLEce and MLE
#'
#' @description \code{computeTime} performs a benchmark of MLEce and MLE on a given dataset.
#'
#' @param data a data set.
#' @param distribution a character string of a distribution assuming that data set comes from.
#' @param coef_out if TRUE, estimated parameters are printed. Default is False.
#'
#' @examples
#' dat <- rBiWei(n=30, c(4,3,3,4,0.6))
#' computeTime(dat, "BiWei")
#'
#' @return a numeric matrix. This matrix include estimated parameters and time.
#'
#' @export
computeTime <- function(data, distribution, coef_out = FALSE){
    distls <- c("BiWei", "BiGam", "Dirichlet")

    if (!is.element(distribution, distls)) {
        stop("unsupported distribution")
    }
    else {
        initb <- FALSE
    }


    if(distribution == "BiWei"){
        est_mat = matrix(NA, nrow = 5, ncol = 2)
        time_vec = c(NA, NA)

        # time for MLEce
        ptm = Sys.time()
        ce_res = BiWei_CE(data)
        time_vec[1] = as.numeric(Sys.time() - ptm)
        ce_est = ce_res$estimation
        est_mat[,1] = ce_est


        # time for MLE
        cme = ce_res$cme
        ptm = Sys.time()
        mle_est = BiWei_ML(data, cme)
        time_vec[2] = as.numeric(Sys.time() - ptm)
        est_mat[,2] = mle_est


        rownames(est_mat) = names(ce_est); colnames(est_mat) = c("MLEce", "MLE")
        results = rbind(est_mat, time_vec)
    }

    if(distribution == "BiGam"){
        est_mat = matrix(NA, nrow = 3, ncol = 2)
        time_vec = c(NA, NA)
        gam_mme = BiGam_MME(data)

        # time for MLEce
        ptm = Sys.time()
        ce_est = BiGam_CE(gam_mme, data)$estimation
        time_vec[1] = as.numeric(Sys.time() - ptm)
        est_mat[,1] = ce_est


        # time for MLE
        ptm = Sys.time()
        mle_est = BiGam_MLE(gam_mme, data)
        time_vec[2] = as.numeric(Sys.time() - ptm)
        est_mat[,2] = mle_est


        rownames(est_mat) = c("shape1", "shape2", "scale"); colnames(est_mat) = c("MLEce", "MLE")
        results = rbind(est_mat, time_vec)
    }

    if(distribution == "Dirichlet"){
        dn = dim(data)[2]
        est_mat = matrix(NA, nrow = dn, ncol = 2)
        time_vec = c(NA, NA)
        est_mme = Diri_MME(data)

        # time for MLEce
        ptm = Sys.time()
        ce_est = Diri_CE(data)$estimation
        time_vec[1] = as.numeric(Sys.time() - ptm)
        est_mat[,1] = ce_est


        # time for MLE
        ptm = Sys.time()
        mle_est = Diri_MLE(data)
        time_vec[2] = as.numeric(Sys.time() - ptm)
        est_mat[,2] = mle_est


        rownames(est_mat) <- paste0("par",seq(1, length(ce_est))); colnames(est_mat) = c("MLEce", "MLE")
        results = rbind(est_mat, time_vec)
    }

    if (coef_out == TRUE){
        return(results)
    } else {
        time_mat = matrix(time_vec)
        rownames(time_mat) = c("MLEce", "MLE")
        return(time_mat)
    }

}


#' The flood events data of the Madawaska basin.
#'
#' It represents the flood events data of the Madawaska basin from 1911 to 1995. (Yue, 2001).
#'
#' @docType data
#' @usage data(flood, package = "MLEce")
#' @format A 2 variables data frame with 77 observations.
#'
#' @references Yue, S. (2001). <doi:10.1002/hyp.259>.
"flood"



#' The counts data of the frequency of occurrence of different kinds of pollen grains.
#'
#' It represents the counts data of the frequency of occurrence of different kinds of pollen grains. (Mosimann, 1962).
#'
#' @docType data
#' @usage data(fossil_pollen, package = "MLEce")
#' @format A 4 variables data frame with 73 observations.
#'
#' @references Mosimann, James E. (1962) <doi:10.1093/biomet/49.1-2.65>.
"fossil_pollen"
