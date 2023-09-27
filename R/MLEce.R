#' Calculating MLEces for three different distributions
#'
#' The closed-form estimators (MLEces) are calculated for three distributions: bivariate gamma, bivariate Weibull and multivariate Dirichlet. 
#' 
#' @details
#' Based on root n-consistent estimators, the closed-form estimators (MLEces) are calculated for the parameters in bivariate gamma, bivariate Weibull and multivariate Dirichlet distributions whose maximum likelihood estimators (MLEs) are not in closed forms. The MLEces are strong consistent and asymptotic normally like the corresponding MLEs, but their calculation are much faster than MLEs. For the bivariate gamma and multivariate Dirichlet distribution, their root n-consistent estimators are the corresponding method of moments estimators (MMEs). The correlation-based estimators (CMEs) are applied as root n-consistent estimators in the bivariate Weibull distribution.
#' @references Kim, H.-M., Jang, Y.-H., Arnold, B. C. and Zhao, J. (2023) New efficient estimators for the Weibull distribution. \emph{Communications in Statistics - Theory and Methods}, 1-26. 
#' @references Jang, Y.-H., Zhao, J., Kim, H.-M., Yu, K., Kwon, S.and Kim, S. (2023) New closed-form efficient estimator for the multivariate gamma distribution. \emph{Statistica Neerlandica}, 1–18.
#' @references Chang, J. H., Lee, S. K. and Kim, H.-M. (2023) An asymptotically efficient closed–form estimator for the multivariate Dirichlet distribution. submitted.
#' @param data a numeric matrix.
#' @param distname a character indicating which distribution to be fitted. \code{BiGam} stands for the bivariate gamma distribution, \code{BiWei} stands for the bivariate Weibull distribution, and \code{Dirichlet} stands for the multivariate Dirichlet distribution.
#' @return \code{MLEce} returns an object of class \code{"MLEce"}. The object class \code{"MLEce"} is a list containing the following components.
#' \item{distribution }{a character string of a distribution assuming that data set comes from and the data was fitted to.}
#' \item{estimation}{the estimated values of parameters in assigned distribution.}
#' @import stats
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom graphics title contour par points
#' @importFrom grDevices dev.interactive devAskNewPage
#' @importFrom graphics legend rug
#' @examples
#' #bivariate gamma distribution
#' data_BiGam <-  rBiGam(100, c(1,4,5))
#' res_BiGam  <- MLEce(data_BiGam, "BiGam")
#' print(res_BiGam)
#' @examples
#' data(flood)
#' est_BiGam <- MLEce(flood, "BiGam")
#' print(est_BiGam)
#' @examples
#' #bivariate Weibull distribution
#' data_BiWei <- rBiWei(n=30, c(4,3,3,4,0.6))
#' res_BiWei <- MLEce(data_BiWei, "BiWei")
#' print(res_BiWei)
#' @examples 
#' #real data example
#' data(airquality)
#' air_data <- airquality[ ,3:4]
#' air_data[ ,2] <- air_data[ ,2]*0.1
#' est_BiWei <- MLEce(air_data, "BiWei")
#' print(est_BiWei)
#' @examples
#' #Dirichlet distribution
#' data_Diri <- LaplacesDemon::rdirichlet(n=60, c(1,2,3))
#' res_Diri <- MLEce(data_Diri, "Dirichlet")
#' print(res_Diri)
#' @examples
#' data(fossil_pollen) 
#' #real data example
#' fossil_data <- fossil_pollen/rowSums(fossil_pollen)
#' eps <- 1e-10
#' fossil_data <- (fossil_data +eps)/(1+2*eps)
#' est_fossil <- MLEce(fossil_data, "Dirichlet")
#' print(est_fossil)
#' @export
MLEce <-  function(data, distname){
   distlist <- c("BiGam", "BiWei",  "Dirichlet")
 if (!is.element(distname, distlist))
   {
    stop("unsupported distribution")
   }


  #check argument data
  if (!( (is.matrix(data) & is.numeric(data)) | is.data.frame(data) ) )
  { stop("data must be a numeric matrix or dataframe") }
   
  n = dim(data)[1]
  
  if (distname == "BiGam") { 
    est <- BiGam_CE( data   )
    est_values=est$estimation
    hessian_mat=est$Hessian
    asym.var = -qr.solve(hessian_mat)
    par.se = sqrt(diag(asym.var)/n)
   }
  else if(distname == "BiWei") { 
    est=BiWei_CE( data )
    est_values=BiWei_CE( data )$estimation
    hessian_mat=BiWei_info(est$cme, data, type = "hessian")$hessian
    #asym.var = -qr.solve(hessian_mat)
    par.se = c("calculated by boostrap method")
   }

  else if(distname == "Dirichlet") {
   est = Diri_CE_bt( data )
   est_values=est$estimation
   Invhess = est$InvHes
   hessian_mat = Invhess
   #asym.var=-Invhess #beta tilde value
   par.se = c("standard error of Dirichlet is not available")
  }
  stat_summary <- t(apply(data, 2, quantile))
  est_results <- list(distribution = distname, estimation=est_values, 
                      hessian_mat=hessian_mat, samplesize = n,
                      se=par.se, stat_summary=stat_summary, data=data)

  class(est_results)<- "MLEce"
  est_results
}

#' @method print MLEce
#' @export
print.MLEce <- function(x, digits = max(4, getOption("digits") - 4), ...){

    if (x$distribution == "BiGam"){
        cat("\nDistribution: Bivariate Gamma\n")
        est = round(x$estimation, digits = digits)
        #colnames(est)=c("shape1", "shape2", "scale")
        cat("\nEstimate:\n")
        #cat("shape1 = ", est[1], ", ", "shape2 = ", est[2], ", ", "scale = ", est[3])
        print(est)
    }

    if(x$distribution == "BiWei"){
        cat("\nDistribution: Bivariate Weibull\n")
        est = round(x$estimation, digits = digits)
        #colnames(est)=c("alpha1", "beta1", "alpha2","beta2", "delta")
        cat("\nEstimate:\n")
        print(est)
    }

    if (x$distribution == "Dirichlet"){
        cat("\nDistribution: Dirichlet\n")
        est = round(x$estimation, digits = digits)
        #termno <- length(x$estimation)
        #colnames(est) <- c(paste("alpha", 1:termno, sep=""))
        cat("\nEstimate:\n")
        print(est)
    }
}
#' Performing closed-form estimators against other methods
#' @param data a numeric matrix.
#' @param distname a character indicating which distribution to be fitted. \code{BiGam} stands for the bivariate gamma distribution, \code{BiWei} stands for the bivariate Weibull distribution, and \code{Dirichlet} stands for the Dirichlet distribution.
#' @param methods a vector of methods: two characters among \code{"MLEce"} (efficient closed-form estimator), \code{"MLE"}, \code{"MME"} and \code{"CME"}. \code{MLEce} stands for efficient closed-form estimators, , \code{"MLE"} (maximum likelihood estimator), \code{"MME"} (method of moments estimator) and \code{"CME"} (correlation based method  estimator for the bivariate Weibull distribution). 
#' @return A matrix with estimate and time in seconds per method for assigned distributions.
#' @import mvtnorm
#' @import sirt
#' @examples
#' #bivariate gamma distribution
#' data_BiGam= rBiGam(100, c(1,4,5))
#' benchMLEce(data_BiGam, distname="BiGam", methods=c("MLEce","MME"))
#' @examples
#' #bivariate Weibull distribution
#' data_BiWei <- rBiWei(n=50, c(4,3,3,4,0.6))
#' benchMLEce(data_BiWei, distname="BiWei", methods=c("MLE","CME"))
#' @examples
#' #multivariate Dirichlet distribution
#' data_Diri <- LaplacesDemon::rdirichlet(80, c(3,4,1,3,4))
#' benchMLEce(data_Diri, distname="Dirichlet", methods=c("MLEce","MLE"))
#' @export
benchMLEce <-  function(data, distname, methods){
   distlist <- c( "BiGam", "BiWei", "Dirichlet")
  if (!is.element(distname, distlist))
    {
     stop("unsupported distribution")
    }

   #check argument methods
   methods <- sapply(methods, function(x) match.arg(x, c("MLEce","MLE","MME","CME")))

   restime <- numeric(length(methods))
   resestimate <- NULL
 
   #check argument data
   if (!( is.matrix(data) & is.numeric(data) ) )
     { stop("data must be a numeric matrix") }

   for(i in 1:length(methods))
   {
     if (distname == "BiGam") { 
       if(methods[i] == "MLEce")
      {
        timebeg <- Sys.time()
        est_values=BiGam_CE( data   )$estimation
        timeend <- Sys.time()
        restime[i] <- timeend - timebeg
        resestimate <- c(resestimate, list(est_values))
      }

      if(methods[i] == "MLE")
     {
       timebeg <- Sys.time()
       est_values=BiGam_MLE( data   )
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
    
      if(methods[i] == "MME")
     {
       timebeg <- Sys.time()
       est= BiGam_MME( data   )
       est_values=c(est$alpha1,est$alpha2,est$beta)
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
   }
  if (distname == "BiWei") { 
      if(methods[i] == "MLEce")
     {
       timebeg <- Sys.time()#proc.time()
       est_values=BiWei_CE( data   )$estimation
       timeend <- Sys.time()#proc.time()
       
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }

      if(methods[i] == "MLE")
     {
       timebeg <- Sys.time()
       est_values=BiWei_MLE( data   )
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
    
      if(methods[i] == "CME")
     {
       timebeg <- Sys.time()
       est_values=BiWei_CME( data   )
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
   }

  if (distname == "Dirichlet") { 
      if(methods[i] == "MLEce")
     {
       timebeg <- Sys.time()
       est_values=Diri_CE( data   )$estimation
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }

      if(methods[i] == "MLE")
     {
       timebeg <- Sys.time()
       est_values <- as.vector(Diri_MLE( data   ))
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
    
      if(methods[i] == "MME")
     {
       timebeg <- Sys.time()
       est_values <- as.vector(Diri_MME( data   ) )
       timeend <- Sys.time()
       restime[i] <- timeend - timebeg
       resestimate <- c(resestimate, list(est_values))
     }
   }
  }

  resestimate <- simplify2array(resestimate)
  if(is.vector(resestimate))
  {
    resestimate <- rbind(resestimate)
    rownames(resestimate) <- colnames(resestimate)[1]
  }
  colnames(resestimate) <- methods
  results <- list(distribution=distname, estresults=rbind(resestimate, "time"=restime))
  class(results) <- "benchMLEce"
  return(results)
}

#' @method print benchMLEce
#' @export
print.benchMLEce <- function(x, digits = max(3, getOption("digits") - 3), ...){

    if (x$distribution == "BiGam"){
        cat("\nDistribution: Bivariate Gamma\n")
        rownames(x$estresults) <- c("shape1","shape2","scale","time")
        print(x$estresults)
    }

    if(x$distribution == "BiWei"){
        cat("\nDistribution: Bivariate Weibull\n")
        rownames(x$estresults) <- c("alpha1", "beta1", "alpha2", "beta2","delta","time")
        print(x$estresults)
    }

    if (x$distribution == "Dirichlet"){
        termno <- dim(x$estresults)[1]
        cat("\nDistribution: Dirichlet\n")
        rownames(x$estresults) <- c(paste("alpha", 1:(termno-1), sep=""), "time")
        print(x$estresults)
    }
}



#' Getting estimated values of efficient closed-form estimators
#' @param object an object of class \code{"MLEce"} made by the function \code{MLEce}.
#' @param digits a numeric number of significant digits.
#' @param ... not used, but exists because of the compatibility.
#' @return a numeric vector or a list, containing assigned distribution and estimated values, is given.
#' @examples 
#' data_BiGam = rBiGam(100, c(1,4,5))
#' res_BiGam = MLEce(data_BiGam, "BiGam")
#' coef(res_BiGam)
#' @examples
#' data_BiWei = rBiWei(n=50, c(4,3,3,4,0.6))
#' est_BiWei <-MLEce(data_BiWei, "BiWei")
#' coef(est_BiWei)
#' @examples
#' data_Diri <- LaplacesDemon::rdirichlet(n=60, c(3,1,2,4))
#' est_Diri <- MLEce(data_Diri, "Dirichlet")
#' coef(est_Diri)
#' @export
coef.MLEce <- function(object, digits = max(3, getOption("digits") - 3), ...){
  if (object$distribution == "BiGam"){
   distrname=c("Bivariate Gamma")
   est = round(object$estimation, digits = digits)
  }
  if (object$distribution == "BiWei"){
    distrname=c("Bivariate Weibull")
    est = round(object$estimation, digits = digits)
  }
  if (object$distribution == "Dirichlet"){
    distrname=c("Dirichlet")
    est = round(object$estimation, digits = digits)
  }
  result = list(distribution=distrname,coef.out = est)
  return(result)
}

#' Getting confidence intervals for efficient closed-form estimators
#' @details
#' The confidence interval is obtained by bootstrap method for the estimated parameters in the assigned distribution.
#' @param object an object of class \code{"MLEce"} made by the function \code{MLEce}.
#' @param bootsize a numeric value for the steps in the bootstrap method; default value is 1,000.
#' @param level a numeric value between 0 and 1 for controlling the significance level of confidence interval; default value is 0.95.
#' @return a numeric a list is given, containing assigned distribution, confidence intervals and alpha which is equal to one minus the significance level.
#' @examples
#' data(flood)
#' est_BiGam <- MLEce(flood, "BiGam")
#' confCI(est_BiGam)
#' @examples
#' datt = rBiWei(n=50, c(4,3,3,4,0.6))
#' est_BiWei <-MLEce(datt, "BiWei")
#' confCI(est_BiWei)
#' @examples
#' data_Diri <- LaplacesDemon::rdirichlet(n=60, c(3,1,2,4))
#' est_Diri <- MLEce(data_Diri, "Dirichlet")
#' confCI(est_Diri)
#' @export
confCI <- function(object, bootsize = 1000, level = 0.95) {
    ce_est = object$estimation
    #hess = object$hessian_mat
    n = object$samplesize
    CI.alpha = 1-level
    boots = bootsize
  if (object$distribution == "BiGam"){


        distrname=c("Bivariate Gamma")
        ce_pb = matrix(0,nrow=boots,ncol=3)
        for (i in 1:boots){
            rg_dat = rBiGam(n, ce_est)
            ce_pb[i,] = BiGam_CE(rg_dat)$estimation
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)
    }

  if (object$distribution == "BiWei"){

        distrname=c("Bivariate Weibull")
        ce_pb = matrix(0,nrow=boots,ncol=5)
        for (i in 1:boots){
            rg_dat = rBiWei(n, ce_est)
            ce_pb[i,] = BiWei_CE(rg_dat)$estimation
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)
     }

  if (object$distribution == "Dirichlet"){
     
     ce_pb = matrix(0,nrow=boots, ncol = length(ce_est))
     distrname=c("Dirichlet")
        for (i in 1:boots){
          
          stbz <- function(x, eps = 1e-10) {
            (x + eps) / (1 + 2 * eps)
          }
            rg_dat = stbz(rdirichlet(n, ce_est))
            ce_pb[i,] = Diri_CE_bt(rg_dat)$estimation
        }
        CI = round(apply(ce_pb,2,quantile,c(CI.alpha/2, 1- CI.alpha/2)),3)
    }
  CIresults = list(distribution=distrname, CI=CI, alpha=CI.alpha)
  return(CIresults)
}

#' Goodness-of-fit test for the efficient closed-form estimators
#' @details Generalized Cramer-von Mises test (chiu and Liu, 2009) is applied to do the goodness-of-fit test for multivariate distributions. For the bivariate gamma and Dirichlet distributions, the L2-symmetric discrepancy (SD2) statistics are applied. But the L2-centred discrepancy (CD2) statistics are applied in the bivariate Weibull distribution.
#' @param x an object of class "MLEce" made by the function \code{MLEce}.
#' @param digits a numeric number of significant digits.
#' @param ...  additional arguments affecting the goodness-of-fit test.
#' @references chiu, S. N. and Liu, K. I. (2009) Generalized Cramer-Von Mises goodness-of-fit tests for multivariate distributions. \emph{Computational Statistics and Data Analysis}, 53, 3817-3834.
#' @examples
#' data_BiGam <-  rBiGam(100, c(1,4,5))
#' res_BiGam  <- MLEce(data_BiGam, "BiGam")
#' gof(res_BiGam)
#' @examples
#' datt = rBiWei(n=50, c(4,3,3,4,0.6))
#' est_BiWei <-MLEce(datt, "BiWei")
#' gof(est_BiWei)
#' @examples
#' data_Diri <- LaplacesDemon::rdirichlet(n=60, c(3,1,2,4))
#' est_Diri <- MLEce(data_Diri, "Dirichlet")
#' gof(est_Diri)
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
      
      stbz <- function(x, eps = 1e-10) {
        (x + eps) / (1 + 2 * eps)
      }
      
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

  if (x$distribution == "BiGam"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-symmetric discrepancy (SD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the bivariate gamma distribution\n")
    cat(paste0('\np.value : ', x$p.GOF, '\n'))
  }
  if (x$distribution == "BiWei"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-centred discrepancy (CD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the bivariate Weibull distribution\n")
    cat(paste0('\np.value : ',x$p.GOF, '\n'))
  }
  if (x$distribution == "Dirichlet"){
    cat("\nGoodness of fit test\n")
    cat("generalized Cramer-Von Mises test based on L2-symmetric discrepancy (SD2 statistics)\n")
    cat("\nalternative hypothesis : data not follows the Dirichlet distribution\n")
    cat(paste0('\np.value : ', x$p.GOF, '\n'))
  }
}



#' Summarizing effective closed-form estimation function
#'
#' @description \code{summary} method for a class "MLEce".
#' @param object an object of class "MLEce" made by the function \code{MLEce}.
#' @param ... not used, but exists because of the compatibility.
#' @param x an object of class "summary.MLEce".
#' @param digits a numeric number of significant digits.
#' @method summary MLEce
#'
#' @return \code{summary} presents information about effective closed-form estimators calculated by \code{MLEce} containing the following components.
#' \item{Distribution}{the distribution assigned to fit the data to.}
#' \item{Quantile}{a numeric vector describing the data set with min, 1st quantile, median, 3rd quantile, and max values.}
#' \item{Correlation}{correlation coefficient between two vectors of the data}
#' \item{Estimation}{estimated values of parameters, standard error and confidence intervals are given.}

#' @examples
#' #bivariate gamma distribution
#' data(flood)
#' est_res1 <- MLEce(flood, "BiGam")
#' summary(est_res1)
#' @examples
#' #bivariate Weibull distribution
#' datt = rBiWei(n=50, c(2,3,3,4,0.4))
#' est_res2 <-MLEce(datt, "BiWei")
#' summary(est_res2)
#' @examples
#' #Dirichilet distribution
#' data(fossil_pollen) 
#' fossil_data <- fossil_pollen/rowSums(fossil_pollen)
#' eps <- 1e-10
#' fossil_data <- (fossil_data +eps)/(1+2*eps)
#' est_res3 <- MLEce(fossil_data, "Dirichlet")
#' summary(est_res3)
#' @export
summary.MLEce <- function(object, ...){
    distribution <- object$distribution

    if(distribution == "BiGam"){
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)[1,2]
        ce_est <- object$estimation
        CI_results <- confCI(object)
        CI <- CI_results$CI
        alpha <- CI_results$alpha
        se <- object$se
        un_lab <- paste0(100*(alpha/2), '%'); up_lab <- paste0(100*(1-alpha/2), '%')

        est_mat <- matrix(c(ce_est, se, CI[1,] , CI[2,] ), nrow = 3)
        colnames(est_mat) <- c("MLEce","Std. Error", un_lab, up_lab ); rownames(est_mat) <- names(ce_est)

        result <- list("stat_summary"= stat_summary, "correlation" = rho,
                       "distribution"=distribution, "estimation_matrix"=est_mat)

    }

    if(distribution == "BiWei"){
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)[1,2]
        ce_est <- object$estimation
        CI_results <- confCI(object)
        CI <- CI_results$CI
        alpha <- CI_results$alpha
        
        n = object$samplesize
        ce_pb = matrix(0,nrow=2000,ncol=5)
        for (i in 1:2000){
          rg_dat = rBiWei(n, ce_est)
          ce_pb[i,] = BiWei_CE(rg_dat)$estimation
        }
        se <- c(sd(ce_pb[,1]),sd(ce_pb[,2]),sd(ce_pb[,3]),sd(ce_pb[,4]),sd(ce_pb[,5]))
        
        un_lab <- paste0(100*(alpha/2), '%'); up_lab <- paste0(100*(1-alpha/2), '%')

        est_mat <- matrix(c(ce_est, se, CI[1,] , CI[2,] ), nrow = 5)
        colnames(est_mat) <- c("MLEce", "Std. Error", un_lab, up_lab ); rownames(est_mat) <- names(ce_est)

        result <- list("stat_summary"=stat_summary, "correlation" = rho,
                       "distribution"=distribution, "p.GOF" = object$p.GOF, "estimation_matrix" = est_mat,
                       "mle" = object$mle)
    }


    if(distribution == "Dirichlet"){
        
        stat_summary <- object$stat_summary
        colnames(stat_summary) <- c("Min","1st Qu","Median","3rd Qu","Max")
        rho <- cor(object$data)
        ce_est <- object$estimation
        n = object$samplesize
        ce_pb = matrix(0,nrow=2000, ncol = length(ce_est))
        
        for (i in 1:2000){
          
          stbz <- function(x, eps = 1e-10) {
            (x + eps) / (1 + 2 * eps)
          }
          rg_dat = stbz(rdirichlet(n, ce_est))
          ce_pb[i,] = Diri_CE_bt(rg_dat)$estimation
        }
        
        se <- c()
        for(j in 1:length(ce_est)){se[j] <- sd(ce_pb[,j])}
        
        CI_results <- confCI(object)
        CI <- CI_results$CI
        alpha <- CI_results$alpha
        un_lab <- paste0(100*(alpha/2), '%'); up_lab <- paste0(100*(1-alpha/2), '%')

        est_mat <- matrix(c(ce_est, se, CI[1,] , CI[2,] ), nrow = length(ce_est))
        colnames(est_mat) <- c("MLEce", "Std. Error",un_lab, up_lab )
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


#' Providing some plots for effective closed-form estimators
#' @description \code{plot} method for a class "MLEce".
#'
#' @param x an object of class "MLEce" made by the function \code{MLEce}.
#' @param which if a subset of the plots is required, specify a subset of 1:4.
#' @param ask logical; if TRUE, the user is asked before each plot.
#' @param ... not used, but exists because of the compatibility.
#'
#' @method plot MLEce
#' @details
#' The boxplot for given data is presented first with \code{which=1}. For \code{which=2}, a contour line is drawn by the probability density function of the estimated parameter based on effective closed-form estimators. In the counter plot, the x-axis is the first column of data and the y-axis is the second column of data. For \code{which=3}, a marginally fitted probability density plot is given for the first column of input data. And a fitted line is added for the efficient closed-form estimator. For \code{which=4}, is a marginally fitted probability density plot is given like the former one for the second column of input data. Note that, marginally fitted probability density plots in \code{which=3} and \code{which=4} present comparisons between efficient closed form estimators (MLEces) and correlation based method estimators (CMEs) for the bivariate Weibull distribution. Note that this \code{plot} commend is limited at the bivariate distributions.
#' @import ggplot2 reshape
#' @importFrom graphics boxplot 
#' @importFrom LaplacesDemon ddirichlet
#' @importFrom graphics axis hist lines
#' @examples
#' data(flood)
#' est_BiGam <- MLEce(flood, "BiGam")
#' plot(est_BiGam, c(3))
#' @examples
#' air_data <- airquality[ ,3:4]
#' air_data[ ,2] <- air_data[ ,2]*0.1
#' est_BiWei <- MLEce(air_data, "BiWei")
#' plot(est_BiWei)
#' @examples
#' data(fossil_pollen) 
#' fossil_data <- cbind(fossil_pollen[,1]/100,rowSums(fossil_pollen[,-1]/100))
#' est_fossil <- MLEce(fossil_data, "Dirichlet")
#' plot(est_fossil,c(2))
#' @export
plot.MLEce <- function(x, which=c(1,2,3,4),
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...){
    show <- rep(FALSE, 4)

    dat <- as.matrix(x$data)
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
    
    if(x$distribution == "BiGam"){
      ce_est <- x$estimation
      
      mme = BiGam_MME(dat)
      initial.param <- c(mme$alpha1,mme$alpha2,mme$beta )
      est.mle = BiGam_MLE(dat, initial.param, tol2=mme)
      
      if (is.null(colnames(dat))){xy.name = c("x", "y")}
      else xy.name = colnames(dat)
      
      if(show[1]==TRUE){
        colnames(dat) = c('variable', 'value')
        boxplot(dat, xlab=c('variables'), ylab=c('values'), main=c('Boxplots for Data (Bivariate Gamma)') )
      
      }
      
      if(show[2]==TRUE){
        ##Contour plots.
        gr.x=seq(from=min(dat[,1]*0.9),to=max(dat[,1]*1.1),length.out = 250)
        gr.y=seq(from=min(dat[,2]*0.9),to=max(dat[,2]*1.1),length.out = 250)
        comdata = expand.grid(gr.x, gr.y)
        index = which(comdata[,1]>comdata[,2])
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
    if(x$distribution == "BiWei"){
        cme <- BiWei_CME(x$data)
        ce_est <- x$estimation

        if (is.null(colnames(dat))){xy.name = c("x", "y")}
        else xy.name = colnames(dat)

        if(show[1]==TRUE){
          colnames(dat) = c('variable', 'value')
          boxplot(dat, xlab=c('variables'), ylab=c('values'), main=c('Boxplots for Data (Bivariate Weibull)') )
          
        }

        if(show[2]==TRUE){
            gr.x = seq(from=min(dat[,1]) * 0.9,to=max(dat[,1])*1.1,length.out = 250)
            gr.y = seq(from=min(dat[,2]) * 0.9,to=max(dat[,2])*1.1,length.out = 250)
            pars = expand.grid(gr.x, gr.y)

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
            legend("topleft", c("CME","MLEce"), col=c(2,3), lwd=2,lty=c(3,2), cex=1.1)

        }

        if(show[4]==TRUE){
            main.title = paste0("Marginal plot of ", xy.name[2])
            gr.x = seq(from=min(dat[,2]) * 0.9,to=max(dat[,2])*1.1,length.out = 250)
            plot(gr.x, dweibull(gr.x,shape=cme[3],scale=cme[4]),type="l",lty=3,col=2,xlab=xy.name[2],ylab="", main = main.title)
            points(gr.x,dweibull(gr.x,shape=ce_est[3],scale=ce_est[4]),type="l",lty=2,xlab="", ylab = "", col=3)
            rug(x$data[,2])
            legend("topleft", c("CME","MLEce"), col=c(2,3), lwd=2,lty=c(3,2), cex=1.1)
        }
    }
    

    if(x$distribution == "Dirichlet"){

        ce_est <- x$estimation
        

        if (is.null(colnames(dat))){xy.name = c("x", "y")}
        else xy.name = colnames(dat)

        if(show[1]==TRUE){
          colnames(dat) = c('variable', 'value')
          boxplot(dat, xlab=c('variables'), ylab=c('values'), main=c('Boxplots for Data (Dirichlet)') )
          
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





#' Generating random data for the bivariate gamma distribution with parameters.
#'
#' @details
#' Random generation for the bivariate gamma distribution is presented. The specific generation formulas can be found  in Jang, et al. (2020). 
#' @param n number of observations.
#' @param paras parameters of bivariate gamma distribution (shape1, shape2, scale).
#' @references Jang, Y.-H., Zhao, J., Kim, H.-M., Yu, K., Kwon, S.and Kim, S. (2023) New closed-form efficient estimator for the multivariate gamma distribution. \emph{Statistica Neerlandica}, 1–18.
#' @return \code{rBiGam} generates random deviates. The length of generated data is determined by \code{"n"}.
#' @examples
#' datt = rBiGam(n=50, c(4,3,3))
#' @export
rBiGam = function(n, paras){
    V1 = rgamma(n,shape = paras[1],scale=paras[3])
    V2 = rgamma(n,shape = paras[2],scale=paras[3])
    rgmat = cbind(V1,V1+V2)
    colnames(rgmat) = c('x','y')
    return(rgmat)
}


#' Generating random data for the bivariate Weibull distribution.
#'
#' @details
#' \code{rBiWei} generates random number data for bivariate Weibull distribution.
#' @param paras parameters of bivariate Weibull distribution (alpha1, beta1, alpha2, beta2, delta).
#' @param n number of observations.
#' @return \code{rBiWei} generates random deviates.The length of generated data is determined by \code{"n"}
#' @examples
#' datt = rBiWei(n=50, c(4,3,3,4,0.6))
#' @importFrom nleqslv nleqslv
#' @import stats
#' @export
rBiWei = function(n,paras){
    a=paras[1];b=paras[2];aa=paras[3];bb=paras[4];del=paras[5]
    y = rweibull(n,aa,bb)
    bstar = 1+del*(y/bb)^aa
    p = 1-del/bstar
    x.index = apply(matrix(c(p,1-p),ncol=2),1,function(pp){sample(1:2,1,replace = T,prob=pp)})
    x = (b*(rgamma(n,1,1)/bstar)^(1/a))*(x.index==1)+(b*(rgamma(n,2,1)/bstar)^(1/a))*(x.index==2)
    return(cbind(x,y))
}

#' The flood events data of the Madawaska basin.
#'
#' The data is a subset of the flood events data of the Madawaska basin which is located in the province of Qebec, Canada. Daily streamflow data from 1919 to 1995 are available from HYDAT CD (1998) and Yue (2001).
#' @docType data
#' @usage data(flood, package = "MLEce")
#' @format A dataframe with 2 variables and 77 observations as follows:
#' \describe{
#' \item{Vnorm}{daily average flood volume}
#' \item{Q}{ flood peak}
#' }
#' @references Environment Canada (1998), HYDAT CD-ROM Version 98-1.05.8: Surface water and sediment data.
#' @references Yue. S. (2001) A bivariate gamma distribution for use in multivariate flood frequency analysis. \emph{Hydrological Processes}, 15, 1033–1045.
"flood"



#' The counts data of the frequency of occurrence of different kinds of fossil pollen grains.
#' 
#' The data is about the counts of the frequency of occurrence of different kinds of fossil pollen grains and is available in Mosimann (1962).
#'
#' @docType data
#' @usage data(fossil_pollen, package = "MLEce")
#' @format A dataframe with 73 observations and 4 variables: pinus, abies, quercus and alnus pollens.
#' @references Mosimann, J. E. (1962) On the compound multinomial distribution, the multivariate beta distribution, and correlations among proportions. \emph{Biometrika}, 49, 65-82.
"fossil_pollen"



#Some function required in former
#-------------------------------------------------------------------------------------------------
#Bivariate gamma distribution------------------------------------------------------
#-------------------------------------------------------------------------------------------------

BiGam_CE = function( data ){
  
  if (all(data[,1] > data[,2]) | all(data[,1]==data[,2]) ){ stop("In bivariate gamma distribution, 
     the first vector must be smaller than the second vector") }
  
  dat1 = data[,1] ; dat2 = data[,2]
  MMEest = BiGam_MME(data)
  a1 = MMEest$alpha1 ; a2 = MMEest$alpha2; b = MMEest$beta ; n=length(dat1)
  
  di.a1 = digamma(a1);      di.a2 = digamma(a2)
  log.b = log(b);           sum_data2 <- sum(dat2)
  tri.a1 = trigamma(a1);    tri.a2 = trigamma(a2)
  
  
  l1 = -n*di.a1  -n*log.b + sum(log(dat1))
  l2 = -n*di.a2  -n*log.b + sum(log(dat2-dat1))
  l3 = -n*(a1+a2)/b + sum_data2/b^2
  
  l11 = -n*tri.a1
  l22 = -n*tri.a2
  l33 = n*(a1+a2)/b^2-2*sum_data2/b^3
  
  l12 = l21 = 0
  l13 = l31 = l23 = l32 = -n/b
  
  J = matrix(c(l11,l12,l13,  l21,l22,l23,
               l31,l32,l33),ncol=3,byrow=T)
  Score = c(l1,l2,l3)
  
  
  J_log = matrix(c(l11*a1^2+l1*a1, l12*a1*a2     , l13*a1*b,
                   l21*a2*a1     , l22*a2^2+l2*a2, l23*a2*b,
                   l31*b*a1      , l32*b*a2      , l33*b^2+l3*b),ncol=3,byrow=T)
  Score_log = c(l1*a1,l2*a2,l3*b)
  pars_current = log(c(a1,a2,b))
  
  OBS_I = -qr.solve(J_log)
  otp = exp(c(pars_current+OBS_I%*%Score_log))
  names(otp) = c("alpha1","alpha2","beta")
  est_results = list(estimation = otp,  Hessian= J, Score=Score)
  return(est_results)
  
}


BiGam_MME = function( data ){
  
  if (!( (is.matrix(data) & is.numeric(data)) | is.data.frame(data) ) )
  { stop("data must be a numeric matrix or dataframe") }
  
  if (all(data[,1] > data[,2]) | all(data[,1]==data[,2]) ){ stop("In bivariate gamma distribution, 
     the first vector must be smaller than the second vector") }
  
  dat1 = data[,1] ; dat2 = data[,2]
  mean_dat1 <- mean(dat1); mean_dat2 <- mean(dat2)
  a1 = mean_dat1^2/mean((dat1-mean_dat1)^2)
  a2 = mean_dat1*(mean_dat2-mean_dat1)/mean((dat1-mean_dat1)^2)
  b1 = mean((dat1-mean(dat1))^2)/mean_dat1
  b2 = mean((dat2-mean_dat2)^2)/mean_dat2
  result = list(alpha1=a1,alpha2=a2,beta=b1)
  return(result)
  
}



BiGam_MLE = function(data,initial.param = NULL,tol=10e-6,fail_num=50,tol2=c(1.2,1.2,0.45),re_sd=diag(c(0.1,0.1,0.01))){
  if(is.null(initial.param)){
    Ini_MME <- BiGam_MME(data)
    initial.param <- c(Ini_MME$alpha1,Ini_MME$alpha2,Ini_MME$beta)
  }
  starting=initial.param
  dist.tol=10 ; start_save = log(initial.param) ; re_iter =0
  while(dist.tol>tol){
    tmpval = log(BiGam_CE(data)$estimation)
    distt = abs(tmpval-start_save)#abs(exp(tmpval)-exp(start_save))
    dist.tol = max(distt)
    start_save = tmpval
    if(sum(distt>tol2)>0){
      start_save = log(starting + abs(c(mvtnorm::rmvnorm(1,c(0,0,0),re_sd) ) ))
      re_iter = re_iter+1
    }
    if(re_iter>=fail_num){
      return(c(0,0,0))
    }
  }
  return(exp(tmpval))
}
#log-likelihood function for the bivariate gamma distribution
dBiGam = function(pars,dat1,dat2,log=TRUE){
  a1=pars[1];a2=pars[2];b=pars[3]
  result =  -lgamma(a1)-lgamma(a2)-(a1+a2)*log(b)+(a1-1)*log(dat1)+(a2-1)*log(dat2-dat1)-dat2/b
  if(log==TRUE){return(result)}
  if(log==FALSE){return(exp(result))}
}

#-------------------------------------------------------------------------------------------------
#Bivariate Weibull distribution------------------------------------------------------
#-------------------------------------------------------------------------------------------------

BiWei_CE = function(data){
  
  #Starting values
  marginal_old = BiWei_CME(data)
  
  #Estimation
  result = BiWei_info(marginal_old,dat=data,type="MLECE")
  est = result$estimation
  names(est) = c("alpha1","beta1","alpha2","beta2","delta")
  estls = list(estimation = est, cme = marginal_old)
  return(estls)
}

delta_score_probit = function(par_vec,dat){
  del = par_vec[5] ; a1 = par_vec[1] ; b1 = par_vec[2] ; a2 = par_vec[3] ; b2 = par_vec[4]
  dat1 = dat[,1] ; dat2 = dat[,2]
  n = length(dat1)
  tmp1 = (dat1/b1)^a1 ; tmp2 = (dat2/b2)^a2
  dg1 = (1+del*tmp1) ;  dg2 = (1+del*tmp2)
  g = dg1*dg2-del ; g5 = tmp1*dg2 + tmp2*dg1 - 1
  return(  (-sum(tmp1*tmp2) + sum(g5/g))*dnorm(qnorm(par_vec[5])) )
}

BiWei_info = function(par_vec,dat,type){ #par_vec in original scale.
  a1 = par_vec[1] ; b1 = par_vec[2]  ;a2 = par_vec[3] ; b2 = par_vec[4]
  del = par_vec[5] ;
  a1Db1 = a1/b1;  a2Db2 = a2/b2 
  dat1 = dat[,1] ; dat2 = dat[,2] ; n=length(dat1)
  
  tmp1 = (dat1/b1)^a1 ; tmp2 = (dat2/b2)^a2
  ltmp1 = log(dat1/b1) ; ltmp2 = log(dat2/b2)
  h1 = tmp1*ltmp1 ; h2 = tmp2*ltmp2 ; ha1 = h1*ltmp1 ; ha2 = h2*ltmp2
  hb1 = -(1/b1)*tmp1*(a1*ltmp1+1) ; hb2 = -(1/b2)*tmp2*(a2*ltmp2+1)
  dg1 = (1+del*tmp1) ; dg2 = (1+del*tmp2)
  
  g = dg1*dg2-del ; g1 = del*h1*dg2  ;g3 = del*h2*dg1
  g2 = -del*tmp1*dg2*(a1Db1) ; g4 = -del*tmp2*dg1*(a2Db2) ; g5 = tmp1*dg2 + tmp2*dg1 - 1
  
  g11 = del*ha1*dg2 ; g33 = del*ha2*dg1 ; g22 = g2*(-a1-1)/b1 ;
  g44 = g4*(-a2-1)/b2 ; g55 = 2*tmp1*tmp2
  
  g13 = g31 = del^2*h1*h2            ; g24 = g42 = del^2*a1Db1*a2Db2*tmp1*tmp2
  g12 = g21 = del*dg2*hb1            ; g34 = g43 = del*dg1*hb2
  g14 = g41 = -del^2*h1*(a2/b2)*tmp2 ; g23 = g32 = -del^2*h2*(a1Db1)*tmp1
  g15 = g51 = h1*(1+2*del*tmp2)      ; g35 = g53 = h2*(1+2*del*tmp1)
  g25 = g52 = -(a1Db1)*tmp1*(1+2*del*tmp2) ; g45 = g54 = -(a2Db2)*tmp2*(1+2*del*tmp1)
  
  l1 = (-sum(h1)-del*sum(h1*tmp2)+ n/a1 + sum(ltmp1)+sum(g1/g))*a1
  l3 = (-sum(h2)-del*sum(h2*tmp1)+ n/a2 + sum(ltmp2)+sum(g3/g))*a2
  l2 = ((a1Db1)*sum(tmp1) + del*(a1Db1)*sum(tmp1*tmp2)-n*(a1Db1)+sum(g2/g))*b1
  l4 = ((a2Db2)*sum(tmp2) + del*(a2Db2)*sum(tmp2*tmp1)-n*(a2Db2)+sum(g4/g))*b2
  l5 = (-sum(tmp1*tmp2) + sum(g5/g))*dnorm(qnorm(par_vec[5]))
  score = c(l1,l2,l3,l4,l5)
  if(type=="score"){ return(score)}
  
  l11 = (-sum(ha1)-del*sum(tmp2*ha1)-n/a1^2+sum((g11*g-g1^2)/g^2) )*(a1^2) + l1
  l33 = (-sum(ha2)-del*sum(tmp1*ha2)-n/a2^2+sum((g33*g-g3^2)/g^2) )*(a2^2) + l3
  l22 = (-((a1+1)/b1)*(a1Db1)*sum(tmp1) - del*((a1+1)/b1)*(a1Db1)*sum(tmp1*tmp2)+n*a1/(b1^2) + sum((g22*g-g2^2)/g^2))*(b1^2) + l2
  l44 = (-((a2+1)/b2)*(a2Db2)*sum(tmp2) - del*((a2+1)/b2)*(a2Db2)*sum(tmp2*tmp1)+n*a2/(b2^2) + sum((g44*g-g4^2)/g^2))*(b2^2) + l4
  l55 = (sum((g55*g-g5^2)/g^2))*dnorm(qnorm(par_vec[5]))^2+l5*dnorm(qnorm(par_vec[5]))*(-qnorm(par_vec[5]))#(1-2*del)
  if(type=="hessian_del"){ return(hessian=l55)}
  
  l13 = l31 = (-del*sum(h1*h2) + sum((g13*g-g1*g3)/g^2))*a1*a2
  l24 = l42 = (-del*(a1Db1*a2Db2)*sum(tmp1*tmp2)+ sum( (g24*g-g2*g4)/g^2 ))*b1*b2
  l12 = l21 = (-sum(hb1) -del*sum(hb1*tmp2)-n/b1 + sum( (g12*g-g1*g2)/g^2 ))*a1*b1
  l34 = l43 = (-sum(hb2) -del*sum(hb2*tmp1)-n/b2 + sum( (g34*g-g3*g4)/g^2 ))*a2*b2
  l14 = l41 = (del*(a2Db2)*sum(h1*tmp2)+sum( (g14*g-g1*g4)/g^2  ))*a1*b2
  l23 = l32 = (del*(a1Db1)*sum(h2*tmp1)+sum( (g23*g-g2*g3)/g^2 ))*a2*b1
  l15 = l51 = (-sum(h1*tmp2) + sum( (g15*g-g1*g5)/g^2))*dnorm(qnorm(par_vec[5]))*a1
  l35 = l53 = (-sum(h2*tmp1) + sum( (g35*g-g3*g5)/g^2))*dnorm(qnorm(par_vec[5]))*a2
  l25 = l52 = ((a1Db1)*sum(tmp1*tmp2) +   sum( (g25*g-g2*g5)/g^2 ))*dnorm(qnorm(par_vec[5]))*b1
  l45 = l54 = ((a2Db2)*sum(tmp2*tmp1) +   sum( (g45*g-g4*g5)/g^2 ))*dnorm(qnorm(par_vec[5]))*b2
  
  J = matrix(c(l11,l12,l13,l14,l15,
               l21,l22,l23,l24,l25,
               l31,l32,l33,l34,l35,
               l41,l42,l43,l44,l45,
               l51,l52,l53,l54,l55),ncol=5,byrow=T)
  if(type=="mar"){return(list(score[1:4],J[1:4,1:4]))}
  if(type=="del"){return(list(score[5],J[5,5]))}
  if(type=="hessian"){return(list(score=score,hessian=J))}
  
  par_vec_tmp = c(log(par_vec[1:4]),qnorm(par_vec[5]))
  par_vec_tmp = par_vec_tmp - score%*%solve(J)
  
  est = c(exp(par_vec_tmp[1:4]) , pnorm(par_vec_tmp[5]))
  
  if(type=="MLECE")return(list(estimation = est, hessian = J, score = score))
}

BiWei_CMEpar = function(dat){
  n = length(dat)
  alphahat = -log(2)/log(1-sqrt((n+1)/(n-1)/3)*(sd(dat)*sqrt((n-1)/(n))/mean(dat))*cor(dat,rank(dat)))
  return( c(alphahat,(mean(dat^alphahat))^(1/alphahat)) )
}

BiWei_CME = function(data){
  estpar1 = c(BiWei_CMEpar(data[,1]),BiWei_CMEpar(data[,2]))
  estpar2 = pnorm( nleqslv::nleqslv(0,function(delta){delta_score_probit(c(estpar1,delta),data)})$x)
  return( c(estpar1,estpar2) )
}

BiWei_MLE=function(data,initial.param=NULL,tol=1e-7){
  if(is.null(initial.param)){
    initial.param <- BiWei_CME(data)
  }
  
  tmp_par = initial.param
  dist=100
  while(dist>tol){
    
    par = BiWei_info(tmp_par,data,type = "hessian")
    new_par = c( c(log(tmp_par[1:4]),qnorm(tmp_par[5])) - par[[1]]%*%qr.solve(par[[2]]))
    new_par = c( exp(new_par[1:4]), pnorm(new_par[5]))
    
    if(new_par[5]>0.9999){new_par[5]=0.9999}
    if(new_par[5]<0.0001){new_par[5]=0.0001}
    dist = max(abs( tmp_par-new_par))
    tmp_par = new_par
  }
  return(c(new_par))
}

#Evaluating log-likelihood value of bivariate Weibull distribution of Gumbel-type
dBiWei = function(par_vec,dat,log=TRUE){
  a=par_vec[1];b=par_vec[2];aa=par_vec[3];bb=par_vec[4];del=par_vec[5]
  x=dat[,1];y=dat[,2]
  n=length(x); xDb <- x/b; yDbb <- y/bb
  result =  n*log(a) -n*log(b) +(a-1)*sum(log(xDb))+
    n*log(aa)-n*log(bb)+(aa-1)*sum(log(yDbb))-
    sum((xDb)^a) - sum((yDbb)^aa) - del*sum( ((xDb)^a)*((yDbb)^aa) )+
    sum(log(  (1+del*(xDb)^a)*(1+del*(yDbb)^aa) -del  ))
  if(log==TRUE){
    return( c(result))
  }else{
    return(c(exp(result)))
  }
}
#-------------------------------------------------------------------------------------------------
#Dirichlet distribution------------------------------------------------------
#-------------------------------------------------------------------------------------------------

Diri_MME <- function(x) {  # calculate MME using data
  m1 <- colMeans(x)
  m2 <- colMeans(x ^ 2)
  return(matrix((m1 * (m1 - m2)) / (m2 - m1 ^ 2) , nrow = 1))
}

Diri_MLE <- function(x, eps = 1e-10, mxit = 1e5) {  # calculate ML estimator using data
  return(matrix(dirichlet.mle(x, eps = eps, maxit = mxit)$alpha ,nrow = 1))  # outputs estimated parameters in a row vector
}

#a function needed to calculate the parametric confidence interval.
Diri_CE_bt <- function(x) {
  m <- ncol(x)
  tldb <- log(Diri_MME(x))
  exp_tldb <- exp(tldb)
  # constants
  c1 <- trigamma(c(sum(exp_tldb), exp_tldb))
  c2 <- digamma(c(sum(exp_tldb), exp_tldb))
  
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
      nrow(x) * exp_tldb[i] * exp_tldb[j] * c1[1]
    }
  }
  D2ll <- Vectorize(D2ll)
  iH <- solve(outer(1:m, 1:m, D2ll))   # inverse Hessian
  ans <- as.vector(tldb - D1ll %*% iH)  # estimation
  names(ans) <- c(paste("alpha", 1:m, sep=""))
  return(list(estimation=exp(ans), InvHess=iH) ) # back-transformation
}
Diri_CE <- function(x) {  # calculate closed-form estimator using data
  m <- ncol(x)
  MME <- Diri_MME(x)
  c1 <- trigamma(c(sum(MME), MME))
  denom <- 1 / c1[1] - sum(1 / c1[-1])
  c2 <- digamma(c(sum(MME), MME))
  iD <- diag(1 / c1[-1])
  iH <-
    iD %*% (diag(1, m) + matrix(1, nrow = m, ncol = m) %*% iD / denom)
  grad <- as.matrix(colMeans(log(x)) - (c2[-1] - c2[1]),ncol=1 )
  ans <- matrix(MME + as.vector(iH %*% grad),nrow=1)
  otp = as.vector(pmax(ans, 1 / nrow(x)))
  names(otp) <- c(paste("alpha", 1:m, sep=""))
  estls = list(estimation = otp)
  return(estls)
}
#----------------------------------------------------------------------------------
# GCVM goodness of fit test---------------------------------------------------
#----------------------------------------------------------------------------------
# Rosen's transformation for GCVM gof test
Rosen = function(dat,pars){
  z1 = pweibull(dat[,2],pars[3:4])
  betastar = (1+pars[5]*(dat[,2]/pars[4])^pars[3])
  p = 1-pars[5]/betastar
  z2=pexp( ((dat[,1]/pars[2])^pars[1])*betastar  )*p +
    pgamma( ((dat[,1]/pars[2])^pars[1])*betastar ,shape=2 ,scale=1 )*(1-p)
  
  z11 = pweibull(dat[,1],pars[1:2])
  betastar = (1+pars[5]*(dat[,1]/pars[2])^pars[1])
  p = 1-pars[5]/betastar
  z22=pexp( ((dat[,2]/pars[4])^pars[3])*betastar  )*p +
    pgamma( ((dat[,2]/pars[4])^pars[3])*betastar ,shape=2 ,scale=1 )*(1-p)
  
  return(list(fisrt = cbind(z1,z2), second=cbind(z11,z22) ))
}

CD2 = function(dat){
  grd = expand.grid(1:nrow(dat),1:nrow(dat))
  xx1 = dat[grd[,1],]
  xx2 = dat[grd[,2],]
  grdtmp = 0;
  return(
    sqrt((13/12)^2-
           2*sum( ( 1+0.5*abs(dat[,1]-0.5)-0.5*abs(dat[,1]-0.5)^2 )*
                    ( 1+0.5*abs(dat[,2]-0.5)-0.5*abs(dat[,2]-0.5)^2 ) )/nrow(dat)+
           sum( (1+0.5*abs(xx1[,1]-0.5)+0.5*abs(xx2[,1]-0.5)-0.5*abs(xx1[,1]-xx2[,1]))*
                  (1+0.5*abs(xx1[,2]-0.5)+0.5*abs(xx2[,2]-0.5)-0.5*abs(xx1[,2]-xx2[,2])))/nrow(dat)^2
    ))
}



# SD2 statistics for GCVM gof test
SD2fun.gam <- function(Data, EST){
  N <- dim(Data)[1]
  
  #Rosenblatt Transformation
  #marginal distribution of #Vnorm~Gamma(alpha_1,beta)
  xdim1 <- pgamma(Data[,1], shape=EST[1],scale=EST[3])
  #conditional dist:Q|Vnorm~Gamma(alpha_2,beta,Vnorm),here Vnorm is the location parameter value
  xdim2minus <- pgamma(Data[,1], shape=EST[2],scale=EST[3])
  xdim2 <- pgamma(Data[,2]+Data[,1], shape=EST[2],scale=EST[3])-xdim2minus
  TransData <- cbind(xdim1,xdim2)
  
  SD2par1 <- (4/3)^2
  SD2par2 <- 2*sum( (1+2*xdim1-2*xdim1^2)*(1+2*xdim2-2*xdim2^2) )/N
  grd <- expand.grid(1:N, 1:N)
  xx1 <- TransData[grd[,1], ]
  xx2 <- TransData[grd[,2], ]
  SD2par3 <- (1-abs(xx1[,1]-xx2[,1]))*(1-abs(xx1[,2]-xx2[,2]))
  
  SD2 <- sqrt(SD2par1-SD2par2 +4*sum(SD2par3)/(N^2) )
  
  return(Statis=SD2)
}

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






