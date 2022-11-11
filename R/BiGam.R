#' Get closed-form estimator for Bivariate gamma
#' @param pars parameters of bivariate gamma (alpha1, alpha2, beta).
#' @param dat data of bivariate gamma
#' @param type output type (MLECE, hessian, score). Default is MLECE
#' @param log log-transformation of data. Default is TRUE
BiGam_CE = function(pars,dat,type="MLECE",log=TRUE){
    dat1 = dat[,1] ; dat2 = dat[,2]
    a1 = pars[1] ; a2 = pars[2]; b = pars[3] ; n=length(dat1)
    di.a1 = digamma(a1)
    di.a2 = digamma(a2)
    log.b = log(b)
    tri.a1 = trigamma(a1)
    tri.a2 = trigamma(a2)

    l1 = -n*di.a1  -n*log.b + sum(log(dat1))
    l2 = -n*di.a2  -n*log.b + sum(log(dat2-dat1))
    l3 = -n*(a1+a2)/b + sum(dat2)/b^2

    l11 = -n*tri.a1
    l22 = -n*tri.a2
    l33 = n*(a1+a2)/b^2-2*sum(dat2)/b^3

    l12 = l21 = 0
    l13 = l31 = l23 = l32 = -n/b


    if(log==FALSE){
        J = matrix(c(l11,l12,l13,
                     l21,l22,l23,
                     l31,l32,l33),ncol=3,byrow=T)
        Score = c(l1,l2,l3)
        pars_current = pars
    }else{
        J = matrix(c(l11*a1^2+l1*a1, l12*a1*a2     , l13*a1*b,
                     l21*a2*a1     , l22*a2^2+l2*a2, l23*a2*b,
                     l31*b*a1      , l32*b*a2      , l33*b^2+l3*b),ncol=3,byrow=T)
        Score = c(l1*a1,l2*a2,l3*b)
        pars_current = log(pars)
    }
    if(type=="MLECE"){
        OBS_I = solve(-J)
        if(log==FALSE){
            otp = c(pars_current+OBS_I%*%Score)
        }else{
            otp = exp(c(pars_current+OBS_I%*%Score))
        }
        names(otp) = c("alpha1","alpha2","beta")
        estls = list(estimation = otp, hess = J)
        return(estls)
    }
    if(type=="hessian"){
        return(J)
    }
    if(type=="score"){
        return(Score)
    }
}


#' log-likelihood for Bivariate gamma
#' @param pars parameters of bivariate gamma (alpha1, alpha2, beta).
#' @param dat1 data of marginal column (univariate gamma).
#' @param dat2 other data of marginal column (univariate gamma).
#' @param log log-transformation of data. Default is TRUE
dBiGam = function(pars,dat1,dat2,log=TRUE){
    a1=pars[1];a2=pars[2];b=pars[3]
    result =  -log(gamma(a1))-log(gamma(a2))-(a1+a2)*log(b)+(a1-1)*log(dat1)+(a2-1)*log(dat2-dat1)-dat2/b
    if(log==TRUE){return(result)}
    if(log==FALSE){return(exp(result))}
}

#' random generation for the Bivariate Gamma distribution with (shape1, shape2, scale).
#'
#' @details
#' These functions implement formulas given in Hyoung-Moon Kim. et al. (2020). (will be revised.)
#'
#' @param pars parameters of BiWeibull (shape1, shape2, scale).
#' @param n number of observations.
#'
#' @return \code{rBiGam} generates random deviates.
#'
#' @examples
#' datt = rBiGam(n=50, c(4,3,3))
#'
#' @export
rBiGam = function(n, pars){
    V1 = rgamma(n,shape = pars[1],scale=pars[3])
    V2 = rgamma(n,shape = pars[2],scale=pars[3])
    rgmat = cbind(V1,V1+V2)
    colnames(rgmat) = c('x','y')
    return(rgmat)
}

#' MME for Bivariate gamma
#' @param dat data of bivariate gamma
#' @param scaletype scale type for bivariate gamma MME
BiGam_MME = function(dat,scaletype="first"){
    if(is.matrix(dat)==F & is.data.frame(dat)==F){stop("data should be in matrix.")}

    dat1 = dat[,1] ; dat2 = dat[,2]
    a1 = mean(dat1)^2/mean((dat1-mean(dat1))^2)
    a2 = mean(dat1)*(mean(dat2)-mean(dat1))/mean((dat1-mean(dat1))^2)
    b1 = mean((dat1-mean(dat1))^2)/mean(dat1)
    b2 = mean((dat2-mean(dat2))^2)/mean(dat2)
    result = c(a1,a2,b1,b2)
    names(result) = c("shape1","shape2","scale_ver1","scale_ver2")
    if(scaletype=="first"){return(result[1:3])}
    if(scaletype=="second"){return(result[c(1,2,4)])}
    if(scaletype=="both"){return(result)}
}

#' @import mvtnorm
BiGam_MLE = function(starting,dat,tol=10e-6,fail_num=50,tol2=c(1.2,1.2,0.45),re_sd=diag(c(0.1,0.1,0.01))){
    dist.tol=10 ; start_save = starting ; re_iter =0
    while(dist.tol>tol){
        tmpval = BiGam_CE(starting,dat)$estimation
        distt = abs(tmpval-starting) ; dist.tol = max(distt)
        starting = tmpval
        if(sum(distt>tol2)>0){
            starting = start_save +c(rmvnorm(1,c(0,0,0),re_sd) )
            re_iter = re_iter+1
        }
        if(re_iter>=fail_num){
            return(c(0,0,0))
        }
    }
    return(tmpval)
}


#' @import mvtnorm
BiGam_MLE_log = function(starting,dat,tol=10e-6,fail_num=50,tol2=c(1.2,1.2,0.45),re_sd=diag(c(0.1,0.1,0.01))){
    dist.tol=10 ; start_save = log(starting) ; re_iter =0
    while(dist.tol>tol){
        tmpval = log(BiGam_CE(exp(start_save),dat,log = TRUE)$estimation)
        distt = abs(tmpval-start_save)#abs(exp(tmpval)-exp(start_save))
        dist.tol = max(distt)
        start_save = tmpval
        if(sum(distt>tol2)>0){
            start_save = log(starting + abs(c(rmvnorm(1,c(0,0,0),re_sd) ) ))
            re_iter = re_iter+1
        }
        if(re_iter>=fail_num){
            return(c(0,0,0))
        }
    }
    return(exp(tmpval))
}

#' SD2 statistics for GCVM gof test
#' @param Data data for gof test
#' @param EST estimates for gof test
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
