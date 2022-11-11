#' Generating random number for the bivariate Weibull distribution with (alpha1, beta1, alpha2, beta2, delta).
#'
#' @details
#' \code{rBiWei} generates random number data for bivariate weibull distribution.
#'
#' @param par_vec parameters of BiWeibull (alpha1, beta1, alpha2, beta2, delta).
#' @param n number of observations.
#'
#' @return \code{rBiWei} generates random deviates.
#'
#' @examples
#' datt = rBiWei(n=50, c(4,3,3,4,0.6))
#'
#' @importFrom nleqslv nleqslv
#' @import stats
#'
#' @export
rBiWei = function(n,par_vec){
    a=par_vec[1];b=par_vec[2];aa=par_vec[3];bb=par_vec[4];del=par_vec[5]
    y = rweibull(n,aa,bb)
    bstar = 1+del*(y/bb)^aa
    p = 1-del/bstar
    x.index = apply(matrix(c(p,1-p),ncol=2),1,function(pp){sample(1:2,1,replace = T,prob=pp)})
    x = (b*(rgamma(n,1,1)/bstar)^(1/a))*(x.index==1)+(b*(rgamma(n,2,1)/bstar)^(1/a))*(x.index==2)
    return(cbind(x,y))
}

#' Evaluating bivariate Weibull distribution of Gumbel-type
#' @param par_vec parameters of bivariate weibull (alpha1, beta1, alpha2, beta2, delta).
#' @param dat data of bivariate weibull
#' @param log log-transformation for data. Default is TRUE.
dBiWei = function(par_vec,dat,log=TRUE){
    a=par_vec[1];b=par_vec[2];aa=par_vec[3];bb=par_vec[4];del=par_vec[5]
    x=dat[,1];y=dat[,2]
    n=length(x)
    result =  n*log(a) -n*log(b) +(a-1)*sum(log(x/b))+
        n*log(aa)-n*log(bb)+(aa-1)*sum(log(y/bb))-
        sum((x/b)^a) - sum((y/bb)^aa) - del*sum( ((x/b)^a)*((y/bb)^aa) )+
        sum(log(  (1+del*(x/b)^a)*(1+del*(y/bb)^aa) -del  ))
    if(log==TRUE){
        return( c(result))
    }else{
        return(c(exp(result)))
    }
}

#' Get root-n consistent estimator by correlation method.
#' @param dat data of bivariate weibull
cor_method = function(dat){
    n = length(dat)
    alphahat = -log(2)/log(1-sqrt((n+1)/(n-1)/3)*(sd(dat)*sqrt((n-1)/(n))/mean(dat))*cor(dat,rank(dat)))
    return( c(alphahat,(mean(dat^alphahat))^(1/alphahat)) )
}

#' Calculating delta score
#' @param par_vec parameters of bivariate weibull (alpha1, beta1, alpha2, beta2, delta).
#' @param dat data of bivariate weibull
delta_score_probit = function(par_vec,dat){
    del = par_vec[5] ; a1 = par_vec[1] ; b1 = par_vec[2] ; a2 = par_vec[3] ; b2 = par_vec[4]
    dat1 = dat[,1] ; dat2 = dat[,2]
    n = length(dat1)
    tmp1 = (dat1/b1)^a1 ; tmp2 = (dat2/b2)^a2
    dg1 = (1+del*tmp1) ;  dg2 = (1+del*tmp2)
    g = dg1*dg2-del ; g5 = tmp1*dg2 + tmp2*dg1 - 1
    return(  (-sum(tmp1*tmp2) + sum(g5/g))*dnorm(qnorm(par_vec[5])) )
}

#' Calculating MLEce for Bivariate weibull
#' @param par_vec parameters of bivariate weibull (alpha1, beta1, alpha2, beta2, delta).
#' @param dat data of bivariate weibull
#' @param type output type (hessian, MLEce, del, mar)
BiWei_info = function(par_vec,dat,type){ #par_vec in original scale.
    a1 = par_vec[1] ; b1 = par_vec[2]  ;a2 = par_vec[3] ; b2 = par_vec[4]
    del = par_vec[5] ;

    dat1 = dat[,1] ; dat2 = dat[,2] ; n=length(dat1)

    tmp1 = (dat1/b1)^a1 ; tmp2 = (dat2/b2)^a2
    ltmp1 = log(dat1/b1) ; ltmp2 = log(dat2/b2)
    h1 = tmp1*ltmp1 ; h2 = tmp2*ltmp2 ; ha1 = h1*ltmp1 ; ha2 = h2*ltmp2
    hb1 = -(1/b1)*tmp1*(a1*ltmp1+1) ; hb2 = -(1/b2)*tmp2*(a2*ltmp2+1)
    dg1 = (1+del*tmp1) ; dg2 = (1+del*tmp2)

    g = dg1*dg2-del ; g1 = del*h1*dg2  ;g3 = del*h2*dg1
    g2 = -del*tmp1*dg2*(a1/b1) ; g4 = -del*tmp2*dg1*(a2/b2) ; g5 = tmp1*dg2 + tmp2*dg1 - 1

    g11 = del*ha1*dg2 ; g33 = del*ha2*dg1 ; g22 = g2*(-a1-1)/b1 ;
    g44 = g4*(-a2-1)/b2 ; g55 = 2*tmp1*tmp2

    g13 = g31 = del^2*h1*h2            ; g24 = g42 = del^2*(a1*a2)/(b1*b2)*tmp1*tmp2
    g12 = g21 = del*dg2*hb1            ; g34 = g43 = del*dg1*hb2
    g14 = g41 = -del^2*h1*(a2/b2)*tmp2 ; g23 = g32 = -del^2*h2*(a1/b1)*tmp1
    g15 = g51 = h1*(1+2*del*tmp2)      ; g35 = g53 = h2*(1+2*del*tmp1)
    g25 = g52 = -(a1/b1)*tmp1*(1+2*del*tmp2) ; g45 = g54 = -(a2/b2)*tmp2*(1+2*del*tmp1)

    l1 = (-sum(h1)-del*sum(h1*tmp2)+ n/a1 + sum(ltmp1)+sum(g1/g))*a1
    l3 = (-sum(h2)-del*sum(h2*tmp1)+ n/a2 + sum(ltmp2)+sum(g3/g))*a2
    l2 = ((a1/b1)*sum(tmp1) + del*(a1/b1)*sum(tmp1*tmp2)-n*(a1/b1)+sum(g2/g))*b1
    l4 = ((a2/b2)*sum(tmp2) + del*(a2/b2)*sum(tmp2*tmp1)-n*(a2/b2)+sum(g4/g))*b2
    l5 = (-sum(tmp1*tmp2) + sum(g5/g))*dnorm(qnorm(par_vec[5]))
    score = c(l1,l2,l3,l4,l5)
    if(type=="score"){ return(score)}

    l11 = (-sum(ha1)-del*sum(tmp2*ha1)-n/a1^2+sum((g11*g-g1^2)/g^2) )*(a1^2) + l1
    l33 = (-sum(ha2)-del*sum(tmp1*ha2)-n/a2^2+sum((g33*g-g3^2)/g^2) )*(a2^2) + l3
    l22 = (-((a1+1)/b1)*(a1/b1)*sum(tmp1) - del*((a1+1)/b1)*(a1/b1)*sum(tmp1*tmp2)+n*a1/(b1^2) + sum((g22*g-g2^2)/g^2))*(b1^2) + l2
    l44 = (-((a2+1)/b2)*(a2/b2)*sum(tmp2) - del*((a2+1)/b2)*(a2/b2)*sum(tmp2*tmp1)+n*a2/(b2^2) + sum((g44*g-g4^2)/g^2))*(b2^2) + l4
    l55 = (sum((g55*g-g5^2)/g^2))*dnorm(qnorm(par_vec[5]))^2+l5*dnorm(qnorm(par_vec[5]))*(-qnorm(par_vec[5]))#(1-2*del)
    if(type=="hessian_del"){ return(hessian=l55)}

    l13 = l31 = (-del*sum(h1*h2) + sum((g13*g-g1*g3)/g^2))*a1*a2
    l24 = l42 = (-del*(a1*a2/b1/b2)*sum(tmp1*tmp2)+ sum( (g24*g-g2*g4)/g^2 ))*b1*b2
    l12 = l21 = (-sum(hb1) -del*sum(hb1*tmp2)-n/b1 + sum( (g12*g-g1*g2)/g^2 ))*a1*b1
    l34 = l43 = (-sum(hb2) -del*sum(hb2*tmp1)-n/b2 + sum( (g34*g-g3*g4)/g^2 ))*a2*b2
    l14 = l41 = (del*(a2/b2)*sum(h1*tmp2)+sum( (g14*g-g1*g4)/g^2  ))*a1*b2
    l23 = l32 = (del*(a1/b1)*sum(h2*tmp1)+sum( (g23*g-g2*g3)/g^2 ))*a2*b1
    l15 = l51 = (-sum(h1*tmp2) + sum( (g15*g-g1*g5)/g^2))*dnorm(qnorm(par_vec[5]))*a1
    l35 = l53 = (-sum(h2*tmp1) + sum( (g35*g-g3*g5)/g^2))*dnorm(qnorm(par_vec[5]))*a2
    l25 = l52 = ((a1/b1)*sum(tmp1*tmp2) +   sum( (g25*g-g2*g5)/g^2 ))*dnorm(qnorm(par_vec[5]))*b1
    l45 = l54 = ((a2/b2)*sum(tmp2*tmp1) +   sum( (g45*g-g4*g5)/g^2 ))*dnorm(qnorm(par_vec[5]))*b2

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

#' Get root-n consistent closed-form estimator by correlation method and MLE.
#' @param data data of bivariate weibull
BiWei_CE = function(data){

    #Starting values
    marginal_old = c(cor_method(data[,1]),cor_method(data[,2]))
    delta_old = pnorm( nleqslv::nleqslv(0,function(delta){delta_score_probit(c(marginal_old,delta),data)})$x)

    #Estimation
    result = BiWei_info(c(marginal_old,delta_old),dat=data,type="MLECE")
    est = result$estimation
    names(est) = c("alpha1","beta1","alpha2","beta2","delta")
    estls = list(estimation = est, cme = c(marginal_old, delta_old))
    return(estls)
}

#' Get MLE for Bivariate weibull
#' @param data data of bivariate weibull
#' @param inits initial values of iterative algorithm for MLE
#' @param tol tolerance for difference.
BiWei_ML=function(data,inits,tol=1e-7){
    tmp_par = inits
    dist=100
    while(dist>tol){

        par = BiWei_info(tmp_par,data,type = "hessian")
        new_par = c( c(log(tmp_par[1:4]),qnorm(tmp_par[5])) - par[[1]]%*%solve(par[[2]]))
        new_par = c( exp(new_par[1:4]), pnorm(new_par[5]))

        if(new_par[5]>0.9999){new_par[5]=0.9999}
        if(new_par[5]<0.0001){new_par[5]=0.0001}
        dist = max(abs( tmp_par-new_par))
        tmp_par = new_par
    }
    return(c(new_par))
}


#' Rosen's transformation for GCVM gof test
#' @param dat data of bivariate weibull
#' @param pars parameters of bivariate weibull (alpha1, beta1, alpha2, beta2, delta).
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

#' CD2 statistics for GCVM gof test
#' @param dat data of bivariate weibull
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

