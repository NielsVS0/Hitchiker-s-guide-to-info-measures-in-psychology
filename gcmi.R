##############################################################
####### Gaussian Copula mutual information estimation   ######
####### Based on the python code by Ince et al. (2018)  ######
##############################################################

library(Rfast)
library(matrixStats)

atleast_2d <- function(x) {
  # x <- atleast_2d assures that the input is at least an array of dimension 2 (or higher)
  # R equivalent to the numpy fucntion atleast_2d
  
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  
  return(x)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

ctransform <- function(x){
  # Copula transformation (empirical CDF)
  # cx = ctransform(x) returns the empirical CDF value along the first
  # axis of x. Data is ranked and scaled within [0 1] (open interval).
  
  x <- atleast_2d(x)
  xr <- rowRanks(x,ties.method = "average")
  cx <- xr / (dim(xr)[length(dim(xr))]+1)
  
  return(cx)
  
}

copnorm <- function(x){
  # Copula normalization
  # cx = copnorm(x) returns standard normal samples with the same empirical
  # CDF value as the input. Operates along the last axis.
  
  cx <- qnorm(ctransform(x))
  
  return(cx)
}

ent_g <- function(x, biascorrect=TRUE) {
  # Entropy of a Gaussian variable in bits
  # H = ent_g(x) returns the entropy of a (possibly multidimensional)
  # Gaussian variable x with bias correction.
  # Columns of x correspond to samples, rows to dimensions/variables
  # (Samples last axis)
  
  x <- atleast_2d(x)
   
  if (length(dim(x)) > 2) {
    stop("x must be at most 2d")
  }
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  
  # demean the data
  x <- sweep(x,1,rowMeans(x))
  # covariance of the data
  C <- x%*%t(x)/(Ntrl -1)
  chC <- cholesky(C)
  
  # entropy in nats
  HX <- sum(log(diag(chC)),na.rm=TRUE) + 0.5*Nvarx*(log(2*pi)+1) # cholesky(x) can give NaNs, numerical problem? absent in python version
  
  ln2 <- log(2)
  
  if (biascorrect) {
    psiterms <- digamma((Ntrl - 1:Nvarx)/2) / 2
    dterm <- (ln2 - log(Ntrl-1)) / 2
    HX <- HX - Nvarx*dterm - sum(psiterms)
  }
  
  #convert to bits
  return(HX/ln2)
}

mi_gg <- function(x, y, biascorrect=TRUE, demeaned=FALSE) {
  # Mutual information (MI) between two Gaussian variables in bits
  
  # I = mi_gg(x,y) returns the MI between two (possibly multidimensional)
  # Gaussian variables, x and y, with bias correction.
  # If x and/or y are multivariate columns must correspond to samples, rows
  # to dimensions/variables. (Samples last axis) 
  
  # biascorrect : true / false option (default true) which specifies whether
  # bias correction should be applied to the esimtated MI.
  # demeaned : false / true option (default false) which specifies whether the
  # input data already has zero mean (true if it has been copula-normalized)
  
  x <- atleast_2d(x)
  y <- atleast_2d(y)
  if (length(dim(x)) > 2 | length(dim(y)) > 2) {
    stop("x and y must be at most 2d")
  }
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  Nvary <- dim(y)[1]
  Nvarxy <- Nvarx + Nvary
  
  if (dim(y)[2] != Ntrl) {
    stop("Number of trials do not match")
  }
  
  # joint variable
  xy <- rbind(x,y)
  if (!demeaned) {
    xy <- sweep(xy,1,rowMeans(xy))
  }
  Cxy <- xy%*%t(xy) / (Ntrl-1)
  # submatrices of joint covariance
  Cx <- Cxy[1:Nvarx,1:Nvarx]
  Cy <- Cxy[(Nvarx+1):Nvarxy,(Nvarx+1):Nvarxy]
  
  chCxy <- cholesky(Cxy)
  chCx <- cholesky(Cx)
  chCy <- cholesky(Cy) 
  
  # entropies in nats
  # normalizations cancel for mutual information
  
  HX <- sum(log(diag(chCx)),na.rm=TRUE) # + 0.5*Nvarx*(log(2*pi) +1)
  HY <- sum(log(diag(chCy)),na.rm=TRUE) # + 0.5*Nvary*(log(2*pi) +1)
  HXY <- sum(log(diag(chCxy)),na.rm=TRUE) # + 0.5*Nvarxy*(log(2*pi) +1)
  
  ln2 = log(2)
  if(biascorrect) {
    psiterms <- digamma((Ntrl - 1:Nvarxy)/2) / 2
    dterm <- (ln2 - log(Ntrl-1)) / 2
    HX <- HX - Nvarx*dterm - sum(psiterms[1:Nvarx])
    HY <- HY - Nvary*dterm - sum(psiterms[1:Nvary])
    HXY <- HXY - Nvarxy*dterm - sum(psiterms[1:Nvarxy])
  }
  
  # MI in bits
  I <- (HX + HY - HXY) / ln2
  return(I)
}

gcmi_cc <- function(x,y){
  # Gaussian-Copula Mutual Information between two continuous variables.
  
  # I = gcmi_cc(x,y) returns the MI between two (possibly multidimensional)
  # continuous variables, x and y, estimated via a Gaussian copula.
  # If x and/or y are multivariate columns must correspond to samples, rows
  # to dimensions/variables. (Samples first axis) 
  # This provides a lower bound to the true MI value.
  
  x <- atleast_2d(x)
  y <- atleast_2d(y)
  if (length(dim(x)) > 2 | length(dim(y)) > 2) {
    stop("x and y must be at most 2d")
  }
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  Nvary <- dim(y)[1]
  
  if (dim(y)[2] != Ntrl) {
    stop("Number of trials do not match")
  }
  
  for (xi in 1:Nvarx){
    if (length(unique(x[xi,])) / Ntrl < 0.9) {
      warning("Input x has more than 10% repeated values")
      break
    }
  }
  for (yi in 1:Nvary){
    if (length(unique(y[yi,])) / Ntrl < 0.9) {
      warning("Input y has more than 10% repeated values")
      break
    }
  }
  
  # copula normalization
  cx <- copnorm(x)
  cy <- copnorm(y)
  
  # parametric Gaussian MI
  I <- mi_gg(cx,cy,TRUE,TRUE)
  return(I)
}

mi_model_gd <- function(x,y,Ym,biascorrect=TRUE, demeaned=FALSE) {
  # Mutual information (MI) between a Gaussian and a discrete variable in bits
  # based on ANOVA style model comparison.
  
  # I = mi_model_gd(x,y,Ym) returns the MI between the (possibly multidimensional)
  # Gaussian variable x and the discrete variable y.
  # For 1D x this is a lower bound to the mutual information.
  # Columns of x correspond to samples, rows to dimensions/variables. 
  # (Samples last axis)
  # y should contain integer values in the range [0 Ym-1] (inclusive).
  
  # biascorrect : true / false option (default true) which specifies whether
  # bias correction should be applied to the esimtated MI.
  # demeaned : false / true option (default false) which specifies whether the
  # input data already has zero mean (true if it has been copula-normalized)
  
  # See also: mi_mixture_gd
  
  x <- atleast_2d(x)
  y <- drop(y)
  
  if (length(dim(x)) > 2) {
    stop("x must be at most 2d")
  }
  if (!is.vector(y)) {
    stop("only univariate discrete variables supported")
  }
  if (!all(is.wholenumber(y))) {
    stop("y should countain only integers (whole number)")
  }
  if (!is.wholenumber(Ym)) {
    stop("Ym should be an integer (whole number)")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  
  if (length(y) != Ntrl) {
    stop("number of trials do not match")
  }
  
  if(!demeaned) {
    x <- sweep(x,1,rowMeans(x))
  }
  
  # class-conditional entropies
  Ntrl_y <- rep(0,Ym)
  Hcond <- rep(0,Ym)
  c <- 0.5*(log(2*pi)+1)
  for (yi in 0:(Ym-1)){
    idx <- y == yi
    xm <- atleast_2d(x[,idx])
    Ntrl_y[yi+1] <- dim(xm)[length(dim(xm))]
    xm <- sweep(xm,1,rowMeans(xm))
    Cm <- xm%*%t(xm) / (Ntrl_y[yi+1] - 1)
    chCm <- cholesky(Cm)
    Hcond[yi+1] <- sum(log(diag(chCm))) # + c*Nvarx
  }
  # class weights
  w <- Ntrl_y / Ntrl
  
  # unconditional entropy from unconditional Gaussian fit
  Cx <- x%*%t(x) / (Ntrl - 1)
  chC <- cholesky(Cx)
  Hunc <- sum(log(diag(chC)),na.rm=TRUE) # + c*Nvarx
  
  ln2 <- log(2)
  
  if (biascorrect){
    vars <- 1:Nvarx
    
    psiterms <- Digamma((Ntrl - vars) / 2) / 2
    dterm <- (ln2 - log(Ntrl - 1)) / 2
    Hunc <- Hunc - Nvarx*dterm - sum(psiterms)
    
    dterm <- (ln2 - log(Ntrl_y - 1)) / 2
    psiterms <- rep(0,Ym)
    for (vi in vars) {
      idx = Ntrl_y - vi
      psiterms <- psiterms + Digamma(idx/2)
    }
    Hcond <- Hcond - Nvarx*dterm - (psiterms/2)
  }
  
  #MI in bits
  I <- (Hunc - sum(w*Hcond)) / ln2
  return(I)
}

gcmi_model_cd <- function(x,y,Ym) {
  # Gaussian-Copula Mutual Information between a continuous and a discrete variable
  # based on ANOVA style model comparison.
  
  # I = gcmi_model_cd(x,y,Ym) returns the MI between the (possibly multidimensional)
  # continuous variable x and the discrete variable y.
  # For 1D x this is a lower bound to the mutual information.
  # Columns of x correspond to samples, rows to dimensions/variables.
  # (Samples last axis)
  # y should contain integer values in the range [0 Ym-1] (inclusive).
  
  # See also: gcmi_mixture_cd
  
  x <- atleast_2d(x)
  y <- drop(y)
  
  if (length(dim(x)) > 2) {
    stop("x must be at most 2d")
  }
  if (!is.vector(y)) {
    stop("only univariate discrete variables supported")
  }
  if (!all(is.wholenumber(y))) {
    stop("y should countain only integers (whole number)")
  }
  if (!is.wholenumber(Ym)) {
    stop("Ym should be an integer (whole number)")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  
  if (length(y) != Ntrl) {
    stop("number of trials do not match")
  }
  
  # check for repeated values
  for (xi in 1:Nvarx){
    if (length(unique(x[xi,])) / Ntrl < 0.9) {
      warning("Input x has more than 10% repeated values")
      break
    }
  }
  
  # check values of discrete variable
  if (min(y) != 0 | max(y) != (Ym-1)) {
    stop("values of discrete variable are out of bounds")
  }
  
  # copula normalization
  cx <- copnorm(x)
  
  # parametric Gaussian MI
  I <- mi_model_gd(cx,y,Ym,TRUE,TRUE)
  return(I)
  
}

mi_mixture_gd <- function(x,y,Ym) {
  # Mutual information (MI) between a Gaussian and a discrete variable in bits
  # calculated from a Gaussian mixture.
  
  # I = mi_mixture_gd(x,y,Ym) returns the MI between the (possibly multidimensional)
  # Gaussian variable x and the discrete variable y.
  # Columns of x correspond to samples, rows to dimensions/variables.
  # (Samples last axis)
  # y should contain integer values in the range [0 Ym-1] (inclusive).
  
  # See also: mi_model_gd
  
  x <- atleast_2d(x)
  y <- drop(y)
  
  if (length(dim(x)) > 2) {
    stop("x must be at most 2d")
  }
  if (!is.vector(y)) {
    stop("only univariate discrete variables supported")
  }
  if (!all(is.wholenumber(y))) {
    stop("y should countain only integers (whole number)")
  }
  if (!is.wholenumber(Ym)) {
    stop("Ym should be an integer (whole number)")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  
  if (length(y) != Ntrl) {
    stop("number of trials do not match")
  }
  
  # class-conditional entropies
  Ntrl_y <- array(0,Ym)
  Hcond <- array(0,Ym)
  m <- array(0,c(Ym,Nvarx))
  w <- array(0,Ym)
  cc <- 0.5*(log(2*pi)+1)
  C <- array(0,c(Nvarx,Nvarx,Ym))
  chC <- array(0,c(Nvarx,Nvarx,Ym))
  for (yi in 0:(Ym-1)){
    # class conditional data
    idx <- y == yi
    xm <- atleast_2d(x[,idx])
    # class mean
    m[yi+1,] <- rowMeans(xm)
    Ntrl_y[yi+1] <- dim(xm)[length(dim(xm))]
    
    xm <- xm - array(m[yi+1,],c(1,Ntrl_y[yi+1]))
    C[,,yi+1] <- xm%*%t(xm) / (Ntrl_y[yi+1] - 1)
    chC[,,yi+1] <- cholesky(C[,,yi+1])
    Hcond[yi+1] <- sum(log(diag(atleast_2d(chC[,,yi+1]))),na.rm=TRUE) + cc*Nvarx
  }
  
  # class weights
  w <- Ntrl_y / Ntrl

  # mixture entropy via unscented transform
  # See:
  # Huber, Bailey, Durrant-Whyte and Hanebeck
  # "On entropy approximation for Gaussian mixture random vectors"
  # http://dx.doi.org/10.1109/MFI.2008.4648062
  
  # Goldberger, Gordon, Greenspan
  # "An efficient image similarity measure based on approximations of 
  # KL-divergence between two Gaussian mixtures"
  # http://dx.doi.org/10.1109/ICCV.2003.1238387
  
  D <- Nvarx
  Ds <- sqrt(Nvarx)
  Hmix = 0
  for (yi in 0:(Ym-1)) {
    Ps <- Ds* t(chC[,,yi+1])
    thsm <- array(m[yi+1,],c(length(m[yi+1]),1))
    # unscented points for this class
    usc <- cbind(thsm + Ps, thsm - Ps)
    
    # class log-likelihoods at unscented points
    log_lik <- array(0,c(Ym,2*Nvarx))
    for (mi in 0:(Ym-1)) {
      # demean points
      dx <- usc - array(m[mi+1,],c(1,length(usc)))
      # Gaussian likelihood
      log_lik[mi+1,] <- norm_innerv(dx,chC[,,mi+1]) - Hcond[mi+1] + 0.5*Nvarx
    }
    
    # log mixture likelihood for these unscented points
    # sum over classes, axis=0 (columns)
    logmixlik <- colLogSumExps(log_lik + log(array(w,c(length(w),2)))) # second term is to incorporate weights w as scaling factor for exp(log_lik): b*exp(a) = exp(a+log(b))
    # add to entropy estimate (sum over unscented points for this class)
    Hmix <- Hmix + w[yi+1]*sum(logmixlik)
  }
  
  Hmix <- - Hmix / (2*D)
  
  # no bias correct
  I <- (Hmix - sum(w*Hcond)) / log(2)
  return(I)
}

norm_innerv <- function(x, chC) {
  # normalised innervations
  m <- solve(chC,x)
  w <- -0.5 * colsums(m * m)
  return(w)
}

gcmi_mixture_cd <- function(x,y,Ym) {
  # Gaussian-Copula Mutual Information between a continuous and a discrete variable
  # calculated from a Gaussian mixture.
  
  # The Gaussian mixture is fit using robust measures of location (median) and scale
  # (median absolute deviation) for each class.
  # I = gcmi_mixture_cd(x,y,Ym) returns the MI between the (possibly multidimensional)
  # continuous variable x and the discrete variable y.
  # For 1D x this is a lower bound to the mutual information.
  # Columns of x correspond to samples, rows to dimensions/variables.
  # (Samples last axis)
  # y should contain integer values in the range [0 Ym-1] (inclusive).
  
  # See also: gcmi_model_cd
  
  x <- atleast_2d(x)
  y <- drop(y)
  
  if (length(dim(x)) > 2) {
    stop("x must be at most 2d")
  }
  if (!is.vector(y)) {
    stop("only univariate discrete variables supported")
  }
  if (!all(is.wholenumber(y))) {
    stop("y should countain only integers (whole number)")
  }
  if (!is.wholenumber(Ym)) {
    stop("Ym should be an integer (whole number)")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  
  if (length(y) != Ntrl) {
    stop("number of trials do not match")
  }
  
  # check for repeated values
  for (xi in 1:Nvarx){
    if (length(unique(x[xi,])) / Ntrl < 0.9) {
      warning("Input x has more than 10% repeated values")
      break
    }
  }
  
  # check values of discrete variable
  if (min(y) != 0 | max(y) != (Ym-1)) {
    stop("values of discrete variable are out of bounds")
  }
  
  # copula normalise each class
  # shift and rescale to match loc and scale of raw data
  # this provides a robust way to fit the Gaussian mixture
  
  classdat <- c()
  ydat <-  c()
  
  for (yi in 0:(Ym-1)) {
    # class conditional data
    idx <- y==yi
    xm <- atleast_2d(x[,idx])
    cxm <- copnorm(xm)
    
    xmmed <- rowMedians(xm)
    # robust measure of s.d. under Gaussian assumption from median absolute deviation
    xmmad <- rowMedians(abs(xm - xmmed))
    cxmscaled <- cxm * (1.482602218505602*xmmad)
    # robust measure of loc from median
    cxmscaled <- cxmscaled + xmmed
    classdat <- c(classdat, cxmscaled)
    ydat <- c(ydat,yi*array(1,dim(xm)[2]))
  }
  
  cx <- array(classdat,c(1,length(classdat)))
  newy <- array(ydat,c(1,length(ydat)))
  I <- mi_mixture_gd(cx,newy,Ym)
  return(I)
  
}

cmi_ggg <- function(x,y,z, biascorrect=TRUE, demeaned=FALSE) {
  # Conditional Mutual information (CMI) between two Gaussian variables
  # conditioned on a third
  
  # I = cmi_ggg(x,y,z) returns the CMI between two (possibly multidimensional)
  # Gaussian variables, x and y, conditioned on a third, z, with bias correction.
  # If x / y / z are multivariate columns must correspond to samples, rows
  # to dimensions/variables. (Samples last axis)
  
  # biascorrect : true / false option (default true) which specifies whether
  # bias correction should be applied to the esimtated MI.
  # demeaned : false / true option (default false) which specifies whether the
  # input data already has zero mean (true if it has been copula-normalized)
  
  x <- atleast_2d(x)
  y <- atleast_2d(y)
  z <- atleast_2d(z)
  
  if (length(dim(x)) > 2 | length(dim(y)) > 2 | length(dim(z)) > 2)  {
    stop("x, y and z must be at most 2d")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  Nvary <- dim(y)[1]
  Nvarz <- dim(z)[1]
  Nvaryz <- Nvary + Nvarz
  Nvarxy <- Nvarx + Nvary
  Nvarxz <- Nvarx + Nvarz
  Nvarxyz <- Nvarx + Nvaryz
  
  if (dim(y)[2] != Ntrl | dim(z)[2] != Ntrl) {
    stop("Number of trials do not match")
  }
  
  # joint variable
  xyz <- rbind(x,y,z)
  if (!demeaned) {
    xyz <- sweep(xyz,1,rowMeans(xyz))
  }
  Cxyz <- xyz%*%t(xyz) / (Ntrl - 1)
  # submatrices of joint covariance
  Cz <- Cxyz[(Nvarxy+1):Nvarxyz,(Nvarxy+1):Nvarxyz]
  Cyz <- Cxyz[(Nvarx+1):Nvarxyz,(Nvarx+1):Nvarxyz]
  Cxz <- array(0,c(Nvarxz,Nvarxz))
  Cxz[1:Nvarx,1:Nvarx] <- Cxyz[1:Nvarx,1:Nvarx]
  Cxz[1:Nvarx,(Nvarx+1):Nvarxz] <- Cxyz[1:Nvarx,(Nvarxy+1):Nvarxyz]
  Cxz[(Nvarx+1):Nvarxz,1:Nvarx] <- Cxyz[(Nvarxy+1):Nvarxyz,1:Nvarx]
  Cxz[(Nvarx+1):Nvarxz,(Nvarx+1):Nvarxz] <- Cxyz[(Nvarxy+1):Nvarxyz,(Nvarxy+1):Nvarxyz]
  
  chCz <- cholesky(Cz)
  chCxz <- cholesky(Cxz)
  chCyz <- cholesky(Cyz)
  chCxyz <- cholesky(Cxyz)
  
  # entropies in nats
  # normalizations cancel for cmi
  HZ <- sum(log(diag(chCz)),na.rm=TRUE) # + 0.5*Nvarz*(log(2*pi)+1.0)
  HXZ <- sum(log(diag(chCxz)),na.rm=TRUE) # + 0.5*Nvarxz*(log(2*pi)+1.0)
  HYZ <- sum(log(diag(chCyz)),na.rm=TRUE) # + 0.5*Nvaryz*(log(2*pi)+1.0)
  HXYZ <- sum(log(diag(chCxyz)),na.rm=TRUE) # + 0.5*Nvarxyz*(log(2*pi)+1.0)
  
  ln2 <- log(2)
  if (biascorrect){
    psiterms <- Digamma((Ntrl - 1:Nvarxyz) / 2) / 2
    dterm <- (ln2 - log(Ntrl-1)) / 2
    HZ <- HZ - Nvarz*dterm - sum(psiterms[1:Nvarz])
    HXZ <- HXZ - Nvarxz*dterm - sum(psiterms[1:Nvarxz])
    HYZ <- HYZ - Nvaryz*dterm - sum(psiterms[1:Nvaryz])
    HXYZ <- HXYZ - Nvarxyz*dterm - sum(psiterms[1:Nvarxyz])
  }
  
  # MI in bits
  I <- (HXZ + HYZ - HXYZ - HZ) / ln2
  return(I)
  
}

gccmi_ccc <- function(x,y,z){
  # Gaussian-Copula CMI between three continuous variables.
  
  # I = gccmi_ccc(x,y,z) returns the CMI between two (possibly multidimensional)
  # continuous variables, x and y, conditioned on a third, z, estimated via a
  # Gaussian copula.
  # If x and/or y are multivariate columns must correspond to samples, rows
  # to dimensions/variables. (Samples first axis)
  
  x <- atleast_2d(x)
  y <- atleast_2d(y)
  z <- atleast_2d(z)
  
  if (length(dim(x)) > 2 | length(dim(y)) > 2 | length(dim(z)) > 2)  {
    stop("x, y and z must be at most 2d")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  Nvary <- dim(y)[1]
  Nvarz <- dim(z)[1]
  
  if (dim(y)[2] != Ntrl | dim(z)[2] != Ntrl) {
    stop("Number of trials do not match")
  }
  
  # check for repeated values
  for (xi in 1:Nvarx){
    if (length(unique(x[xi,])) / Ntrl < 0.9) {
      warning("Input x has more than 10% repeated values")
      break
    }
  }
  for (yi in 1:Nvary){
    if (length(unique(y[yi,])) / Ntrl < 0.9) {
      warning("Input y has more than 10% repeated values")
      break
    }
  }
  for (zi in 1:Nvarz){
    if (length(unique(z[zi,])) / Ntrl < 0.9) {
      warning("Input z has more than 10% repeated values")
      break
    }
  }
  
  # copula normalization
  cx <- copnorm(x)
  cy <- copnorm(y)
  cz <- copnorm(z)
  # parametric Gaussian CMI
  I <- cmi_ggg(cx,cy,cz,TRUE,TRUE)
  return(I)
}

gccmi_ccd <- function(x,y,z,Zm){
  # Gaussian-Copula CMI between 2 continuous variables conditioned on a discrete variable.
  
  # I = gccmi_ccd(x,y,z,Zm) returns the CMI between two (possibly multidimensional)
  # continuous variables, x and y, conditioned on a third discrete variable z, estimated
  # via a Gaussian copula.
  # If x and/or y are multivariate columns must correspond to samples, rows
  # to dimensions/variables. (Samples first axis)
  # z should contain integer values in the range [0 Zm-1] (inclusive).
  
  x <- atleast_2d(x)
  y <- atleast_2d(y)
  
  if (length(dim(x)) > 2 | length(dim(y)) > 2) {
    stop("x and y must be at most 2d")
  }
  if (!is.vector(z)) {
    stop("only univariate discrete variables supported")
  }
  if (!all(is.wholenumber(z))) {
    stop("z should countain only integers (whole number)")
  }
  if (!is.wholenumber(Zm)) {
    stop("Zm should be an integer (whole number)")
  }
  
  Ntrl <- dim(x)[2]
  Nvarx <- dim(x)[1]
  Nvary <- dim(y)[1]
  
  if (dim(y)[2] != Ntrl | length(z) != Ntrl) {
    stop("Number of trials do not match")
  }
  
  # check for repeated values
  for (xi in 1:Nvarx){
    if (length(unique(x[xi,])) / Ntrl < 0.9) {
      warning("Input x has more than 10% repeated values")
      break
    }
  }
  for (yi in 1:Nvary){
    if (length(unique(y[yi,])) / Ntrl < 0.9) {
      warning("Input y has more than 10% repeated values")
      break
    }
  }
  
  # check values of discrete variable
  if (min(z) != 0 | max(z) != (Zm-1)) {
    stop("values of discrete variable are out of bounds")
  }
  
  # calculate gcmi for each z value
  Icond <- array(0,Zm)
  Pz <- array(0,Zm)
  cx <- c()
  cy <- c()
  for (zi in 0:(Zm-1)){
    idx <- z==zi
    thsx <- copnorm(x[,idx])
    thsy <- copnorm(y[,idx])
    Pz[zi+1] <- sum(idx)
    cx <- c(cx,thsx)
    cy <- c(cy,thsy)
    Icond[zi+1] <- mi_gg(thsx,thsy,TRUE,TRUE)
  }
  
  Pz <- Pz / Ntrl
  # conditional mutual information
  CMI <- sum(Pz*Icond)
  I <- mi_gg(cx,cy,TRUE, FALSE)
  return(array(c(CMI,I),2,dimnames = list(c('CMI','I'))))
}

