##############################################################################################
############ Helper functions to calculate information theoretical measures for   ############
############ A hitchhiker's guide to information-theoretic measures in psychology ############
############               N. Van santen, Y. Roseel, D. Marinazzo                 ############
##############################################################################################


library(progress)

discretize_equal_width <- function(values, bins) {
  discretized <- cut(values, breaks = bins, include.lowest = TRUE, labels = FALSE)
  return(discretized)
}

discretize_quantile <- function(values, bins) {
  # Calculate quantile probabilities
  probs <- seq(0, 1, length.out = bins + 1)
  
  # Get quantile breaks
  breaks <- quantile(values, probs = probs, na.rm = TRUE)
  
  # Handle duplicates by using unique quantiles
  breaks <- unique(breaks)
  
  # If we have fewer unique breaks than requested bins + 1,
  # we effectively have fewer bins
  if (length(breaks) < bins + 1) {
    warning("Number of unique quantiles is less than requested bins. Using ", 
            length(breaks) - 1, " bins instead of ", bins, ".")
  }
  
  # Discretize using the quantile breaks
  discretized <- cut(values, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  return(discretized)
}

discretize <- function(values, bins, method="equal_width"){
  
  switch(method,
         equal_width = {return(discretize_equal_width(values,bins))},
         quantile = {return(discretize_quantile(values,bins))})
}

get.lambda.shrink <- function(n, u, t, verbose) { #from entropy package (Hausser, Strimmer)
  
  # *unbiased* estimator of variance of u
  varu <- u*(1-u)/(n-1)
  
  # misspecification
  msp <- sum( (u-t)^2 )
  
  # estimate shrinkage intensity  
  if (msp == 0)
  {
    #warning("Overshrinkage")
    lambda <- 1
  }
  else
    lambda <- sum( varu ) / msp
  
  if (lambda > 1)
  {
    lambda <- 1 # truncate at 1
    #warning("Overshrinkage")
  }
  
  if (lambda < 0)
  {
    lambda <- 0
    #warning("Undershrinkage")
  }
  
  if (verbose)
  {
    cat(paste("Estimating optimal shrinkage intensity lambda.freq (frequencies):", 
              round(lambda, 4)) , "\n")
  }
  
  return(lambda)
}

freqs.shrink <- function (y, lambda.freqs, verbose = FALSE) { # From entropy package (Hausser, Strimmer)
  
  target <- 1/length(y) # uniform target (note length() works also for matrices)
  n <- sum(y)
  u <- y/n
  
  if (missing(lambda.freqs))
  {
    if (n==1 || n==0)
    {
      lambda.freqs <- 1
    }
    else
    {
      lambda.freqs <- get.lambda.shrink(n, u, target, verbose)
    }
    
  }
  else
  {
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity lambda.freq (frequencies):", 
                round(lambda.freqs, 4)) , "\n")
    }
    
  }
  u.shrink <- lambda.freqs * target + (1 - lambda.freqs) * u
  
  # attr(u.shrink, "lambda.freqs") <- lambda.freqs
  
  return(u.shrink)
}

entropy_hist <- function(data , X, discrete, binning_method, numbins, units, biascorrection){
  # data: dataframe with variables as columns and observations as rows
  # X:  names of variables for which to calculate entropy, must be character type.
  #     NULL leads to joint entropy over entire dataset
  # discrete: Boolean to denote if the variables are discrete or continuous
  # binning:  character type to denote the binning strategy. choices are "L2", "cuberoot"
  # units:  character type to denote units of entropy. Choices are "bits" or "nats"
  
  X <- unique(X)
  
  
  if (is.null(X)){X <- colnames(data)}
  
  if (length(discrete) != 1 &&  length(discrete) != length(X)){
    stop("The parameter 'discrete' needs to be a single boolan value or a vector of boolean values corresponding to the variables in X")
  }
  
  N <- dim(data)[1]
  M <- length(X)
  
  
  if (length(discrete) == 1) {
    discrete <- rep(discrete, M)
  }
  
  
  data_discretized <- data.frame(matrix(nrow=N,ncol= M))
  colnames(data_discretized) <- X
  
  bin_widths <- numeric(M)
  
  for(i in 1:M){
    
    if(discrete[i]){
      
      data_discretized[,i] <- data[,X[i]]
      
    } else {
      
      if(binning_method == "equal_width"){
        
        if(is.double(numbins)){
          data_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=numbins)
          bin_widths[i] <- diff(range(data[,X[i]]))/length(table(data_discretized[,i]))
        } else{
          
          switch(numbins,
                 L2 = {
                   if(IQR(data[,X[i]])==0){
                     data_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=(N^(1/3)))
                     bin_widths[i] <- diff(range(data[,X[i]]))/length(table(data_discretized[,i]))
                   } else {
                     bin_widths[i] <- 2*IQR(data[,X[i]])/N^(1/3)
                     number_of_bins <- round(diff(range(data[,X[i]]))/bin_widths[i])
                     data_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=number_of_bins)
                   }
                 },
                 cuberoot = {
                   data_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=(N^(1/3)))
                   bin_widths[i] <- diff(range(data[,X[i]]))/length(table(data_discretized[,i]))
                 })
        }
      } else if (binning_method =="quantile") {
        
        data_discretized[,i] <- discretize(data[,X[i]],bins = numbins, method="quantile")
        bin_widths[i] <- diff(range(data[,X[i]]))/length(table(data_discretized[,i]))
        
      } 
      
    }
    
  }
  
  #stopifnot(is.data.frame(data_discretized))
  #if (any(sapply(data_discretized, function(col) all(is.na(col))))) {
  #  stop("One or more columns in 'data_discretized' are completely NA. Cannot compute state counts.")
  #}
  #data_discretized[] <- lapply(data_discretized, as.factor)
  
  state_counts <- table(data_discretized)
  
  
  if(all(discrete)){
    
    bin_correction <- 1
    
  } else {
    
    bin_correction <- prod(bin_widths[!discrete])
  }
  
  #probabilities
  switch(biascorrection,
         ml = {
           p_states <- state_counts / N
           K <- 1
         },
         mm = {
           p_states <- state_counts / N
           K <- sum(table(data_discretized) != 0)
         },
         shrink = {
           p_states <- freqs.shrink(state_counts)
           K <- 1
         })
  
  entropy <- -sum(p_states * log2(p_states/bin_correction), na.rm = TRUE) + (K - 1)/(2*N*log(2))
  
  
  if(entropy < 0){
    warning("Entropy is negative")
    #entropy <- 0
  }
  
  switch(units,
         bits = {return(entropy)},
         nats = {return(entropy*log(2))})
  
  
}

entropy_cov <- function(data, X, Y, units, biascorrection){
  
  X <- unique(X)
  Y <- unique(Y)
  
  if (is.null(X)){X <- colnames(data)}
  M <- length(X)
  N <- dim(data)[1]
  
  if(is.null(Y)){
    
    cov_data <- cov(data[X])
    
    if(det(cov_data) == 0){ #singular
      entropy <- entropy_cov(data = data, X = X[1], Y = NULL, units = units, biascorrection = biascorrection)
    } else {
      entropy <- 0.5*M*log(2*pi*exp(1)) + 0.5*log(det(cov_data))
    }
    
    if(biascorrection == 'Gauss'){
      
      psis <- vector(mode = "numeric",length = M)
      for(i in 1:M){
        psis[i] <- digamma((N-i)/2)
      }
      entropy <- entropy - 0.5*(M*log(2/(N-1)) + sum(psis))
    }
    
  } else {
    cov_data <- cov(data)
    P <- length(Y)
    
    cov_XX <- cov_data[X, X, drop = FALSE]
    cov_XY <- cov_data[X, Y, drop = FALSE]
    cov_YY <- cov_data[Y, Y, drop = FALSE]
    
    res.cov <- cov_XX - cov_XY %*% solve(cov_YY) %*% t(cov_XY)
    
    entropy <- 0.5*M*log(2*pi*exp(1)) + 0.5*log(det(res.cov))
    
    
    if(biascorrection == 'Gauss'){
      
      psis <- vector(mode = "numeric",length = P)
      for(i in 1:P){
        psis[i] <- digamma((N-(M+i))/2)
      }
      entropy <- entropy - 0.5*(P*log(2/(N-1)) + sum(psis))
    }
  }
  
  if(entropy < 0){
    warning("Entropy is negative")
    #entropy <- 0
  }
  
  switch(units,
         nats = {return(entropy)},
         bits = {return(entropy*log2(exp(1)))})
}

Entropy <- function(data, X=NULL, Y=NULL, estimator="hist", discrete = T, biascorrection = "ml", units="nats", binning_method = "equal_width", numbins = "L2") {
  
  # data: vector, n x m matrix or dataframe with n samples of m variables.
  # X: name (char) or vector of names of variables to calculate entropy of. If NULL, joint entropy of data is calculated
  # Y: name or names (char) of variables to condition on
  # estimator: method of estimation. "hist" or "cov".
  # discrete: Boolean parameter to signify if data is discrete (T) or not (F). Only relevant when estimator is "hist".
  # biascorrection: method of biascorrection. "ml" (maximum likelyhood or without correction), "mm" (miller madow correction), or "shrink" (James-stein shrinkage estimation)
  
  if (is.vector(data)) {
    data <- data.frame(data)
  }
  
  if (!all(X %in% colnames(data))) {
    stop("Some variables in X are not present in the data.")
  }
  
  if (!all(Y %in% colnames(data))) {
    stop("Some variables in Y are not present in the data.")
  }
  
  if (!numbins %in% c("L2", "cuberoot") && !is.double(numbins)) {
    stop("Invalid binning strategy. Choose 'L2' or 'cuberoot' with binning method 'equal_width' or a number with binning method 'quantile'.")
  }
  
  if (!units %in% c("bits", "nats")) {
    stop("Invalid units. Choose 'bits' or 'nats'.")
  }
  
  if (!estimator %in% c("hist","cov")) {
    stop("Invalid estimator. Choose 'hist' or 'cov'.")
  }
  
  if (!biascorrection %in% c("ml","mm","shrink", "Gauss")) {
    stop("Invalid estimator. Choose 'ml', 'mm', 'shrink' or 'Gauss'.")
  }
  
  if(!binning_method %in% c("equal_width","quantile")) {
    stop("parameter binnin_method needs to be 'equal_width' or 'quantile'.")
  }
  
  switch(estimator,
         hist = {
           if(is.null(Y)){
             return(entropy_hist(data,X,discrete = discrete,biascorrection=biascorrection,units=units,binning_method=binning_method, numbins=numbins))
           } else {
             return(entropy_hist(data,c(X,Y),discrete = discrete,biascorrection=biascorrection,units=units,binning_method=binning_method, numbins=numbins) 
                    - entropy_hist(data,Y,discrete = discrete,biascorrection=biascorrection,units=units,binning_method=binning_method, numbins=numbins))
           }
         },
         cov = {
           return(entropy_cov(data,X,Y,units,biascorrection))
         }
  )
}

mutual_info_hist <- function(data , X, Y, discrete, binning, units, biascorrection){
  
  X <- unique(X)
  Y <- unique(Y)
  
  N <- dim(data)[1]
  M_X <- length(X)
  M_Y <- length(Y)
  
  if (length(discrete) == 1) {
    discrete_X <- rep(discrete, M_X)
    discrete_Y <- rep(discrete, M_Y)
  } else if (length(discrete == 2)){
    discrete_X <- rep(discrete[1], M_X)
    discrete_Y <- rep(discrete[2], M_Y)
  } else {
    discrete_X <- discrete[1:M_X]
    discrete_Y <- discrete[M_X+1:M_X+M_Y]
  }
  
  
  X_discretized <- data.frame(matrix(nrow=N,ncol= M_X))
  colnames(X_discretized) <- X
  
  bin_widths_X <- numeric(M_X)
  
  Y_discretized <- data.frame(matrix(nrow=N,ncol= M_Y))
  colnames(Y_discretized) <- Y
  
  bin_widths_Y <- numeric(M_Y)
  
  for(i in 1:M_X) {
    
    if(discrete_X[i]){
      
      X_discretized[,i] <- data[,X[i]]
      
    } else {
      
      switch(binning,
             L2 = {
               if(IQR(data[,X[i]])==0){
                 X_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=(N^(1/3)))
                 bin_widths_X[i] <- diff(range(data[,X[i]]))/length(table(X_discretized[,i]))
               } else {
                 bin_widths_X[i] <- 2*IQR(data[,X[i]])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,X[i]]))/bin_widths_X[i])
                 X_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=number_of_bins)
               }
             },
             cuberoot = {
               X_discretized[,i] <- discretize_equal_width(data[,X[i]],bins=(N^(1/3)))
               bin_widths_X[i] <- diff(range(data[,X[i]]))/length(table(X_discretized[,i]))
             })
      
    }
    
  }
  
  for(i in 1:M_Y) {
    
    if(discrete_Y[i]){
      
      Y_discretized[,i] <- data[,Y[i]]
      
    } else {
      
      switch(binning,
             L2 = {
               if(IQR(data[,Y[i]])==0){
                 Y_discretized[,i] <- discretize_equal_width(data[,Y[i]],bins=(N^(1/3)))
                 bin_widths_Y[i] <- diff(range(data[,Y[i]]))/length(table(Y_discretized[,i]))
               } else {
                 bin_widths_Y[i] <- 2*IQR(data[,Y[i]])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,Y[i]]))/bin_widths_Y[i])
                 Y_discretized[,i] <- discretize_equal_width(data[,Y[i]],bins=number_of_bins)
               }
             },
             cuberoot = {
               Y_discretized[,i] <- discretize_equal_width(data[,Y[i]],bins=(N^(1/3)))
               bin_widths_Y[i] <- diff(range(data[,Y[i]]))/length(table(Y_discretized[,i]))
             })
      
    }
    
  }
  
  data_discretized <- cbind(X_discretized,Y_discretized)
  
  state_counts_X <- table(data_discretized[X])
  state_counts_Y <- table(data_discretized[Y])
  state_counts_joint <- table(data_discretized[c(X,Y)])
  
  
  switch(biascorrection,
         ml = {
           p_states_X <- state_counts_X / N
           p_states_Y <- state_counts_Y / N
           p_states_joint <- state_counts_joint / N
           
           correction <- 0
         },
         mm = {
           p_states_X <- state_counts_X / N
           p_states_Y <- state_counts_Y / N
           p_states_joint <- state_counts_joint / N
           
           K_X <- sum(table(data_discretized[X]) != 0)
           K_Y <- sum(table(data_discretized[Y]) != 0)
           K_XY <- sum(table(data_discretized[c(X,Y)]) != 0)
           
           correction <- ((K_XY - 1)-(K_X - 1)-(K_Y - 1))/(2*N*log(2))
         },
         shrink = {
           p_states_X <- freqs.shrink(state_counts_X)
           p_states_Y <- freqs.shrink(state_counts_Y)
           p_states_joint <- freqs.shrink(state_counts_joint)
           
           correction <- 0
         })
  
  
  mi <- sum(p_states_joint * log2(p_states_joint / (p_states_X %o% p_states_Y)), na.rm=TRUE) - correction
  
  
  switch(units,
         bits = {return(mi)},
         nats = {return(mi*log(2))})
  
}

mutual_info_cov <- function (data, X, Y, Z, units, biascorrection) {
  
  if(biascorrection == 'Gauss'){
    mi <- entropy(data, X, Z, biascorrection) - entropy_cov(data, X, c(Y,Z), biascorrection)
  } else {
    
    X <- unique(X)
    Y <- unique(Y)
    Z <- unique(Z)
    
    if(is.null(Z)){
      
      #cor_data <- cor(data)
      #cor_XY <- cor_data[X,Y]
      #mi <- -0.5*log(1 - (cor_XY*cor_XY))
      
      # Get the covariance matrix for the full set {X, Y}
      cov_data <- cov(data[, c(X, Y)])
      
      # Compute the determinant of covariance matrices
      det_X <- det(cov(data[, X, drop = FALSE]))  # |Σ_X|
      det_Y <- det(cov(data[, Y, drop = FALSE]))  # |Σ_Y|
      det_XY <- det(cov_data)                     # |Σ_XY|
      
      # Ensure determinants are positive
      if (det_X > 0 && det_Y > 0 && det_XY > 0) {
        mi <- 0.5 * log(det_Y * det_X / det_XY)  # Multivariate MI formula
      } else {
        mi <- 0  # If covariance is ill-conditioned
      }
      
    } else {
      
      cov_data <- cov(data[,c(X,Y,Z)])
      vars <- c(X,Y)
      res.cov <- (cov_data[vars,vars] - cov_data[vars,Z]%*% solve(cov_data[Z,Z])%*%t(cov_data[vars,Z]))
      if(all(diag(res.cov) > 0)) {
        res.cor <- cov2cor(res.cov)[1,2]
        if(abs(res.cor) < 0.999) {
          mi <- -1/2 * log(1 - (res.cor*res.cor))
        } else {print("Mutual information is infinite")}
      }
      else { mi <- 0}
    }
  }
  switch(units,
         nats = {return(mi)},
         bits = {return(mi*log2(exp(1)))})
  
}

mutual_info_gc <- function (data, X, Y, Z, units){
  
  X <- unique(X)
  Y <- unique(Y)
  Z <- unique(Z)
  
  if(is.null(Z)){
    mi <- gcmi_cc(t(data[X]),t(data[Y]))
  } else {
    mi <- gccmi_ccc(t(data[X]),t(data[Y]),t(data[Z]))
  }
  
  switch(units,
         nats = {return(mi*log(2))},
         bits = {return(mi)})
  
}

Mutual_Info <- function(data, X=NULL, Y=NULL, Z=NULL, estimator="hist", discrete = T, biascorrection = "ml", units="nats", binning = "L2"){
  
  
  if(!is.data.frame(data)) {
    warning("data input is not a data frame.")
    data <- data.frame(data)
  }
  
  if(is.null(X) || is.null(Y)){stop("X or Y variables not specified")}
  
  if (!all(X %in% colnames(data))) {
    stop("Some variables in X are not present in the data.")
  }
  
  if (!all(Y %in% colnames(data))) {
    stop("Some variables in Y are not present in the data.")
  }
  
  if (!binning %in% c("L2", "cuberoot")) {
    stop("Invalid binning strategy. Choose 'L2' or 'cuberoot'.")
  }
  
  if (!units %in% c("bits", "nats")) {
    stop("Invalid units. Choose 'bits' or 'nats'.")
  }
  
  if (!estimator %in% c("hist","cov","gc")) {
    stop("Invalid estimator. Choose 'hist', 'cov' or 'gc'.")
  }
  
  if (!biascorrection %in% c("ml","mm","shrink")) {
    stop("Invalid estimator. Choose 'ml', 'mm' or 'shrink'.")
  }
  
  if (length(discrete) != 1 &&  length(discrete) != 2 && length(discrete) != (length(X) + length(Y))){
    stop("The parameter 'discrete' needs to be a single boolan value or a vector of boolean values of length 2 (correspending to a value for the X variables and the Y variables) or a vector with boolean values corresponding to the combined number of variables in X and Y")
  }
  
  switch(estimator,
         hist = {
           if(is.null(Z)){
             return(mutual_info_hist(data, X = X, Y = Y, discrete = discrete, binning = binning, units = units, biascorrection = biascorrection))
           } else {
             return(mutual_info_hist(data, X = X, Y = c(Y,Z), discrete = discrete, binning = binning, units = units, biascorrection = biascorrection)
                    - mutual_info_hist(data, X = X, Y = Z, discrete = discrete, binning = binning, units = units, biascorrection = biascorrection))
           }
         },
         cov = {return(mutual_info_cov(data = data, X = X, Y = Y, Z = Z, units = units, biascorrection = biascorrection))},
         gc = {return(mutual_info_gc(data = data, X = X, Y = Y, Z = Z, units = units))})
  
  
  
}

Co_Info <- function(data, X=NULL, Y=NULL, Z=NULL, estimator="hist", discrete=TRUE, biascorrection = "ml", binning = "L2", units = "nats") {
  
  if(!is.data.frame(data)){data <- data.frame(data)}
  
  if(is.null(X) || is.null(Y) || is.null(Z)){stop("X, Y or Z variables not specified")}
  
  mi <- Mutual_Info(data = data, X = X, Y = Y, estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  mi_cond <- Mutual_Info(data = data, X = X, Y = Y, Z = Z, estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  
  return(mi - mi_cond)
}

Total_Corr <- function(data, X=NULL, estimator="hist", binning_method="equal_width", numbins = "L2", units = "nats", discrete=T, biascorrection = "ml") {
  
  if(is.null(X)){ #if the variables are not specified than the TC over all variables will be returned
    
    X <- colnames(data)
    
  }
  
  entropies <- c()
  
  for (variable in X) {
    entropies <- c(entropies, Entropy(data,X=variable,estimator = estimator, discrete=discrete, biascorrection = biascorrection, binning_method = binning_method, numbins=numbins, units=units))
  }
  joint_entropy <- Entropy(data,X=X,estimator=estimator, discrete=discrete, biascorrection=biascorrection, binning_method=binning_method, numbins=numbins, units=units)
  
  TC <- sum(entropies) - joint_entropy
  
  return(TC)
  
}

Resid_Entropy <- function(data, X, estimator="hist", binning_method="equal_width", numbins = "L2", units = "nats",discrete=discrete, biascorrection="ml") {
  
  environment_X <- colnames(data)[colnames(data) != X]
  
  R_X <- Entropy(data, X=X, Y=environment_X,estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning_method = binning_method, numbins=numbins)
  
  return(R_X)
}

Dual_Total_Corr <- function(data, X=NULL, estimator="hist", binning_method="equal_width", numbins = "L2", units = "nats", discrete=T, biascorrection="ml") {
  
  if(is.null(X)){ #if the variables are not specified than the TC over all variables will be returned
    
    X <- colnames(data)
    
  }
  
  joint_entropy <- Entropy(data,X=X,estimator=estimator, binning_method=binning_method, numbins=numbins, units=units,discrete=discrete,biascorrection = biascorrection)
  
  resid_entropies <- c()
  for (variable in X) {
    resid_entropies <- c(resid_entropies, Resid_Entropy(data[X],X=variable, estimator=estimator, binning_method=binning_method, numbins=numbins, units=units,discrete=discrete,biascorrection=biascorrection))
  }
  
  DTC <- joint_entropy - sum(resid_entropies)
  
  return(DTC)
  
}

O_Info <- function(data, X=NULL, estimator = "hist", binning_method="equal_width", numbins = "L2", units = "nats", discrete=T, biascorrection = "ml") {
  
  if(estimator=='gc'){
    if(is.null(X)){ #if the variables are not specified than the TC over all variables will be returned
      
      X <- colnames(data)
      
    }
    
    terms <- vector(mode = "numeric",length = length(2:(length(X)-1)))
    
    for (i in 2:(length(X)-1)){
      terms[i-1] <- Co_Info(data, X=X[i], Y=X[1:(i-1)], Z=X[(i+1):length(X)],estimator = 'gc',units = units) 
    }
    
    O <- sum(terms)
  } else{
    TC <- Total_Corr(data, X=X, estimator=estimator, binning_method=binning_method, numbins=numbins, units=units, discrete=discrete, biascorrection=biascorrection)
    DTC <- Dual_Total_Corr(data, X=X, estimator=estimator, binning_method=binning_method, numbins=numbins, units=units, discrete=discrete, biascorrection=biascorrection)
    
    O <- TC - DTC
  }
  
  return(O)
}

surr_shuf <- function(x, typeshuf='same') {
  
  if (is.null(dim(x))){
    N <- length(x)
    M <- 1
  }
  else {
    N <- dim(x)[1]
    M <- dim(x)[2]
  }
  
  if (M==1){
    typeshuf = 'same'
  }
  
  xs <- matrix(nrow=N,ncol=M)
  p <- xs
  
  switch(typeshuf,
         same = {
           p <- sample(N)
           if(M==1){
             xs <- x[p]
           }
           else {
             xs <- x[p,] 
           }
         },
         ind = {
           for (m in 1:M){
             p[,m] <- sample(N)
             xs[,m] <- x[p[,m],m]
           } 
         })
  
  xs
}

O_surr <- function(x, nsurr, estimator = "hist", biascorrection="ml", discrete = T, units= "bits") {
  # x should be the data of the relevant variables only
  
  Os <- numeric(nsurr)
  
  nvars <- dim(x)[2]
  
  for (i in 1:nsurr) {
    
    j <- sample.int(nvars,1)
    xx <- x
    
    xx[,j] <- surr_shuf(x[,j])
    Os[i] <- O_Info(xx, estimator = estimator, discrete = discrete , units = units, biascorrection = biascorrection)
    
  }
  
  return(Os)
}

surrSignif <- function(data, vars=NULL, order=3, nsurr=10, alpha=0.05, mtcorrection = "BH", estimator="hist", biascorrection="mm", units="bits", discrete=T, visual = T, showsurrs=T, ordered = F, verbose = T){
  
  N <- nrow(data)  # Number of observations
  
  if(is.null(vars)){
    tuples <- combn(colnames(data),order) # create tuples
  } else {
    tuples <- combn(vars,order) # create tuples
  }
  ncomb <- ncol(tuples)  # Number of tuples
  
  # Store original and surrogate results
  original_O <- numeric(ncomb)
  surr_O <- matrix(0, nrow = ncomb, ncol = nsurr)
  significant <- vector(mode = "logical", length = ncomb)
  p_values <- vector(mode ="numeric", length=ncomb)
  CI <- matrix(,nrow=ncomb,ncol=2)
  overlap <- vector(mode = "logical", length = ncomb)
  
  if(verbose){  
    pb <- progress_bar$new(
      format = "  calculating [:bar] :percent in :elapsed eta: :eta",
      total = ncomb, clear = FALSE, width= 60)
  }
  
  # Helper: CI interval intersection
  intervals_intersect <- function(a, b) max(a[1], b[1]) <= min(a[2], b[2])
  
  for (i in 1:ncomb){
    
    original_O[i] <- O_Info(data,tuples[,i],estimator = estimator,biascorrection=biascorrection, units=units, discrete = discrete)
    surr_O[i,] <- O_surr(data[,tuples[,i]],estimator = estimator, biascorrection = biascorrection, nsurr = nsurr, discrete = discrete, units = units)
    
    #range_surrs <- range(surr_O[i,])
    #ifelse(original_O[i] > range_surrs[1] && original_O[i] < range_surrs[2],significant[i]<-F,significant[i]<-T)
    
    if(original_O[i] > mean(surr_O[i,])) {
      p_values[i] <- (sum(surr_O[i,]>original_O[i]) + 1) / (nsurr + 1)
    } else if (original_O[i] < mean(surr_O[i,])) {
      p_values[i] <- (sum(surr_O[i,]<original_O[i]) + 1) / (nsurr + 1)
    } else {p_values[i] <- 1}
    
    CI[i, ] <- quantile(surr_O[i, ], probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    
    #if(order>3){
    
    #  sub_tuples <- combn(tuples[,i],(order-1),simplify = F)
    #  found_overlap <- FALSE
    
    #  for(S in sub_tuples){
    #    surr_S <- O_surr(
    #      data[, S, drop = FALSE],
    #      estimator = estimator, biascorrection = biascorrection,
    #      nsurr = nsurr, discrete = discrete, units = units)
    
    #    CI_sub <- quantile(surr_S, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    
    
    #    if (intervals_intersect(CI[i,], CI_sub)) {
    #      found_overlap <- TRUE
    #      break
    #    }
    
    #  }
    #  overlap[i] <- found_overlap
    
    #} else{overlap[i] <- F}
    
    if(verbose){
      pb$tick()
    }
    
    
  }
  
  
  p_adjusted <- p.adjust(p_values, method = mtcorrection )
  
  
  significant <- p_adjusted <= alpha #& !overlap
  
  # Initialize empty tables with all variable names
  var_names <- unique(as.vector(tuples))
  redundant_counts <- table(factor(rep(var_names, 0), levels = var_names))
  synergistic_counts <- table(factor(rep(var_names, 0), levels = var_names))
  
  # Separate significant tuples by O-Information sign
  if (any(significant)) {
    sig_idxs <- which(significant)
    
    for (i in sig_idxs) {
      vars <- tuples[, i]
      if (original_O[i] > mean(surr_O[i,])) {
        redundant_counts <- redundant_counts + table(factor(vars, levels = var_names))
      } else if (original_O[i] < mean(surr_O[i,])) {
        synergistic_counts <- synergistic_counts + table(factor(vars, levels = var_names))
      }
      # If original_O[i] == 0, we skip (neutral, neither redundant nor synergistic)
    }
  }
  
  type <- ifelse(original_O > rowMeans(surr_O),"redundant","synergistic")
  
  if(visual){
    
    if(ordered){
      order_index <- order(original_O)
    } else {order_index <- 1:length(original_O)}
    
    x_vals <- 1:ncomb
    
    mean_surr <- apply(surr_O,1,mean)
    sd_surr <- apply(surr_O,1,sd)
    
    meanplussd <- mean_surr + sd_surr
    meanminussd <- mean_surr - sd_surr
    
    y_min <- min(c(original_O, surr_O, meanminussd))
    y_max <- max(c(original_O, surr_O, meanplussd))
    
    plot(x_vals, original_O[order_index], col = ifelse(significant,"green","black")[order_index], pch = 16,
         xlab = "Tuple Index", ylab = "O-Information", main = "Original vs Surrogate O-Information", ylim = c(y_min, y_max))
    
    # Add red vertical lines showing the min-max range of surrogates
    #for (i in 1:ncomb) {
    #  idx <- order_index[i]
    #  lines(rep(x_vals[i], 2), range(surr_O[idx, ]), col = "red", lwd = 2)
    #}
    
    if(showsurrs){
      for (i in 1:nsurr){
        points(x_vals,surr_O[,i][order_index],col=rgb(1,0,0,alpha=0.25))
      }
    }
    
    #points(x_vals,mean_surr[order_index],col="red",type='l')
    #polygon(c(x_vals,rev(x_vals)),c(meanplussd[order_index],meanminussd[order_index]),col=rgb(1,0,0,alpha=0.1),border = NA)
    
    
    grid()
    
    if(showsurrs){  
      legend("topleft", legend = c("Original significant", "Original insignificant", "Surrogates"),
             col = c("green", "black", "red"), pch = c(16, 16, 1), lty = c(NA, NA, NA), lwd = c(NA,NA, NA))
    } else {
      legend("topleft", legend = c("Original significant", "Original insignificant"),
             col = c("green", "black"), pch = c(16, 16), lty = c(NA, NA), lwd = c(NA,NA))
    }
    
    
  }
  
  if(verbose){print(paste0("Number of significant tuples of order ", order, ": ",sum(significant),"/",ncomb))}
  
  return(list(tuples = tuples, 
              original_O = original_O, 
              surr_O = surr_O, 
              CI = CI,
              p_values = p_values,
              p_adjusted = p_adjusted,
              significant = significant,
              type = type,
              redundant_counts = redundant_counts,
              synergistic_counts = synergistic_counts))
  
  
}

MMI <- function(data, S, R, composite= NULL, estimator = "gc", discrete= T, units = "nats", binning = "L2", biascorrection = "ml"){
  
  M <- dim(R)[2]
  if(is.null(M)) {
    M <- length(R)
    dim(R) <- c(1,M)
  }
  
  MI <- c()
  
  if(is.null(composite)){
    for (i in 1:M){
      
      r <- R[,i][!is.na(R[,i])]
      MI <- c(MI,Mutual_Info(data,S,r,estimator = estimator, discrete=discrete, units=units, binning = binning, biascorrection = biascorrection))
    }
  } else {
    for(i in 1:length(composite)){
      MI <- c(MI,Mutual_Info(data,S,R[,(i+(i-1)*(composite[1]-1)):(composite[1]+(i-1)*composite[length(composite)])],estimator = estimator, discrete=discrete, units=units, binning = binning, biascorrection = biascorrection))
      
    }
  }
  
  
  MMI <- min(MI)
  
  return(MMI)
}

cMMI <- function(data, S, R, estimator = "gc", discrete= T, units = "nats", binning = "L2", biascorrection = "ml"){
  
  M <- dim(R)[2]
  if(is.null(M)) {
    M <- length(R)
    dim(R) <- c(1,M)
  }
  
  MI <- c()
  for (i in 1:M){
    
    r <- R[,i][!is.na(R[,i])]
    MI <- c(MI,Mutual_Info(data,S,r,estimator = estimator, discrete=discrete, units=units, binning = binning, biascorrection = biascorrection))
  }
  
  for(i in 1:(M-1)) {
    if(M<2){break}
    for(j in (i+1):M){
      r1 <- R[,i][!is.na(R[,i])]
      r2 <- R[,j][!is.na(R[,j])]
      MI <- c(MI,Mutual_Info(data,r1,r2,estimator = estimator, discrete=discrete, units=units, binning = binning, biascorrection = biascorrection))
    }
  }
  cMMI <- min(MI)
  
  return(cMMI)
}

PID <- function(data, S, R=NULL, composite = F, R1 = NULL, R2 = NULL, redfun=MMI,estimator="gc", discrete = T, units="nats", binning = "L2", biascorrection = "ml") {
  
  if(composite){
    
    atoms <- c(paste0("{",R1,"}{",paste0(R2,collapse=","),"}"),
               paste0("{",R1,"}"),paste0("{",paste0(R2,collapse=","),"}"),
               paste0("{",paste0(c(R1,paste0(R2,collapse=",")),collapse= " "),"}"))
    
    k <- length(R1); l <- length(R2)
    
    I_partial <- c()
    
    I_partial <- c(I_partial, redfun(data,S,c(R1,R2), composite = c(k,l), estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)) # redundancy: I_partial(S;{1}{2})
    I_partial <- c(I_partial, (redfun(data,S,R1,composite=k,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 1: I_partial(S;{1})
    I_partial <- c(I_partial, (redfun(data,S,R2,composite=l,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 2: I_partial(S;{2})
    I_partial <- c(I_partial, (Mutual_Info(data,X=S,Y=c(R1,R2),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)
                               - sum(I_partial))) # synergy: I_partial(S;{12})
  }
  
  else {
    
    if(length(R)> 3 || length(R)<2) {
      stop("The number of source (response) variables should larger than 1 and smaller than 4.")
    }
    else if (length(R) == 2){
      atoms <- c(paste0("{",R[1],"}{",R[2],"}"),
                 paste0("{",R[1],"}"),paste0("{",R[2],"}"),
                 paste0("{",paste0(R,collapse= " "),"}"))
      
      I_partial <- c()
      
      I_partial <- c(I_partial, redfun(data,S,R,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)) # redundancy: I_partial(S;{1}{2})
      I_partial <- c(I_partial, (redfun(data,S,R[1],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 1: I_partial(S;{1})
      I_partial <- c(I_partial, (redfun(data,S,R[2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 2: I_partial(S;{2})
      I_partial <- c(I_partial, (Mutual_Info(data,X=S,Y=R,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)
                                 - sum(I_partial))) # synergy: I_partial(S;{12})
    }
    else {
      atoms <- c(paste0("{",R[1],"}{",R[2],"}{",R[3],"}"),
                 paste0("{",R[1],"}{",R[2],"}"),paste0("{",R[1],"}{",R[3],"}"),paste0("{",R[2],"}{",R[3],"}"),
                 paste0("{",R[1],"}{",paste0(R[2:3],collapse=" "),"}"),paste0("{",R[2],"}{",paste0(R[c(1,3)],collapse=" "),"}"),paste0("{",R[3],"}{",paste0(R[1:2],collapse=" "),"}"),
                 paste0("{",R[1],"}"),paste0("{",R[2],"}"),paste0("{",R[3],"}"),paste0("{",paste0(R[1:2],collapse=" "),"}{",paste0(R[c(1,3)],collapse=" "),"}{",paste0(R[c(2,3)],collapse=" "),"}"),
                 paste0("{",paste0(R[1:2],collapse=" "),"}{",paste0(R[c(1,3)],collapse=" "),"}"),paste0("{",paste0(R[1:2],collapse=" "),"}{",paste0(R[c(2,3)],collapse=" "),"}"),paste0("{",paste0(R[c(1,3)],collapse=" "),"}{",paste0(R[c(2,3)],collapse=" "),"}"),
                 paste0("{",paste0(R[1:2],collapse=" "),"}"),paste0("{",paste0(R[c(1,3)],collapse=" "),"}"),paste0("{",paste0(R[c(2,3)],collapse=" "),"}"),
                 paste0("{",paste0(R,collapse= " "),"}"))
      I_partial <- c()
      
      I_partial <- c(I_partial, redfun(data,S,R,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)) # 1. redundancy: I_partial(S;{1}{2}{3})
      I_partial <- c(I_partial, (redfun(data,S,R[1:2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 2. {1}{2}
      I_partial <- c(I_partial, (redfun(data,S,R[c(1,3)],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 3. {1}{3} 
      I_partial <- c(I_partial, (redfun(data,S,R[2:3],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 4. {2}{3}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],NA),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[1:3]))) # 5. {1}{23}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[2],NA),c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1,2,4)]))) # 6. {2}{13}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[3],NA),c(R[1],R[2])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1,3,4)]))) # 7. {3}{12}
      I_partial <- c(I_partial, (redfun(data,S,R[1],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,2,3,5)]))) # 8. {1}
      I_partial <- c(I_partial, (redfun(data,S,R[2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,2,4,6)]))) # 9. {2}
      I_partial <- c(I_partial, (redfun(data,S,R[3],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,3,4,7)]))) # 10. {3}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[2]),c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial))) # 11. {12}{13}{23}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[2]),c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:8,11)]))) # 12. {12}{13}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[2]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:7,9,11)]))) # 13. {12}{23}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:7,10,11)]))) # 14. {13}{23}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[2])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:9,11,12,13)]))) # 15. {12}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:8,10,11,12,14)]))) # 16. {13}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:7,9,10,11,13,14)]))) # 17. {23}
      I_partial <- c(I_partial, (redfun(data,S,cbind(c(R[1],R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial))) # 18. {123}
      
    }
  }
  
  pid <- data.frame(atoms, I_partial)
  
  return(pid)
}

MI_syn_red <- function(data,X,Y,kmax, synergy=TRUE, redundancy=TRUE, estimator='cov',units="nats", binning = "L2", discrete = T, biascorrection = "ml") { # A decomposition of pairwise MI
  
  variables <- colnames(data)
  
  mi_without <- Mutual_Info(data, X=X, Y=Y, estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  
  
  # redundancy
  if(redundancy){
    r <- kmax
    drivers <- variables[! variables %in% c(X,Y)]
    mi_red <- data.frame(Condition = character(), MutualInfo = numeric(), stringsAsFactors = FALSE)
    mi_red <- rbind(mi_red, data.frame(Condition = " ", MutualInfo = mi_without))
    
    drivers_red <- vector(mode = "character", length = 0)
    
    while (length(drivers)>0 && r>0) {
      
      mi <- vector(mode = "numeric", length = 0)
      
      for (var in drivers){
        conditionals <- c(drivers_red,var)
        mi <- c(mi,Mutual_Info(data,X=X,Y=Y,Z=conditionals,estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning))
      }
      
      min <- min(mi)
      condition <- drivers[which.min(mi)]
      
      drivers_red <- c(drivers_red, condition)
      
      mi_red <- rbind(mi_red, data.frame(Condition = paste(drivers_red,collapse = ", "), MutualInfo = min))
      
      drivers <- drivers[! drivers %in% condition]
      r <- r-1
    }
  } else {mi_red <- data.frame(Condition = character(), MutualInfo = numeric(), stringsAsFactors = FALSE)}
  
  # synergy
  if(synergy){
    s <- kmax
    drivers <- variables[! variables %in% c(X,Y)]
    mi_syn <- data.frame(Condition = character(), MutualInfo = numeric(), stringsAsFactors = FALSE)
    mi_syn <- rbind(mi_syn, data.frame(Condition = " ", MutualInfo = mi_without))
    
    drivers_syn <- vector(mode = "character", length = 0)
    
    while (length(drivers)>0 && s>0) {
      
      mi <- vector(mode = "numeric", length = 0)
      
      for (var in drivers){
        conditionals <- c(drivers_syn,var)
        mi <- c(mi,Mutual_Info(data,X=X,Y=Y,Z=conditionals,estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning))
      }
      
      max <- max(mi)
      condition <- drivers[which.max(mi)]
      
      drivers_syn <- c(drivers_syn, condition)
      
      mi_syn <- rbind(mi_syn, data.frame(Condition = paste(drivers_syn,collapse = ", "), MutualInfo = max))
      
      drivers <- drivers[! drivers %in% condition]
      
      s <- s - 1
    }
  } else {mi_syn <- data.frame(Condition = character(), MutualInfo = numeric(), stringsAsFactors = FALSE)}
  
  mi_red_syn <- list("redundancy"=mi_red, "synergy"=mi_syn)
  
  return(mi_red_syn)
  
}

MI_decomp_surr <- function(x, drivers_red, drivers_syn, synergy=TRUE, redundancy=TRUE, nsurr, kmax, surrtype="iaafft", estimator="cov", units="nats", discrete, binning, biascorrection,verbose=T) {
  
  # 1. redundancy
  if(redundancy){
    mi_red_surr <- matrix(nrow=nsurr,ncol=kmax)
    
    
    if(verbose){
      pb_red <- progress_bar$new(
        format = "  calculating redundant surrogates [:bar] :percent in :elapsed eta: :eta",
        total = kmax, clear = FALSE, width= 100)
    }
    
    
    for (dr in 3:(kmax+2)) {
      xx <- x[,drivers_red[1:(dr-1)]]
      for (k in 1:nsurr) {
        switch(surrtype,
               iaafft = {ys <- surr_iaafft(x[,drivers_red[dr]])},
               shuf = {ys <- surr_shuf(x[,drivers_red[dr]])})
        
        
        xxx <- cbind(xx,ys)
        colnames(xxx) <- drivers_red[1:dr]
        
        mi <- Mutual_Info(xxx, X=drivers_red[1],Y=drivers_red[2],Z=drivers_red[3:dr], estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
        
        mi_red_surr[k,dr-2] <- mi
      }
      if(verbose){
        pb_red$tick()
      }
    }
  } else{mi_red_surr <- matrix(nrow=nsurr,ncol=kmax)}
  
  # 2. synergy
  if(synergy){
    mi_syn_surr <- matrix(nrow=nsurr,ncol=kmax)
    
    if(verbose){
      pb_syn <- progress_bar$new(
        format = "  calculating synergistic surrogates [:bar] :percent in :elapsed eta: :eta",
        total = kmax, clear = FALSE, width= 100)
    }
    
    for (dr in 3:(kmax+2)) {
      xx <- x[,drivers_syn[1:(dr-1)]]
      for (k in 1:nsurr) {
        switch(surrtype,
               iaafft = {ys <- surr_iaafft(x[,drivers_syn[dr]])},
               shuf = {ys <- surr_shuf(x[,drivers_syn[dr]])})
        
        xxx <- cbind(xx,ys)
        colnames(xxx) <- drivers_syn[1:dr]
        
        mi <- Mutual_Info(xxx, X=drivers_syn[1],Y=drivers_syn[2],Z=drivers_syn[3:dr], estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
        
        mi_syn_surr[k,dr-2] <- mi
      }
      if(verbose){
        pb_syn$tick()
      }
    }
  } else{mi_syn_surr <- matrix(nrow=nsurr,ncol=kmax)}
  
  return(list(red_sur=mi_red_surr, syn_sur=mi_syn_surr))
}

MI_decomp <- function(data, target, driver, synergy=TRUE, redundancy=TRUE,estimator ="gc", nsurr = 10, alpha=0.05, surrtype="iaafft",kmax=NULL, units = "nats", binning = "L2", discrete = T, biascorrection = "ml",verbose=T,visual=F) {
  
  names <- colnames(data)
  
  if(is.null(kmax) || kmax>(length(names)-2)) {
    kmax <- length(names) - 2
    warning("kmax must be between 2 and number of variables - 2. Kmax is set at max.")
  }
  
  mi_red_syn <- MI_syn_red(data, X=target,Y=driver, synergy=synergy, redundancy=redundancy, kmax=kmax, estimator = estimator, units=units, binning = binning, discrete = discrete, biascorrection = biascorrection)
  
  syn <- mi_red_syn$synergy
  red <- mi_red_syn$redundancy
  
  labels_syn <- c("pairwise",unlist(strsplit(syn$Condition[length(names)-1],split = ", ")))
  labels_red <- c("pairwise",unlist(strsplit(red$Condition[length(names)-1],split = ", ")))
  
  drivers_red <- c(target,driver,unlist(strsplit(red$Condition[dim(red)[1]],split = ", ")))
  drivers_syn <- c(target,driver,unlist(strsplit(syn$Condition[dim(syn)[1]],split = ", ")))
  
  if(kmax == 0) {
    U <- red[1,2]; R <- 0; S <- 0
    atoms <- data.frame(cbind(U,R,S)) #total redundancy and synergy
    colnames(atoms) <- c("U","R","S")
    
    return(list(mi_red_syn=mi_red_syn,atoms=atoms))
  } else {
    
    mi_red_syn_surr <- MI_decomp_surr(data, drivers_red = drivers_red, drivers_syn = drivers_syn, synergy=synergy, redundancy=redundancy, nsurr = nsurr, surrtype=surrtype, kmax = kmax ,estimator = estimator,units=units, discrete = discrete, binning = binning, biascorrection = biascorrection,verbose=verbose)
    
    #calculating 'atoms' of the Mutual info decomposition
    U <- 0; R <- 0; S <- 0
    
    red_surr <- mi_red_syn_surr$red_sur
    syn_surr <- mi_red_syn_surr$syn_sur
    
    if(redundancy){
      
      
      p_red <- numeric(length = ncol(red_surr))
      
      for(i in 1:length(p_red)){
        p_red[i] <- (sum(red_surr[,i]<red[i+1,2]) + 1) / (nsurr+1)
      }
      
      p_red_adjusted <- p.adjust(p_red, method = "BH" )
      
      pvals_red <- data.frame(rbind(p_red,p_red_adjusted))
      colnames(pvals_red) <- setdiff(drivers_red,c(driver,target))
      
      #p_red <- (colSums(sweep(red_surr, 2, red[2:(ncol(red_surr)+1), 2], "<")) + 1) / nsurr
      
      
      for(i in 1:(kmax+1)){
        
        U <- red[i,2]
        R <- c(R,red[1,2] - U)
        if(isTRUE(red[i+1,2]>red[i,2]) || isTRUE(p_red_adjusted[i]>alpha) || i == (kmax + 1)){ # check value i because p-values start at first to condition on. There is no p-value for the pairwise MI
          multiplet_red <- drivers_red[1:(i+1)]
          break
        }
      }
      
      diff_red <- c(0,diff(R))
      red_tot <- sum(diff_red)
      contributions_red <- data.frame(cbind(multiplet_red,diff_red))
      
    } else {
      p_red <- numeric(length = ncol(red_surr))
      p_red_adjusted <- numeric(length = ncol(red_surr))
      pvals_red <- data.frame(rbind(p_red,p_red_adjusted))
      contributions_red <- data.frame()
      red_tot <- 0
    }
    
    if(synergy){
      
      p_syn <- numeric(length = ncol(syn_surr))
      
      for(i in 1:length(p_syn)){
        p_syn[i] <- (sum(syn_surr[,i]>syn[i+1,2]) + 1) / (nsurr+1)
      }
      
      p_syn_adjusted <- p.adjust(p_syn, method = "BH" )
      
      pvals_syn <- data.frame(rbind(p_syn,p_syn_adjusted))
      colnames(pvals_syn) <- setdiff(drivers_syn,c(driver,target))
      
      for(i in 1:(kmax+1)){
        S <- c(S,syn[i,2] - syn[1,2])
        
        if(isTRUE(syn[i+1,2]<syn[i,2]) || isTRUE(p_syn_adjusted[i]>alpha) || i == (kmax + 1)){
          multiplet_syn <- drivers_syn[1:(i+1)]
          break
        }
      }
      
      diff_syn <- c(0,diff(S))
      
      syn_tot <- sum(diff_syn)
      contributions_syn <- data.frame(cbind(multiplet_syn,diff_syn))
    } else {
      p_syn <- numeric(length = ncol(syn_surr))
      p_syn_adjusted <- numeric(length = ncol(syn_surr))
      pvals_syn <- data.frame(rbind(p_syn,p_syn_adjusted))
      contributions_syn <- data.frame()
      syn_tot <- S
    }
    
    atoms <- data.frame(cbind(U,red_tot,syn_tot)) #total redundancy and synergy
    colnames(atoms) <- c("U","R","S")
    
    if(visual){
      ymin_red <- min(red$MutualInfo,mi_red_syn_surr$red_sur)
      ymax_red <- max(red$MutualInfo,mi_red_syn_surr$red_sur)
      
      ymin_syn <- min(syn$MutualInfo,mi_red_syn_surr$syn_sur)
      ymax_syn <- max(syn$MutualInfo,mi_red_syn_surr$syn_sur)
      
      par(mfrow = c(2, 1))
      plot(red$MutualInfo, ylab = expression(I[gc]), xaxt='n', pch=19, xlab="", main=paste0('Mutual information: ',driver,' and ',target),ylim=c(ymin_red,ymax_red))
      for (i in 1:nsurr){
        points(2:(kmax+1),mi_red_syn_surr$red_sur[i,],col=rgb(1,0,0,alpha=0.25))
      }
      axis(1, at=seq(1,(kmax+1),by=1), labels=labels_red, las=2)
      #title(xlab="Variables conditioned on", line=4)
      title(sub = "Redundancy", line=-15)
      abline(v=1:(kmax+1),col = "lightgray", lty = "dotted")
      legend("topright",legend=c("Real data", "Surrogate data"),pch = c(19,1),
             col=c("black", "red"))
      
      plot(syn$MutualInfo, ylab = expression(I[gc]), xaxt='n', pch=19, xlab="",ylim=c(ymin_syn,ymax_syn) )
      for (i in 1:nsurr){
        points(2:(kmax+1),mi_red_syn_surr$syn_sur[i,],col=rgb(1,0,0,alpha=0.25))
      }
      axis(1, at=seq(1,(kmax+1),by=1), labels=labels_syn, las=2)
      title(xlab="Variables conditioned on", line=4)
      title(sub = "Synergy", line=-15)
      abline(v=1:(kmax+1),col = "lightgray", lty = "dotted")
    }
    
    return(list(mi_red_syn=mi_red_syn,mi_red_syn_surr=mi_red_syn_surr,atoms=atoms,pvals_red=pvals_red,pvals_syn=pvals_syn,contributions_red=contributions_red,contributions_syn=contributions_syn))
  }
}
