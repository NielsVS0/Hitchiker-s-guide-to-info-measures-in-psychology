# Calculation of various information theoretical measures using both a histogram based (HB) and covariance based (cov) esitmator.

source("C:path-to-script/gcmi.R")

discretize_equal_width <- function(values, bins) {
  discretized <- cut(values, breaks = bins, include.lowest = TRUE, labels = FALSE)
  return(discretized)
}

entropy_hist <- function(data , X, discrete, binning, units, biascorrection){
  
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
      
      switch(binning,
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
    
  }
  
  state_counts <- table(data_discretized)
  
  
  if(all(discrete)){
    
    bin_correction <- 1
    
    } else {
      
      bin_correction <- prod(bin_widths[!discrete])
    }
  
  #probabilities
  p_states <- state_counts / N
  
  if(biascorrection){
    K <- length(state_counts)
    entropy <- -sum(p_states * log2(p_states/bin_correction), na.rm = TRUE) + (K - 1)/(2*N*log(2))
  } else {
    entropy <- -sum(p_states * log2(p_states/bin_correction), na.rm = TRUE)
  }
  
  if(entropy < 0){
    warning("Entropy is negative")
    entropy <- 0
    }
  
  switch(units,
         bits = {return(entropy)},
         nats = {return(entropy*log(2))})
  
  
}

entropy_cov <- function(data, X, Y, units){
  
  X <- unique(X)
  Y <- unique(Y)
  
  if (is.null(X)){X <- colnames(data)}
  n <- length(X)
  
  if(is.null(Y)){
    
    cov_data <- cov(data[X])
    
    if(det(cov_data) == 0){ #singular
      entropy <- entropy_cov(data = data, X = X[1], Y = NULL, units = units)
    } else {
      entropy <- 0.5*n*log(2*pi*exp(1)) + 0.5*log(det(cov_data))
    }
    
  } else {
    cov_data <- cov(data)
    
    
    cov_XX <- cov_data[X, X, drop = FALSE]
    cov_XY <- cov_data[X, Y, drop = FALSE]
    cov_YY <- cov_data[Y, Y, drop = FALSE]
    
    res.cov <- cov_XX - cov_XY %*% solve(cov_YY) %*% t(cov_XY)
    
    entropy <- 0.5*n*log(2*pi*exp(1)) + 0.5*log(det(res.cov))
  }
  
  if(entropy < 0){
    warning("Entropy is negative")
    entropy <- 0
    }
  
  switch(units,
         nats = {return(entropy)},
         bits = {return(entropy*log2(exp(1)))})
}

Entropy <- function(data, X=NULL, Y=NULL, estimator="hist", discrete = T, biascorrection = T, units="nats", binning = "L2") {
  
  if (is.vector(data)) {
    data <- data.frame(data)
  }
  
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
  
  if (!estimator %in% c("hist","cov")) {
    stop("Invalid estimator. Choose 'hist' or 'cov'.")
  }
  
  switch(estimator,
         hist = {
           if(is.null(Y)){
             return(entropy_hist(data,X,discrete = discrete,biascorrection=biascorrection,units=units,binning=binning))
             } else {
               return(entropy_hist(data,c(X,Y),discrete = discrete,biascorrection=biascorrection,units=units,binning=binning) 
                      - entropy_hist(data,Y,discrete = discrete,biascorrection=biascorrection,units=units,binning=binning))
         }
           },
         cov = {
           return(entropy_cov(data,X,Y,units))
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
 
  
  p_states_X <- state_counts_X / N
  p_states_Y <- state_counts_Y / N
  p_states_joint <- state_counts_joint / N
  
  
  mi <- p_states_joint * log2(p_states_joint / (p_states_X %o% p_states_Y))
  
  
  if(biascorrection) {
    
    #correction <- - ((length(state_counts_X) - 1) + (length(state_counts_Y) - 1) - (length(state_counts_joint) - 1))/(2*N*log(2))
    
    correction <- (length(state_counts_X) - 1)*(length(state_counts_Y) - 1)/(2*N*log(2))
      
  } else {correction <- 0}
  

  
  switch(units,
         bits = {return(sum(mi, na.rm = TRUE) - correction)},
         nats = {return((sum(mi, na.rm = TRUE) - correction)*log(2))})
  
}
  
mutual_info_cov <- function (data, X, Y, Z, units) {
  
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
  
Mutual_Info <- function(data, X=NULL, Y=NULL, Z=NULL, estimator="hist", discrete = T, biascorrection = T, units="nats", binning = "L2"){
  
  
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
         cov = {return(mutual_info_cov(data = data, X = X, Y = Y, Z = Z, units = units))},
         gc = {return(mutual_info_gc(data = data, X = X, Y = Y, Z = Z, units = units))})
  
  
  
}

Co_Info <- function(data, X=NULL, Y=NULL, Z=NULL, estimator="hist", discrete=TRUE, biascorrection = T, binning = "L2", units = "nats") {
  
  data <- data.frame(data)
  
  if(is.null(X) || is.null(Y) || is.null(Z)){stop("X, Y or Z variables not specified")}
  
  mi <- Mutual_Info(data = data, X = X, Y = Y, estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  mi_cond <- Mutual_Info(data = data, X = X, Y = Y, Z = Z, estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  
  return(mi - mi_cond)
}

Total_Corr <- function(data, X=NULL, estimator="hist", binning = "L2", units = "nats", discrete=T, biascorrection = T) {
  
  if(is.null(X)){ #if the variables are not specified than the TC over all variables will be returned
    
    X <- colnames(data)
    
  }
  
  entropies <- c()
  
  for (variable in X) {
    entropies <- c(entropies, Entropy(data,X=variable,estimator = estimator, discrete=discrete, biascorrection = biascorrection, binning = binning, units=units))
  }
  joint_entropy <- Entropy(data,X=X,estimator=estimator, discrete=discrete, biascorrection=biascorrection, binning=binning, units=units)
  
  TC <- sum(entropies) - joint_entropy
  
  return(TC)
  
}

Resid_Entropy <- function(data, X, estimator="hist", binning = "L2", units = "nats",discrete=discrete, biascorrection=biascorrection) {
  
  environment_X <- colnames(data)[colnames(data) != X]
  
  R_X <- Entropy(data, X=X, Y=environment_X,estimator = estimator, discrete = discrete, biascorrection = biascorrection, units = units, binning = binning)
  
  return(R_X)
}

Dual_Total_Corr <- function(data, X=NULL, estimator="hist", binning = "L2", units = "nats", discrete=T, biascorrection=T) {
  
  if(is.null(X)){ #if the variables are not specified than the TC over all variables will be returned
    
    X <- colnames(data)
    
  }
  
  joint_entropy <- Entropy(data,X=X,estimator=estimator, binning=binning, units=units,discrete=discrete,biascorrection = biascorrection)
  
  resid_entropies <- c()
  for (variable in X) {
    resid_entropies <- c(resid_entropies, Resid_Entropy(data[X],X=variable, estimator=estimator, binning=binning, units=units,discrete=discrete,biascorrection=biascorrection))
  }
  
  DTC <- joint_entropy - sum(resid_entropies)
  
  return(DTC)
  
}

O_Info <- function(data, X=NULL, estimator = "hist", binning = "L2", units = "nats", discrete=T, biascorrection = T) {
  
  TC <- Total_Corr(data, X=X, estimator=estimator, binning=binning, units=units, discrete=discrete, biascorrection=biascorrection)
  DTC <- Dual_Total_Corr(data, X=X, estimator=estimator, binning=binning, units=units, discrete=discrete, biascorrection=biascorrection)
  
  O <- TC - DTC
  
  return(O)
}

IID <- function(data, order=3, estimator = "hist", binning = "L2", units = "nats", discrete=T, biascorrection=T, visual=F) {
  
  if (!is.numeric(order) || length(order) != 1 || order != as.integer(order)) {
    stop("'order' must be a single integer.")
  }
  
  if (order < 3) {
    stop("'order' must be at least 3.")
  }
  
  if (order > ncol(data)) {
    stop(paste0("'order' cannot be larger than the number of columns in 'data' (", ncol(data), ")."))
  }
  
  OI <- c()
  variables <- colnames(data)
  tuples <- combn(variables,order)
  for (i in 1:dim(tuples)[2]){
    OI <- c(OI,O_Info(data,X=tuples[,i],estimator=estimator, binning=binning, units=units, discrete = discrete, biascorrection = biascorrection)) 
  }
  
  tuple_string <- c()
  for (i in 1:dim(tuples)[2]){
    tuple_string <- c(tuple_string,paste0(tuples[,i],collapse = " "))
  }
  hoistructure <- data.frame(tuple_string, OI)
  
  if (visual) {
    barplot(hoistructure[order(hoistructure$OI),]$OI)
  }
  
  return(hoistructure)
}

MMI <- function(data, S, R, estimator = "gc", discrete= T, units = "nats", binning = "L2", biascorrection = T){
  
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
  MMI <- min(MI)
  
  return(MMI)
}

PID_MMI <- function(data, S, R, estimator="gc", discrete = T, units="nats", binning = "L2", biascorrection = T) {
  
  if(length(R)> 3 || length(R)<2) {
    stop("The number of source (response) variables should larger than 1 and smaller than 4.")
  }
  else if (length(R) == 2){
    atoms <- c(paste0("{",R[1],"}{",R[2],"}"),
               paste0("{",R[1],"}"),paste0("{",R[2],"}"),
               paste0("{",paste0(R,collapse= " "),"}"))
    
    I_partial <- c()
  
    I_partial <- c(I_partial, MMI(data,S,R,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)) # redundancy: I_partial(S;{1}{2})
    I_partial <- c(I_partial, (MMI(data,S,R[1],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 1: I_partial(S;{1})
    I_partial <- c(I_partial, (MMI(data,S,R[2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # Unique 2: I_partial(S;{2})
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
    
    I_partial <- c(I_partial, MMI(data,S,R,estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)) # 1. redundancy: I_partial(S;{1}{2}{3})
    I_partial <- c(I_partial, (MMI(data,S,R[1:2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 2. {1}{2}
    I_partial <- c(I_partial, (MMI(data,S,R[c(1,3)],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 3. {1}{3} 
    I_partial <- c(I_partial, (MMI(data,S,R[2:3],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-I_partial[1])) # 4. {2}{3}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],NA),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[1:3]))) # 5. {1}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[2],NA),c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1,2,4)]))) # 6. {2}{13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[3],NA),c(R[1],R[2])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1,3,4)]))) # 7. {3}{12}
    I_partial <- c(I_partial, (MMI(data,S,R[1],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,2,3,5)]))) # 8. {1}
    I_partial <- c(I_partial, (MMI(data,S,R[2],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,2,4,6)]))) # 9. {2}
    I_partial <- c(I_partial, (MMI(data,S,R[3],estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1,3,4,7)]))) # 10. {3}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial))) # 11. {12}{13}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:8,11)]))) # 12. {12}{13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:7,9,11)]))) # 13. {12}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection)-sum(I_partial[c(1:7,10,11)]))) # 14. {13}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:9,11,12,13)]))) # 15. {12}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:8,10,11,12,14)]))) # 16. {13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial[c(1:7,9,10,11,13,14)]))) # 17. {23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2],R[3])),estimator=estimator,discrete=discrete,units=units,binning=binning,biascorrection=biascorrection) - sum(I_partial))) # 18. {123}
    
  }

  pid_mmi <- data.frame(atoms, I_partial)
  
  return(pid_mmi)
}

