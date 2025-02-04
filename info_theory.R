# Calculation of various information theoretical measures using both a histogram based (HB) and covariance based (cov) esitmator.
# This depends partly on the use of thepackage infotheo
library(infotheo)

source("C:path-to-script/gcmi.R")

discretize_equal_width <- function(values, bins) {
  discretized <- cut(values, breaks = bins, include.lowest = TRUE, labels = FALSE)
  return(discretized)
}

#histogram based estimator that can be used to estimate differential entropy
entropy_hist <- function(data , X, discrete=T, binning="L2", units = "bits"){
  
  N <- dim(data)[1]
  M <- length(X)
  
  if(!discrete){
    data_discretized <- data.frame(matrix(nrow=N,ncol= 1))
    colnames(data_discretized) <- X
    switch(binning,
           L2 = {
             if(IQR(data[,X])==0){
                 data_discretized[,1] <- discretize_equal_width(data[,X],bins=(N^(1/3)))
             } else {
                bin_width <- 2*IQR(data[,X])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,X]))/bin_width)
                 data_discretized[,1] <- discretize_equal_width(data[,X],bins=number_of_bins)
             }
           },
           cuberoot = {
             data_discretized[,1] <- discretize_equal_width(data[,X],bins=(N^(1/3)))
             bin_width <- diff(range(data[,X]))/length(table(data_discretized[,1]))
           })
    
    state_counts <- table(data_discretized[X])
    bin_widths <- rep(bin_width,length(state_counts))
    
  } else {
    state_counts <- table(data[X]) # Count occurrences of each unique state
    bin_widths <- rep(1,length(state_counts))
  }
  
  
  #print(bin_widths)
  
  
  #probabilities
  
  p_states <- state_counts / N
  
  entropy <- -sum(p_states * log2(p_states/bin_widths), na.rm = TRUE) #- (length(state_counts) - 1)/(2*N*log(2))
  
  
  switch(units,
         bits = {return(entropy)},
         nats = {return(entropy*log(2))})
  
  
}

# estimators below should only be used on discrete data
Entropy <- function(data, X, estimator="shrink", discretized=TRUE, binning = "L2", units= "nats", bins = NA) {
  # X can be a single or multiple strings for univariate or multivariate entropy respectively
  # possible estimators:  "cov" = covariance based estimation
  #                       "shrink" = James-stein shrinkage estimator
  #                       "emp" = ML estimator
  
  M <- length(X)
  N <- dim(data)[1]
  if(estimator == "cov") {
    switch(units, nats = {base = exp(1)}, bits = {base = 2})
    return(0.5*log(det(cov(data[X]))*(2*pi*exp(1))^M,base=base)) #wrong: base change only affects log and not the 0.5 factor
  }
  else {
    if(!discretized){
      data_discretized <- data.frame(matrix(nrow=N,ncol= M))
      colnames(data_discretized) <- X
      switch(binning,
             L2 = {
               for (i in 1:M){
                 if(IQR(data[,X[i]])==0){
                   data_discretized[,i] <- discretize(data[,X[i]])
                 } else {
                   bin_width <- 2*IQR(data[,X[i]])/N^(1/3)
                   number_of_bins <- round(diff(range(data[,X[i]]))/bin_width)
                   data_discretized[,i] <- discretize(data[,X[i]],nbins=number_of_bins)
                   }
              }
            },
            cuberoot = {
               for (i in 1:M){
                 data_discretized[,i] <- discretize(data[,X[i]])
              }
            },
            custom = {
              for (i in 1:M){
                data_discretized[,i] <- discretize(data[,X[i]],nbins=bins)
              }
            })
    } else {data_discretized <- data[,X]}
    switch(units, 
           nats = {return(entropy(data_discretized, method=estimator))}, 
           bits = {return(natstobits(entropy(data_discretized, method=estimator)))})
    
    
  }
}

Mutual_Info <- function(data, X, Y, estimator="shrink", discretized=TRUE, binning = "L2", units = "nats") {
  
  MX <- length(X)
  MY <- length(Y)
  N <- dim(data)[1]
  
  if(estimator == "cov") {
    switch(units, nats = {base = exp(1)}, bits = {base = 2})
    X_cov <- cov(data[X])
    Y_cov <- cov(data[Y])
    X_pcov_Y <- X_cov - cov(data[X],data[Y])%*%solve(Y_cov)%*%cov(data[Y],data[X]) #partial covariance 
    return(0.5*log(det(X_cov)/det(X_pcov_Y),base=base))
  }
  else if (estimator == "gc"){
    switch(units, 
           nats = {return(log(2)*gcmi_cc(t(data[X]),t(data[Y])))}, 
           bits = {return(gcmi_cc(t(data[X]),t(data[Y])))})
  }
  else {
    if(!discretized){
      X_discretized <- data.frame(matrix(nrow=N,ncol= MX))
      colnames(X_discretized) <- X
      Y_discretized <- data.frame(matrix(nrow=N,ncol= MY))
      colnames(Y_discretized) <- Y
      switch(binning,
             L2 = {
              for (i in 1:MX){
                bin_width <- 2*IQR(data[,X[i]])/N^(1/3)
                number_of_bins <- round(diff(range(data[,X[i]]))/bin_width)
                X_discretized[,i] <- discretize(data[X[i]],nbins=number_of_bins)
              }
              for (i in 1:MY){
                 bin_width <- 2*IQR(data[,Y[i]])/N^(1/3)
                number_of_bins <- round(diff(range(data[,Y[i]]))/bin_width)
                Y_discretized[,i] <- discretize(data[Y[i]],nbins=number_of_bins)
              }
            },
             cuberoot = {
              for (i in 1:MX){
                 X_discretized[,i] <- discretize(data[X[i]])
               }
              for (i in 1:MY){
                Y_discretized[,i] <- discretize(data[Y[i]])
               }
             })
    } else {
      X_discretized <- data[,X]
      Y_discretized <- data[,Y]
      }
    switch(units, 
           nats = {return(mutinformation(X_discretized, Y_discretized, method=estimator))}, 
           bits = {return(natstobits(mutinformation(X_discretized, Y_discretized, method=estimator)))})
  }
}

Co_Info <- function(data, X, Y, Z, estimator="shrink", discretized=TRUE, binning = "L2", units = "nats") {
  
  MX <- length(X)
  MY <- length(Y)
  MZ <- length(Z)
  N <- dim(data)[1]
  
  if(estimator == "cov") {
    switch(units, nats = {base = exp(1)}, bits = {base = 2})
    X_cov <- cov(data[X])
    Y_cov <- cov(data[Y])
    Z_cov <- cov(data[Z])
    YZ_cov <- cov(data[c(Y,Z)])
    X_pcov_Y <- X_cov - cov(data[X],data[Y])%*%solve(Y_cov)%*%cov(data[Y],data[X]) #partial covariance X|Y
    X_pcov_Z <- X_cov - cov(data[X],data[Z])%*%solve(Z_cov)%*%cov(data[Z],data[X]) #partial covariance X|Z
    X_pcov_YZ <- X_cov - cov(data[X],data[c(Y,Z)])%*%solve(YZ_cov)%*%cov(data[c(Y,Z)],data[X]) #partial covariance X|YZ
    return(0.5*log((det(X_cov)*det(X_pcov_YZ))/(det(X_pcov_Y)*det(X_pcov_Z)),base=base))
  }
  else if (estimator == "gc"){
    switch(units, 
           nats = {return(log(2)*(gcmi_cc(t(data[X]),t(data[Y]))-gccmi_ccc(t(data[X]),t(data[Y]),t(data[Z]))))}, 
           bits = {return(gcmi_cc(t(data[X]),t(data[Y]))-gccmi_ccc(t(data[X]),t(data[Y]),t(data[Z])))})
  }
  else {
    if (!discretized){
      X_discretized <- data.frame(matrix(nrow=N,ncol= MX))
      colnames(X_discretized) <- X
      Y_discretized <- data.frame(matrix(nrow=N,ncol= MY))
      colnames(Y_discretized) <- Y
      Z_discretized <- data.frame(matrix(nrow=N,ncol= MZ))
      colnames(Z_discretized) <- Z
      switch(binning,
             L2 = {
               for (i in 1:MX){
                 bin_width <- 2*IQR(data[,X[i]])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,X[i]]))/bin_width)
                 X_discretized[,i] <- discretize(data[X[i]],nbins=number_of_bins)
               }
               for (i in 1:MY){
                 bin_width <- 2*IQR(data[,Y[i]])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,Y[i]]))/bin_width)
                 Y_discretized[,i] <- discretize(data[Y[i]],nbins=number_of_bins)
               }
               for (i in 1:MZ){
                 bin_width <- 2*IQR(data[,Z[i]])/N^(1/3)
                 number_of_bins <- round(diff(range(data[,Z[i]]))/bin_width)
                 Z_discretized[,i] <- discretize(data[Z[i]],nbins=number_of_bins)
               }
               #data_discretized <- cbind.data.frame(X_discretized,Y_discretized,Z_discretized)
             },
             cuberoot = {
               for (i in 1:MX){
                 X_discretized[,i] <- discretize(data[X[i]])
               }
               for (i in 1:MY){
                 Y_discretized[,i] <- discretize(data[Y[i]])
               }
               for (i in 1:MZ){
                 Z_discretized[,i] <- discretize(data[Z[i]])
               }
               #data_discretized <- cbind.data.frame(X_discretized,Y_discretized,Z_discretized)
            })
    } else {
      X_discretized <- data[,X]
      Y_discretized <- data[,Y]
      Z_discretized <- data[,Z]
    }
    
    CI <- mutinformation(X_discretized,Y_discretized,method=estimator) + mutinformation(X_discretized,Z_discretized,method=estimator) - mutinformation(X_discretized,cbind(Y_discretized,Z_discretized),method=estimator)
    switch(units, 
           nats = {return(CI)}, 
           bits = {return(natstobits(CI))})
  }
}

Total_Corr <- function(data, X, estimator="shrink", binning = "L2", units = "nats",discrete=T) {
  
  entropies <- c()
  for (variable in X) {
    entropies <- c(entropies, Entropy(data,variable,estimator = estimator, binning = binning, units=units,discretized=discrete))
  }
  joint_entropy <- Entropy(data,X,estimator=estimator, binning=binning, units=units,discretized=discrete)
  
  TC <- sum(entropies) - joint_entropy
  
  return(TC)
  
}

Resid_Entropy <- function(data, X, estimator="shrink", binning = "L2", units = "nats",discrete=discrete) {
  
  whole <- colnames(data) # whole system, i.e. all variables
  environmentX <- whole[whole != X] # environment of X, i.e. all variables besides X
  
  joint_entropy <- Entropy(data,whole,estimator=estimator, binning=binning, units=units,discretized=discrete)
  environment_entropy <- Entropy(data,environmentX,estimator=estimator, binning=binning, units=units,discretized=discrete)
  
  R_X <- joint_entropy - environment_entropy
  
  return(R_X)
}

Dual_Total_Corr <- function(data, X, estimator="shrink", binning = "L2", units = "nats",discrete=T) {
  
  joint_entropy <- Entropy(data,X,estimator=estimator, binning=binning, units=units,discretized=discrete)
  
  resid_entropies <- c()
  for (variable in X) {
    resid_entropies <- c(resid_entropies, Resid_Entropy(data[X],variable, estimator=estimator, binning=binning, units=units,discrete=discrete))
  }
  
  DTC <- joint_entropy - sum(resid_entropies)
  
  return(DTC)
  
}

O_Info <- function(data, X, estimator = "shrink", binning = "L2", units = "nats", discrete=T) {
  
  TC <- Total_Corr(data,X,estimator=estimator, binning=binning, units=units,discrete=discrete)
  DTC <- Dual_Total_Corr(data,X,estimator=estimator, binning=binning, units=units,discrete=discrete)
  
  O <- TC - DTC
  
  return(O)
}

IID <- function(data, order, estimator = "shrink", binning = "L2", units = "nats", barplot=F) {
  
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
    OI <- c(OI,O_Info(data,tuples[,i],estimator=estimator, binning=binning, units=units)) 
  }
  
  tuple_string <- c()
  for (i in 1:dim(tuples)[2]){
    tuple_string <- c(tuple_string,paste0(tuples[,i],collapse = " "))
  }
  hoistructure <- data.frame(tuple_string, OI)
  
  if (barplot) {
    barplot(hoistructure[order(hoistructure$OI),]$OI)
  }
  
  return(hoistructure)
}

MMI <- function(data, S, R, estimator = "gc", discretized= T, units = "nats"){
  M <- dim(R)[2]
  if(is.null(M)) {
    M <- length(R)
    dim(R) <- c(1,M)
    }
  MI <- c()
  for (i in 1:M){
    r <- R[,i][!is.na(R[,i])]
    MI <- c(MI,Mutual_Info(data,S,r,estimator = estimator, discretized=discretized, units=units))
  }
  MMI <- min(MI)
  
  return(MMI)
}

PID_MMI <- function(data,S,R,estimator="gc", discretized = T,units="nats") {
  if(length(R)> 3 || length(R)<2) {
    stop("The number of source variables should larger than 1 and smaller than 4.")
  }
  else if (length(R) == 2){
    atoms <- c(paste0("{",R[1],"}{",R[2],"}"),
               paste0("{",R[1],"}"),paste0("{",R[2],"}"),
               paste0("{",paste0(R,collapse= " "),"}"))
    I_partial <- c()
  
    I_partial <- c(I_partial, MMI(data,S,R,estimator=estimator,discretized=discretized,units=units)) # redundancy: I_partial(S;{1}{2})
    I_partial <- c(I_partial, (MMI(data,S,R[1],estimator=estimator,discretized=discretized,units=units)-I_partial[1]))
    I_partial <- c(I_partial, (MMI(data,S,R[2],estimator=estimator,discretized=discretized,units=units)-I_partial[1]))
    I_partial <- c(I_partial, (Mutual_Info(data,S,R,estimator=estimator,discretized=discretized,units=units) - sum(I_partial)))
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
    
    I_partial <- c(I_partial, MMI(data,S,R,estimator=estimator,discretized=discretized,units=units)) # 1. redundancy: I_partial(S;{1}{2}{3})
    I_partial <- c(I_partial, (MMI(data,S,R[1:2],estimator=estimator,discretized=discretized,units=units)-I_partial[1])) # 2. {1}{2}
    I_partial <- c(I_partial, (MMI(data,S,R[c(1,3)],estimator=estimator,discretized=discretized,units=units)-I_partial[1])) # 3. {1}{3} 
    I_partial <- c(I_partial, (MMI(data,S,R[2:3],estimator=estimator,discretized=discretized,units=units)-I_partial[1])) # 4. {2}{3}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],NA),c(R[2],R[3])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[1:3]))) # 5. {1}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[2],NA),c(R[1],R[3])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[c(1,2,4)]))) # 6. {2}{13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[3],NA),c(R[1],R[2])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[c(1,3,4)]))) # 7. {3}{12}
    I_partial <- c(I_partial, (MMI(data,S,R[1],estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1,2,3,5)]))) # 8. {1}
    I_partial <- c(I_partial, (MMI(data,S,R[2],estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1,2,4,6)]))) # 9. {2}
    I_partial <- c(I_partial, (MMI(data,S,R[3],estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1,3,4,7)]))) # 10. {3}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discretized=discretized,units=units)-sum(I_partial))) # 11. {12}{13}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[1],R[3])),estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1:8,11)]))) # 12. {12}{13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2]),c(R[2],R[3])),estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1:7,9,11)]))) # 13. {12}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[3]),c(R[2],R[3])),estimator=estimator,discretized=discretized,units=units)-sum(I_partial[c(1:7,10,11)]))) # 14. {13}{23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[c(1:9,11,12,13)]))) # 15. {12}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[3])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[c(1:8,10,11,12,14)]))) # 16. {13}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[2],R[3])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial[c(1:7,9,10,11,13,14)]))) # 17. {23}
    I_partial <- c(I_partial, (MMI(data,S,cbind(c(R[1],R[2],R[3])),estimator=estimator,discretized=discretized,units=units) - sum(I_partial))) # 18. {123}
    
  }
  
  pid_mmi <- data.frame(atoms, I_partial)
  
  return(pid_mmi)
}

