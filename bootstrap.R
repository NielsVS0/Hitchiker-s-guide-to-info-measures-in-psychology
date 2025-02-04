# bootstrapping

create_boots_data <- function(X,nboot,cov=FALSE) {
  
  # input:  X = n x m dataset
  #         nboot = number of bootstrap iterations
  #         cov = return covariance matrix
  
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  databoots <- vector("list", nboot)
  
  if (cov){ 
    
    for (i in 1:nboot) {
      idx <- sample.int(n,replace=TRUE)
      Xb <- X[idx,]
      
      databoots[[i]] <- cov(Xb)
    }
  }
  else{ 
    
    for (i in 1:nboot) {
      idx <- sample.int(n,replace=TRUE)
      
      databoots[[i]] <- X[idx,]
      
    }
  }
  
  return(databoots)
 
}

bootstrap_O_information <- function(data, var=NULL, n, nboot = 1000, alpha = 0.05, correction=F, estimator="mm", discrete=T, info=TRUE) {
  
  N <- nrow(data)  # Number of observations
  
  if(is.null(var)){
    tuples <- combn(colnames(data),n) # create tuples
  } else {
    tuples <- combn(var,n) # create tuples
  }
  ncomb <- ncol(tuples)  # Number of tuples
  
  # Store original and bootstrap results
  original_O <- numeric(ncomb)
  bootstrap_O <- matrix(0, nrow = ncomb, ncol = nboot)
  
  #create bootstrap samples
  
  boots_data <- create_boots_data(data,nboot)
  
  # Calculate O-information for each triple
  for (i in 1:ncomb) {
    if(info) {print(ncomb-i)}
    vars <- tuples[,i ]
  
    original_O[i] <- O_Info(data,vars,estimator = estimator,discrete = discrete)
    
    for (b in 1:nboot) {
      # Calculate O-information on bootstrap sample
      bootstrap_O[i, b] <- O_Info(boots_data[[b]],vars,estimator = estimator,discrete = discrete)
    }
  }
  
  p_values <- sapply(1:nrow(bootstrap_O), function(i) {
    boot_dist <- bootstrap_O[i, ]  # Bootstrap distribution for row i
    obs <- original_O[i]           # Corresponding observed value
    
    if (obs > 0) {
      # One-tailed test for positive observed values
      p_value <- (1 + sum(boot_dist < 0)) / (nboot + 1)
    } else if (obs < 0) {
      # One-tailed test for negative observed values
      p_value <- (1 + sum(boot_dist > 0)) / (nboot + 1)
    } else {
      # Observed value is exactly zero
      p_value <- 1
    }
    return(p_value)
  })
    
  
  if(correction){
  # adjust using FDR
    p <- p.adjust(p_values, method = "BH")
  }else{p <- p_values}
  # Determine significance
  significant <- p < alpha
  # Compute confidence intervals
  lower_bound <- apply(bootstrap_O, 1, function(boot_dist) {
    quantile(boot_dist, probs = alpha / 2)
  })
  upper_bound <- apply(bootstrap_O, 1, function(boot_dist) {
    quantile(boot_dist, probs = 1 - alpha / 2)
  })
  
  # Check if confidence interval includes zero
  includes_zero <- (lower_bound <= 0) & (upper_bound >= 0)
  
  # Return results
  return(data.frame(
    Tuple = apply(tuples, 2, paste0, collapse = " "),
    Original_O = original_O,
    P_Value = p,
    Significant = significant,
    Lower_CI = lower_bound,
    Upper_CI = upper_bound,
    Includes_Zero = includes_zero
  ))
}