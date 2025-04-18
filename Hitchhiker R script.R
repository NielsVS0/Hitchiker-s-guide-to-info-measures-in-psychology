##############################################################################################
############            Code to reproduce pictures and examples from              ############
############ A hitchhiker's guide to information-theoretic measures in psychology ############
############               N. Van santen, Y. Roseel, D. Marinazzo                 ############
##############################################################################################

##### initialization #####
setwd("your working directory")

# necessary package; instal with install.packages("package") if necessary
library(mvtnorm)
library(qgraph)
library(ppcor)

# load the dataset (Empathy data can be found at https://data.mendeley.com/datasets/b8d5n5h3gc/1)
data <- read.csv(file="C:\\path-to-data\\data.csv",header = TRUE,sep = ";")

# names
names<- c("1FS", "2EC", "3PT_R", "4EC_R", "5FS", "6PD", "7FS_R", 
          "8PT","9EC", "10PD", "11PT", "12FS_R", "13PD_R", "14EC_R", "15PT_R", 
          "16FS", "17PD", "18EC_R", "19PD_R", "20EC", "21PT", "22EC", "23FS", 
          "24PD", "25PT", "26FS", "27PD", "28PT")

# column names
colnames(data) <- names

# load scripts with the functions we need to calculate the info measures
source("C:/path-to-script/info_theory.R") # information theoretic measures
  

##### Uncertainty and information #####

#### Figure 1: Entropy for a binary variable ####
p <- seq(from = 0, to=1,by=0.01)

entropy <- -p*log(p,base=2) - (1-p)*log(1-p,base=2) # entropy for a single binary variable in bits

plot(p,entropy,pch=19,ylab="H(p) (bits)",xlab="p",cex.lab=1.3, cex.axis=1.5)

#### Figure 2: comparison between entropy, MAD, and variance ####

# Initialize dataframe
df_entropy<- data.frame(
  Variable = names,
  Mean = sapply(names, function(var) mean(data[[var]])),
  Variance = sapply(names, function(var) var(data[[var]])),
  MAD = sapply(names, function(var) mean(abs(data[[var]] - median(data[[var]])))),
  Entropy_Hist = sapply(names, function(var) Entropy(data, var, discrete = TRUE, estimator = "hist")),
  Entropy_par = sapply(names, function(var) Entropy(data, var, estimator = "cov")),
  stringsAsFactors = FALSE
)

# Plot MAD vs Entropy_Hist
plot(df_entropy$MAD, df_entropy$Entropy_Hist, 
     ylab = "Entropy (nats)", xlab = "Mean Absolute Deviation", 
     cex.lab = 1.5, cex.axis = 1.3, pch = 19, 
     xlim = range(df_entropy$MAD))

highlight_idx <- c(4, 12, 18, 25)
text(df_entropy$MAD[highlight_idx], df_entropy$Entropy_Hist[highlight_idx], 
     labels = df_entropy$Variable[highlight_idx], pos = c(4, 3, 2, 3), cex = 1.5)

# Overlay Variance vs Entropy_Hist
par(new = TRUE)
plot(df_entropy$Variance, df_entropy$Entropy_Hist, col = "red", 
     ann = FALSE, yaxt = "n", xaxt = "n", pch = 19)
text(df_entropy$Variance[highlight_idx], df_entropy$Entropy_Hist[highlight_idx], 
     labels = df_entropy$Variable[highlight_idx], pos = c(2, 4, 4, 3), cex = 1.5, col = "red")
       
axis(3, at = pretty(df_entropy$Variance), labels = round(pretty(df_entropy$Variance), 2), col = "red", col.axis = "red")
mtext("Variance", side = 3, line = 3, cex = 1.5, col = "red")
#legend("topleft", pch = 19, col = c("black", "red"), 
#       legend = c("Entropy vs MAD", "Entropy vs Variance"), cex = 1.5)

# Histograms for 12FS_R and 25PT
par(mfrow=c(2,1))
hist(data[["12FS_R"]],breaks=seq(-0.5, 4.5, by=1), main="12FS_R",xlab='',cex.lab=1.5, cex.axis=1.3)
abline(v = c(median(data[["12FS_R"]]),mean(data[["12FS_R"]])), col=c('blue','red'), lwd = 2, lty = 'dashed')
legend("topleft", lty='dashed', col=c('blue','red'),legend=c('median','mean'),cex=1.5)
hist(data[["25PT"]],breaks=seq(-0.5, 4.5, by=1), main="25PT",xlab='scores',cex.lab=1.5, cex.axis=1.3)
abline(v = c(median(data[["25PT"]]),mean(data[["25PT"]])), col=c('blue','red'), lwd = 2, lty = 'dashed')

#### Figure 3: Entropy of different distributions ####

# Set seed for reproducibility
set.seed(123)

# Number of samples
n <- 100000

# Generate distributions
dist1 <- c(rnorm(n / 2, mean = -2, sd = 0.6577446), rnorm(n / 2, mean = 2, sd = 0.6577446)) # Bimodal
dist2 <- c(rnorm(n / 2, mean = -1.5, sd = 0.6577446), rnorm(n / 2, mean = 1.5, sd = 0.6577446)) # Less Bimodal
dist3 <- c(rnorm(n / 2, mean = -1, sd = 0.6577446), rnorm(n / 2, mean = 1, sd = 0.6577446)) # mimick uniform with normal
dist4 <- c(rnorm(n / 2, mean = -0.5, sd = 0.6577446), rnorm(n / 2, mean = 0.5, sd = 0.6577446)) # Near Gaussian
dist5 <- rnorm(n, mean = 0, sd = 0.6577446) # Gaussian
distT <- rt(n, df = 3) * 0.38 # Heavy-tailed
distE <- rexp(n, rate = 1 / 0.6577446) # Skewed

# Store distributions in a data frame
distr <- data.frame(dist1, dist2, dist3, dist4, dist5, distT, distE)

# Compute variances
distVars <- sapply(distr, var)

# Compute entropies
distEnt <- sapply(names(distr), function(var) Entropy(distr, var, discrete = FALSE, units = "nats", estimator = "hist"))

# Define function for histogram plotting with dynamic labels
plot_hist <- function(data, breaks, var_val, ent_val) {
  hist(data, breaks = breaks, xlim = c(-4, 4), col = "skyblue", 
       main = bquote(S^2 * " = " * .(round(var_val, 2)) * "; " ~ H(X) * " = " * .(round(ent_val, 2))), 
       xlab = "", cex.lab = 1.4, cex.main = 1.9, cex.axis = 1.3)
}

# Plot distributions with dynamically generated expressions
par(mfrow = c(7, 1), mar = c(4, 4, 2, 1))

plot_hist(dist1, 63, distVars["dist1"], distEnt["dist1"])
plot_hist(dist2, 82, distVars["dist2"], distEnt["dist2"])
plot_hist(dist3, 112, distVars["dist3"], distEnt["dist3"])
plot_hist(dist4, 132, distVars["dist4"], distEnt["dist4"])
plot_hist(dist5, 148, distVars["dist5"], distEnt["dist5"])
plot_hist(distT, 2578, distVars["distT"], distEnt["distT"])
plot_hist(distE, 225, distVars["distE"], distEnt["distE"])

# Set x-axis label for the last plot
mtext("Value", side = 1, line = 2, cex = 1.4)

# Print variance and entropy results
distVars
distEnt

dev.off()

#### Figure 5: MI versus cor ####

# plot theoretical relationship between MI and cor for normally distributed variables
cor_theory <- seq(from=0,to=1,by=0.0001)
MI_cov_theory <- -0.5*log(1-(cor_theory*cor_theory))
plot(cor_theory,MI_cov_theory,type = "line",ylab = "I(X;Y) in nats",xlab=expression(rho),cex.lab=1.4,cex.axis=1.3)
abline(v = 1, col = "red", lty = 2)
text(0.55, 4, labels = expression(I(X*";"*Y) %->% infinity ~ "as" ~ rho %->% 1), pos = 4, cex=2)

# plot cor vs MI in empathy data

n <-length(names)  # Automatically detects number of variables

mi_gc <- matrix(NA, nrow = n, ncol = n)
mi_HB <- matrix(NA, nrow = n, ncol = n)
mi_cov <- matrix(NA, nrow = n, ncol = n)

for (i in seq_along(names)) {
  for (j in seq_along(names)) {
    mi_HB[i, j] <- Mutual_Info(data, names[i], names[j], estimator = "hist", discrete = TRUE)
    mi_cov[i, j] <- Mutual_Info(data, names[i], names[j], estimator = "cov")
    mi_gc[i, j] <- Mutual_Info(data, names[i], names[j], estimator = "gc")
  }
}

# Compute Spearman correlation matrix
cor_matrix <- cor(data, method = "spearman")

# Assign row and column names to MI matrices
colnames(mi_gc) <- rownames(mi_gc) <- colnames(cor_matrix)
colnames(mi_HB) <- rownames(mi_HB) <- colnames(cor_matrix)
colnames(mi_cov) <- rownames(mi_cov) <- colnames(cor_matrix)

# Convert matrices to data frames and reshape
reshape_mi <- function(mi_matrix) {
  df <- as.data.frame(as.table(mi_matrix))
  df <- df[df$Var1 != df$Var2, ]  # Remove self-correlations
  df <- df[!duplicated(t(apply(df[, 1:2], 1, sort))), ]  # Remove duplicate pairs
  return(df)
}

cor_df <- reshape_mi(cor_matrix)
mi_gc_df <- reshape_mi(mi_gc)
mi_HB_df <- reshape_mi(mi_HB)
mi_cov_df <- reshape_mi(mi_cov)

# Merge data frames
merged_df <- Reduce(function(x, y) merge(x, y, by = c("Var1", "Var2")), 
                    list(cor_df, mi_gc_df, mi_HB_df, mi_cov_df))

colnames(merged_df) <- c("Variable1", "Variable2", "Correlation", "MI_gc", "MI_HB", "MI_cov")


# Display the result
print(merged_df)

# Plot MI vs Correlation
plot(merged_df$Correlation, merged_df$MI_gc, xlab = expression(rho), 
     ylab = "I(X;Y) in nats", col = 'green', cex.lab = 1.4, cex.axis = 1.3, 
     pch = 19, ylim = range(merged_df[4:6]))

points(merged_df$Correlation, merged_df$MI_cov, col = rgb(1, 0, 0, alpha = 0.5), pch = 19)
points(merged_df$Correlation, merged_df$MI_HB, col = rgb(0, 0, 1, alpha = 0.5), pch = 19)

legend("topleft", pch = 19, col = c("green", "blue", "red"), 
       legend = c("Gaussian Copula", "Histogram-based", "Covariance-based"), cex = 2)

#### Figure 6: Comparing MI with correlation for different scatter plots ####

# Function to generate nonlinear data patterns
generate_nonlinear <- function(pattern, N = 10000, sd = 0.1) {
  set.seed(123)
  
  result <- switch(pattern,
                   quadratic = {
                     x <- seq(-4, 4, length.out = N) 
                     y <- x^2 + rnorm(N, sd = sd)
                     cbind(x, y)
                   },
                   sine_squared = {
                     x <- seq(-pi, pi, length.out = N) + rnorm(N, sd = sd)
                     y <- sin(x)^2 + rnorm(N, sd = sd)
                     cbind(x, y)
                   },
                   circle = {
                     theta <- seq(0, 4 * pi, length.out = N)
                     x <- 4 * cos(theta) + rnorm(N, sd = sd)
                     y <- 4 * sin(theta) + rnorm(N, sd = sd)
                     cbind(x, y)
                   },
                   diamond_filled = {
                     x <- runif(N, -4, 4) + rnorm(N, sd = sd)
                     y <- runif(N, -4, 4) + rnorm(N, sd = sd)
                     inside <- abs(x) + abs(y) <= 4
                     cbind(x[inside], y[inside])
                   },
                   four_clusters = {
                     centers <- matrix(c(-2, 2, 2, 2, -2, -2, 2, -2), ncol = 2, byrow = TRUE)
                     generate_cluster <- function(center, n, sd = 0.5) {
                       matrix(rnorm(2 * n, mean = rep(center, each = n), sd = sd), ncol = 2)
                     }
                     do.call(rbind, lapply(1:4, function(i) generate_cluster(centers[i,], round(N / 4))))
                   })
  
  colnames(result) <- c("V1", "V2")
  as.data.frame(result)
}

# Function to generate linear data with specified correlation
generate_linear <- function(correlation, N = 10000, var1 = 1, var2 = 1) {
  sigma <- matrix(c(var1, correlation * sqrt(var1 * var2), correlation * sqrt(var1 * var2), var2), 
                  nrow = 2, byrow = TRUE)
  result <- mvtnorm::rmvnorm(n = N, mean = c(0, 0), sigma = sigma)
  colnames(result) <- c("V1", "V2")
  as.data.frame(result)
}

### Data Generation ###

# Linear datasets
linear_correlations <- c(1, 0.5, 0, -0.5, -1)
linear_datasets <- lapply(linear_correlations, generate_linear, N = 10000)
names(linear_datasets) <- paste0("lin", c(1, "05", "0", "min05", "min1"))

# Nonlinear datasets
nonlinear_patterns <- c("quadratic", "sine_squared", "circle", "four_clusters", "diamond_filled")
nonlinear_datasets <- lapply(nonlinear_patterns, generate_nonlinear, N = 10000)
names(nonlinear_datasets) <- c("quadratic", "sinesquared", "circle", "four_clusters", "diamondf")

# Combine all datasets
shapes_data <- c(linear_datasets, nonlinear_datasets)
shapes_names <- names(shapes_data)

### Compute Correlation and Mutual Information ###

compute_assoc <- function(data) {
  list(
    correlation = cor(data$V1, data$V2),
    MI_HB = Mutual_Info(data, "V1", "V2",discrete = F),
    MI_COV = Mutual_Info(data, "V1", "V2", estimator = "cov"),
    MI_GC = Mutual_Info(data, "V1", "V2", estimator = "gc")
  )
}

# Compute for all datasets
assoc <- data.frame(relation = shapes_names, t(sapply(shapes_data, compute_assoc)))

### Plot Results ###

dev.off()
par(mfrow = c(2, 5))

for (i in 1:length(shapes_names)) {
  plot(shapes_data[[i]], ann = FALSE, axes = FALSE)
  
  # Extract numeric values safely
  rho <- round(as.numeric(assoc$correlation[i]), 2)
  I_HB <- round(as.numeric(assoc$MI_HB[i]), 2)
  I_COV <- ifelse(is.infinite(as.numeric(assoc$MI_COV[i])), "∞", round(as.numeric(assoc$MI_COV[i]), 2))
  I_GC <- round(as.numeric(assoc$MI_GC[i]), 2)
  
  
  # Add labels correctly
  mtext(paste("ρ =", rho, ";  I_HB =", I_HB), side = 3, cex = 1.5)
  mtext(paste("I_COV =", I_COV, ";  I_GC =", I_GC), side = 1, cex = 1.5)
}





##### Higher order information #####

#### Figure 7; tables C1, C2, C3: co-information and partial correlations ####

compute_info_and_pcor <- function(x, y, z) { #conditioning happens on the third variable
  df <- data.frame(x, y, z)
  co_inf <- Co_Info(df,"x","y","z",estimator="cov")
  
  cov_matrix <- cov(df)
  cor_matrix <- cov2cor(cov_matrix)
  pcor_matrix <- pcor(df)
  
  list(co_inf = co_inf, cor_matrix = cor_matrix, pcor_matrix = pcor_matrix)
}

# Generate data and compute metrics for each scenario
set.seed(42)  # For reproducibility

# Independency 
x_i = rnorm(10000)
y_i = x_i + rnorm(10000)
z_i = rnorm(10000)
data_i <- list(x=x_i,y=y_i,z=z_i)

info_i <- compute_info_and_pcor(data_i$x,data_i$y,data_i$z) # table C1

# Redundancy
z_r = rnorm(10000)
x_r = z_r + rnorm(10000)
y_r = 0.5 * (x_r + 2*z_r) + rnorm(10000)
data_r <- list(x=x_r,y=y_r,z=z_r)

info_r <- compute_info_and_pcor(data_r$x,data_r$y,data_r$z) # table C2

# Synergy
x_s = rnorm(10000)
y_s = 0.5*x_s + rnorm(10000)
z_s = x_s + y_s + 0.5*rnorm(10000)
data_s <- list(x=x_s,y=y_s,z=z_s)

info_s <- compute_info_and_pcor(data_s$x, data_s$y, data_s$z) #  table C3

# Generate plots (figure 7)
#pdf("triplet_cor.pdf", width = 12, height = 3)
par(mfrow = c(1, 3))
qgraph(info_i$cor_matrix, graph="cor", labels=c("X", "Y", "Z"))
qgraph(info_r$cor_matrix, graph="cor", labels=c("X", "Y", "Z"))
qgraph(info_s$cor_matrix, graph="cor", labels=c("X", "Y", "Z"))
dev.off()

#pdf("triplet_pcor.pdf", width = 12, height = 3)
par(mfrow = c(1, 3))
qgraph(info_i$pcor_matrix$estimate, graph="pcor", labels=c("X", "Y", "Z"))
qgraph(info_r$pcor_matrix$estimate, graph="pcor", labels=c("X", "Y", "Z"))
qgraph(info_s$pcor_matrix$estimate, graph="pcor", labels=c("X", "Y", "Z"))
dev.off()

#### Figure 8; tables 1, 2, C4: Higher-order information for n=3,4 ####

## IID (n=3) ##

options <- combn(names,3) # all combinations of triplets

OI_hb <- vector(mode="numeric",length=dim(options)[2]) 
OI_cov <- vector(mode="numeric",length=dim(options)[2])

for (i in 1:dim(options)[2]){
  print(dim(options)[2] - i)
  OI_hb[i] <- O_Info(data,options[,i],discrete=T,estimator="hist")
  OI_cov[i] <- O_Info(data,options[,i],estimator="cov")
}
option_string <- c()
for (i in 1:dim(options)[2]){
  option_string <- c(option_string,paste0(options[,i],collapse = " "))
}
Threewaystructure <- data.frame(option_string, OI_hb, OI_cov)
#write.table(Threewaystructure, file="briganti_3waystructure.txt", append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)

Threeway_sorted <- Threewaystructure[order(Threewaystructure$OI_hb),]

# figure 8
barplot(Threeway_sorted$OI_hb, border = rgb(0, 0, 1, alpha = 0.4), col=rgb(0, 0, 1, alpha = 0.3),ylab="O-information (nats)",cex.lab=1.3,cex.ax=1.5)
barplot(Threeway_sorted$OI_cov, border = rgb(1, 0, 0, alpha = 0.4), col=rgb(1, 0, 0, alpha = 0.3), add=TRUE, axes=F)
legend("topleft", legend = c("Histogram-based", "covariance-based"), fill = c(rgb(0, 0, 1, alpha = 0.5), rgb(1, 0, 0, alpha = 0.5)),cex=1.5)


## Calculate O-inf with bootstrap significance testing

source("C:/Users/nvsanten/OneDrive - UGent/Documenten/R Scripts/Information theory/bootstrap.R")


#for most synergistic and redundant triplets according to OI_hb (table 1)
syn_triple <- names[c(12,19,23)]

syn_hb <- bootstrap_O_information(data,var=syn_triple,n=3,nboot=1000,alpha=0.01,discrete = T,correction=T,estimator = "hist")

red_triple <- names[c(16,23,26)]

red_hb <- bootstrap_O_information(data,var=red_triple,n=3,nboot=1000,alpha=0.01,discrete= T, correction=T,estimator = "hist")

## IID (n=4) ##

options4 <- combn(names,4)
OI_hb <- vector(mode="numeric",length=dim(options4)[2]) 
OI_cov <- vector(mode="numeric",length=dim(options4)[2])
for (i in 1:dim(options4)[2]){
  print(dim(options4)[2] - i)
  OI_hb[i] <- O_Info(data,options4[,i],discrete=T, estimator= "hist")
  OI_cov[i] <- O_Info(data,options4[,i],estimator="cov")
}

option4_string <- c()
for (i in 1:dim(options4)[2]){
  option4_string <- c(option4_string,paste0(options4[,i],collapse = " "))
}
Fourwaystructure <- data.frame(option4_string, OI_hb, OI_cov)
#write.table(Fourwaystructure, file="briganti_4waystructure.txt", append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)

Fourway_sorted <- Fourwaystructure[order(Fourwaystructure$OI_hb),]

# figure 8
barplot(Fourway_sorted$OI_hb, border = rgb(0, 0, 1, alpha = 0.4), col=rgb(0, 0, 1, alpha = 0.3),ylab="O-information (nats)", ylim=c(min(Fourway_sorted[c("OI_hb","OI_cov")]),max(Fourway_sorted[c("OI_hb","OI_cov")])),cex.lab=1.3,cex.ax=1.5)
barplot(Fourway_sorted$OI_cov, border = rgb(1, 0, 0, alpha = 0.4), col=rgb(1, 0, 0, alpha = 0.3), add=TRUE, axes=F)

#bootstrap checking for highest synergistic and redundant quadruplet (table 1)

syn_quad <- names[c(2,5,15,16)]

syn_quad_hb <- bootstrap_O_information(data,var=syn_quad,n=4,nboot=1000,alpha=0.01,discrete = T,correction=T,estimator = "hist")

red_quad <- names[c(5,16,23,26)]

red_quad_cov <- bootstrap_O_information(data,var=red_quad,n=4,nboot=1000,alpha=0.01, correction=T,estimator = "cov")

## PID (triplets) ## (table 2)

# S=12FS_R, R={19PD_R,23FS} (highest synergy dominated triplet)
PID_MMI(data,"12FS_R",c("19PD_R","23FS"),estimator = "hist",discrete = T)

# S=16FS, R={23FS,26FS} (highest redundancy dominated triplet)
PID_MMI(data,"16FS",c("23FS","26FS"),estimator = "hist",discrete = T)

## PID (quadruplets) ## (table C4)

# S=2EC, R={5FS,15PT_R,16FS} (highest synergy dominated quadruplet)
PID_MMI(data,"2EC",c("5FS","15PT_R","16FS"),estimator = "hist",discrete=T)

# S=x1, R={x6,x8,x9} (highest redundancy dominated quadruplet)
PID_MMI(data,"5FS",c("16FS","23FS","26FS"),estimator = "hist",discrete = T)
