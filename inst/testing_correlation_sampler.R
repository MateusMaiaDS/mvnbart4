rm(list=ls())
# Loading the mian function
Rcpp::sourceCpp("src/demo.cpp")

# function to create correlation matrix ouf of correlation coefficients
makeSigma <- function(sigma, d){
     Sigma <- diag(d)

     indices <- expand.grid(1:d, 1:d)
     indices <- indices[indices[,1] < indices[2],]
     indices <- indices[order(indices[,1],indices[,2]),]

     for (i in 1:d){
          for (j in 1:d){
               Sigma[i,j] <- ifelse(
                    i == j,
                    1,
                    ifelse(
                         i < j,
                         sigma[which(indices[,1] == i & indices[,2] == j)],
                         sigma[which(indices[,1] == j & indices[,2] == i)]
                    )
               )
          }
     }
     return(Sigma)
}

# Simluating Residuals
n <- 100
d <- 2

sigma_true <- c(-0.8) # must be of length (d^2 - d)/2
Sigma_true <- makeSigma(sigma_true, d)
det(Sigma_true) # must be > 0
Sigma_true_chol <- t(chol(Sigma_true))
resid <- matrix(NA, nrow = n, ncol = d)

for (i in 1:n){
     resid[i,] <- Sigma_true_chol %*% rnorm(d)
}


# Setting the initial values
sigma0 <- rep(0.0, (d^2 - d)/2)
df_ <- 5

#Just to auxiliar I gonna define y_mat as the residuals and y_hat as the zero matrix
y_mat_ <- resid
y_hat_ <- matrix(0,nrow = nrow(y_mat_),ncol = ncol(y_mat_))
n_mcmc <- 2000

sigma_post <- sigma_sampler(nmcmc = n_mcmc,
              d = d,
              sigma_0 = sigma0,
              y_mat = y_mat_,
              y_hat = y_hat_,df = df_)

# Plottin for the 2d case
if(d==2){
     plot(c(sigma_post),type = "l")
     abline(h = sigma_true, lty = 2, col = "blue", lwd = 2)
}
