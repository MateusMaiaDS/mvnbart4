libs=c("truncnorm","msm","inline","Rcpp","RcppArmadillo","rbenchmark")
if( sum(!(libs %in% .packages(all.available = TRUE)))>0){ install.packages(libs[!(libs %in% .packages(all.available = TRUE))])}
for(i in 1:length(libs)) {library(libs[i],character.only = TRUE,quietly=TRUE)}


#needed for openMP parallel
# Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
# Sys.setenv("PKG_LIBS"="-fopenmp")

#no of cores for openMP version
cores = 1

#surce code from same dir
Rcpp::sourceCpp('src/tn_normal_sampler.cpp')


#sample size
nn=1000


bb= 100
aa=-100
microbenchmark::microbenchmark( rnorm_trunc(upper = 0,lower = -Inf , mu = 0, sigma = 1)   )
rnorm_trunc(mu = 0,sigma = 1,lower = -Inf,upper = 10)
