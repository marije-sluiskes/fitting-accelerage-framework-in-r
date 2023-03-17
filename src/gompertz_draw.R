
################################################################
# This script contains functions to draw event times from a Gompertz distribution
# parameterized either as an AFT model or as a PH model.
# h_0(t) = a * exp(b * t) in the PH case
# h_0(t) = tau/sigma * exp(t / sigma) in the AFT case (sigma = 1/b, tau = a/b)
################################################################

# draw from Gompertz AFT 
rgompertz_aft <- function(n, sigma, tau, linpred){
  u <- runif(n)
  return( (sigma * log( -log(u)/tau + 1) ) / exp(linpred) )
}

# draw from condtional Gompertz AFT (using rejection sampling)
rcgompertz_aft <- function(n, sigma, tau, linpred, s){
  
  all_vals <- vector(length = n)
  count <- 1
  while(count < n + 1){
    
    new_val <-  (sigma * log( -log(runif(1))/tau + 1) ) / exp(linpred)
    
    if(new_val >= s){
      all_vals[count] <- new_val
      count <- count + 1
    }
    
  }
  return(all_vals)
}

# draw from Gompertz PH
rgompertz_ph <- function(n, a, b, linpred){
  u <- runif(n)
  a <- a * exp(linpred)
  return( log( 1 - log(u) * b / a) / b )
}

# draw from conditional Gompertz PH
rcgompertz_ph <- function(n, a, b, linpred, s){
  u <- runif(n)
  a <- a * exp(linpred)
  return( log( exp(b * s) - log(u) * b / a) / b )
}

# get Gompertz baseline survival probability 
gomp_baseline_surv <- function(t, a, b){
  exp( (a/b) * (1 - exp(b*t)))
}

# get inverse of Gompertz baseline survival probability
inverse_gomp_baseline_surv <- function(s, a, b){
  (1/b) * log (1 - (b/a) * log(s))
}

# returns gomp_baseline_surv ^ exp(linpred)
gomp_PH_surv <- function(t, a, b, linpred){
  gomp_baseline_surv(t, a, b) ^ exp(linpred)
}

