# Load necessary library
library(dplyr)

# Example probability distributions
P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**4  # Distribution Q

# Ensure the distributions are valid (non-negative and sum to 1)
P <- P / sum(P)
Q <- Q / sum(Q)







# Define the tail weight function (e.g., exponential weight for tail emphasis)
tail_weight_function <- function(index, lambda = 1, n) {
  # Normalize lambda based on the sample size
  lambda_adjusted <- lambda / length(index)
  # Tail weight function based on the index of the probability
  exp(-lambda_adjusted * abs(index - n))
}



# Calculate tail-weighted KL divergence for discrete distributions
tail_weighted_kl_divergence <- function(P, Q, tail_weight_function, lambda = 1) {
  # Ensure that the distributions are the same length
  if (length(P) != length(Q)) {
    stop("Distributions P and Q must be of the same length.")
  }
  
  n <- length(P)
  
  # Calculate the tail-weighted KL divergence
  kl_divergence <- sum(sapply(1:n, function(i) {
    p_x <- P[i]
    q_x <- Q[i]
    # Avoid division by zero and log(0)
    if (p_x == 0 || q_x == 0) {
      return(0)
    }
    weight <- tail_weight_function(i, lambda, n)
    return(weight * p_x * log(p_x / q_x))
  }))
  
  return(kl_divergence)
}

# Compute the tail-weighted Jensen-Shannon divergence
tail_weighted_js_divergence <- function(P, Q, tail_weight_function, lambda = 1) {
  # Compute the average distribution
  M <- (P + Q) / 2
  
  # Compute the tail-weighted KL divergences
  KL_PM <- tail_weighted_kl_divergence(P, M, tail_weight_function, lambda)
  KL_QM <- tail_weighted_kl_divergence(Q, M, tail_weight_function, lambda)
  
  # Average the KL divergences
  js_divergence <- (KL_PM + KL_QM) / 2
  
  return(js_divergence)
}

# Compute the tail-weighted Jensen-Shannon divergence
lambda <- 1  # Tail weight parameter
js_divergence <- tail_weighted_js_divergence(P, Q, tail_weight_function, lambda)
print(paste("Tail-Weighted Jensen-Shannon Divergence:", js_divergence))

