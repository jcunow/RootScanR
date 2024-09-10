#### Compare similarity between distributions


# # Example distributions
# P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
# Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q
#
# # Ensure the distributions are valid (non-negative and sum to 1)
# P <- P / sum(P)
# Q <- Q / sum(Q)
#
# plot(Q);points(P,col="red")








# Define the tail weight function (e.g., exponential weight for tail emphasis)
# lambda = shape parameter (0 = constant weighting), x0 = curve offset (= inflexion point )
tail_weight_function <- function(x, parameter = list(lambda = 0.2,x0=30), method = "sigmoid",baseline.weight = 0.1) {

  ## equal weighting along index
  if(method == "constant"){
    method = "exp"
    lambda = 0
  }
  # exponential weighting function
  if(method == "exp"){
    # sample size
    n = max(x,na.rm=TRUE)
    # Normalize lambda based to make it index scale invariant
    lambda_adjusted <- parameter$lambda / n
    # Tail weight function based on the index of the probability
    val = exp(-lambda_adjusted * abs(x - n))*(1-baseline.weight) + baseline.weight
  }
  # sigmoid weighting function
  if(method == "sigmoid"){
    # k = steepness of curve, x0 = inflection point
   val =  1 / (1 + exp(-parameter$lambda * (x - parameter$x0))) *(1-baseline.weight) + baseline.weight
  }
  # linerarily increasing weights
  if(method == "linear"){
    # Normalize linear weighting based on the sample size
    val = (x / length(x))*(1-baseline.weight) + baseline.weight
  }
  # gompertz weighting
  if(method == "gompertz"){
      val = (1/1+baseline.weight) * exp(-parameter$x0 * exp(-parameter$lambda * x)) *(1-baseline.weight) + baseline.weight
  }

  return(val)
}




# Calculate tail-weighted KL divergence for discrete distributions
tail_weighted_kl_divergence <- function(P, Q, tail_weight_function, parameter = list(lambda = 0.2,x0=30)) {
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
    weight <- tail_weight_function(i,  parameter)
    return(weight * p_x * log(p_x / q_x))
  }))

  return(kl_divergence)
}


# Compute the tail-weighted Jensen-Shannon divergence
tail_weighted_js_divergence <- function(P, Q, tail_weight_function, parameter = list(lambda = 0.2,x0=30)) {
  # Compute the average distribution
  M <- (P + Q) / 2

  # Compute the tail-weighted KL divergences
  KL_PM <- tail_weighted_kl_divergence(P, M, tail_weight_function, parameter)
  KL_QM <- tail_weighted_kl_divergence(Q, M, tail_weight_function, parameter)

  # Average the KL divergences
  js_divergence <- (KL_PM + KL_QM) / 2
  js_divergence <- round(js_divergence,4)

  return(js_divergence)
}


# Compute the tail-weighted Jensen-Shannon divergence
parameter = list(lambda = 0.2,x0=30)  # Tail weight parameter
js_divergence <- tail_weighted_js_divergence(P, Q, tail_weight_function, parameter)
print(paste("Tail-Weighted Jensen-Shannon Divergence:", js_divergence))

