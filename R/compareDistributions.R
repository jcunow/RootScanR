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








#' Tail weight function (e.g., exponential weight for tail emphasis)
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "linear, "exponential", "sigmoid", "gompertz")
#' @param baseline.weight  minimal weight between 0-1
#' @keywords internal
#' @export
#' @return weights between 0-1 along the index
#'
#'
#'
#'
#'
tail_weight_function <- function(x, parameter = list(lambda = 0.2,x0=5),
                                 method = "sigmoid",baseline.weight = 0,inverse =FALSE) {

  x =x / sum(x,na.rm=TRUE)
  ## equal weighting along index
  if(method == "constant"){
    method = "exponential"
    parameter$lambda = 0
  }
  # exponential weighting function
  if(method == "exponential"){
    # sample size
    n = length(x)
    # Normalize lambda based to make it index scale invariant
    lambda_adjusted <- parameter$lambda / n
    # Tail weight function based on the index of the probability
    weights = exp(-lambda_adjusted * abs(x - n))*(1-baseline.weight) + baseline.weight
  }
  # sigmoid weighting function
  if(method == "sigmoid"){
    # k = steepness of curve, x0 = inflection point
    weights =  1 / (1 + exp(-parameter$lambda * (x - parameter$x0))) *(1-baseline.weight) + baseline.weight
  }
  # linerarily increasing weights
  if(method == "linear"){
    # Normalize linear weighting based on the sample size
    weights = (x / length(x))*(1-baseline.weight) + baseline.weight
  }
  # gompertz weighting
  if(method == "gompertz"){
    weights = (1/1+baseline.weight) * exp(-parameter$x0 * exp(-parameter$lambda * x)) *(1-baseline.weight) + baseline.weight
  }
  if(inverse == TRUE){
    weights = rev(weights)

  }

  weights = weights / sum(weights,na.rm=TRUE)
  return(weights)
}




#
#' Calculate tail-weighted KL divergence for discrete distributions
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "linear, "exponential", "sigmoid", "gompertz")
#'
#' @return KL divergence, not symmetrical - changing the input order will chenge the result
#' @export
#'
tail_weighted_kl_divergence <- function(P, Q, parameter = list(lambda = 0.2,x0=30),inverse=FALSE,method = "constant") {
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
    weight <- tail_weight_function(i,  parameter,method = method, inverse=inverse)
    return(weight * p_x * log(p_x / q_x))
  }))

  return(kl_divergence)
}





#' Compute the tail-weighted Jensen-Shannon divergence
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "linear, "exponential", "sigmoid", "gompertz")
#'
#' @return Jensen-Shannon Divergence - symmetric version of KL
#' @export
#'
#' @examples
#' P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
#' Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q
#'
#' # Ensure the distributions are valid (non-negative and sum to 1)
#' P <- P / sum(P)
#' Q <- Q / sum(Q)
#'
#' tail_weighted_js_divergence(P,Q,tail_weight_function(P,method="sigmoid),parameter = list(lambda = 0.2,x0=30))
tail_weighted_js_divergence <- function(P, Q,  parameter = list(lambda = 0.2,x0=30),method="constant",inverse=FALSE) {
  # Compute the average distribution
  M <- (P + Q) / 2

  # Compute the tail-weighted KL divergences
  KL_PM <- tail_weighted_kl_divergence(P, M, tail_weight_function(inverse = inverse,method=method), parameter)
  KL_QM <- tail_weighted_kl_divergence(Q, M, tail_weight_function(inverse = inverse,method=method), parameter)

  # Average the KL divergences
  js_divergence <- (KL_PM + KL_QM) / 2
  js_divergence <- round(js_divergence,4)

  return(js_divergence)
}


# Compute the tail-weighted Jensen-Shannon divergence
#parameter = list(lambda = 0.2,x0=30)  # Tail weight parameter
#js_divergence <- tail_weighted_js_divergence(P, Q, tail_weight_function, parameter)
#print(paste("Tail-Weighted Jensen-Shannon Divergence:", js_divergence))





#' A tailweighted Version of 1 dimensional Wasserstein distance betwwen two probability vectors
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "linear, "exponential", "sigmoid", "gompertz")
#'
#' @return
#' @export
#'
#' @examples
#' P <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)  # Distribution P
#' Q <- c(0.025,0.05,0.1,0.15, 0.2, 0.3,0.4, 0.5,0.3,0.1)**6  # Distribution Q
#'
#' # Ensure the distributions are valid (non-negative and sum to 1)
#' P <- P / sum(P)
#' Q <- Q / sum(Q)
#'
#' tail_weighted_wasserstein_distance(P,Q,inverse=F,method="constant",parameter = list(lambda = 0.2,x0=3))
tail_weighted_wasserstein_distance = function(Q,P,inverse=F,parameter = list(lambda = 0.2,x0=3),method= "constant"){
  # tail weighting
  weighted.P = tail_weight_function(P,inverse=inverse,parameter = parameter,method = method )
  weighted.Q = tail_weight_function(Q,inverse=inverse,parameter = parameter,method = method)
  # distance computation
  transport::wasserstein1d(a=weighted.P,b=weighted.Q)
}
