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
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point? )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param baseline.weight  minimal weight between 0-1
#'
#' @keywords internal
#' @export
#' @return weights between 0-1 along the index
#'
#'
#'
#'
#'
tail_weight_function <- function(index = NULL, parameter = list(lambda = 0.2,x0=5), index.spacing = "equal",
                                 method = "sigmoid",baseline.weight = 0,inverse =FALSE) {

  if( index.spacing == "equal"){
    # sample size defines index. equally spaced
    n = length(index)
    index = 1:n
  }else{
    # specifying depth allows unequal spacing
    n = max(index,na.rm=T)
    index = index
  }


  ## equal weighting along index
  if(method == "constant"){
    weights = 1
  }
  # asymptotic
  if(method == "asympotic"){
    weights = ((-exp( parameter$lambda / n * abs((index) - n)) + max(exp( parameter$lambda / n * abs((index) - n)))) /
                 max(exp( parameter$lambda / n * abs((index) - n)))) * (1-baseline.weight) + baseline.weight

  }
  # exponential weighting function
  if(method == "exponential"){

    # Tail weight function based on the index of the probability
    weights = exp(- parameter$lambda / n * abs((index) - n)) * (1-baseline.weight) + baseline.weight
  }
  # sigmoid weighting function
  if(method == "sigmoid"){
    # k = steepness of curve, x0 = inflection point
    weights =  1 / (1 + exp(-parameter$lambda * (index - parameter$x0))) *(1-baseline.weight) + baseline.weight
  }
  # linerarily increasing weights
  if(method == "linear"){
    # Normalize linear weighting based on the sample size
    weights = ((index) / n)  *(1-baseline.weight) + baseline.weight
  }
  # gompertz weighting
  if(method == "gompertz"){
    weights =  exp(-parameter$x0 * exp(-parameter$lambda * index)) *(1-baseline.weight) + baseline.weight
  }

  # step weighting
  if(method == "step"){
    weights = index
    weights[index <= parameter$x0] <- 0
    weights[index > parameter$x0] <- 1
    weights = weights * (1-baseline.weight) + baseline.weight
  }

  if(inverse == TRUE){
    weights = rev(weights)
  }

  out = (weights) / sum((weights))
  return(out)
}




#
#' Calculate tail-weighted KL divergence for discrete distributions
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param alignPQ if TRUE, index end values will be cut off in case of unequal length of P & Q so that length of P & Q is equal
#' @param cut if FALSE, 0 will be added to the shorter vector. If TRUE, the longer vector will be shortened at the end.
#'
#' @details
#' Kullback-Leibler Divergence
#'
#'
#' @return KL divergence, not symmetrical - changing the input order will change the result
#' @export
#'
tail_weighted_kl_divergence <- function(P, Q, index= 1:min(c(length(Q),length(P))),index.spacing = "equal",  parameter = list(lambda = 0.2,x0=30),cut = FALSE,inverse=FALSE,method = "constant", alignPQ = TRUE) {
  # Ensure that the distributions are the same length
  if (length(P) != length(Q) & alignPQ != TRUE) {
    stop("Distributions P and Q must be of the same length.")
  }

  n = min(c(length(Q),length(P)))
  m = max(c(length(Q),length(P)))
  if(alignPQ == TRUE){
    if(cut == TRUE){
      P = P[1:n]
      Q = Q[1:n]
    }
    if(cut == FALSE){
      vl1 = abs(length(P) - m)
      vl2 = abs(length(Q) - m)
      P = c(P[1:m-vl1],rep(0,vl1))
      Q = c(Q[1:m-vl2],rep(0,vl2))
    }
  }

  weight <- tail_weight_function(index = index, index.spacing = index.spacing,  parameter =  parameter,method = method, inverse=inverse)

  P = (P * weight) / sum(P * weight)
  Q = (Q * weight) / sum(Q * weight)

  # Calculate the tail-weighted KL divergence
  kl_divergence <- sum(sapply(index, function(i) {
    p_x <- P[i]
    q_x <- Q[i]
    # Avoid division by zero and log(0)
    if (p_x <= 0 || q_x <= 0) {
      return(0)
    }

    return( p_x * log(p_x / q_x))
  }))

  return(kl_divergence)
}





#' Compute the tail-weighted Jensen-Shannon divergence
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "asymptotic", "linear, "exponential", "sigmoid", "gompertz","step")
#' @param alignPQ if TRUE, index end values will be cut off in case of unequal length of P & Q so that length of P & Q is equal
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
#' tail_weighted_js_divergence(P,Q,parameter = list(lambda = 0.2,x0=30))
tail_weighted_js_divergence <- function(P, Q,  parameter = list(lambda = 0.2,x0=30),method="constant",inverse=FALSE, alignPQ = TRUE,index= 1:min(c(length(Q),length(P))),index.spacing = "equal") {

  # Ensure that the distributions are the same length
  if (length(P) != length(Q)) {
    stop("Distributions P and Q must be of the same length.")
  }
  # Compute the average distribution
  M <- (P + Q) / 2

  # Compute the tail-weighted KL divergences
  KL_PM <- tail_weighted_kl_divergence(P, M, index = index, index.spacing = index.spacing,  parameter =  parameter,method = method, inverse=inverse, alignPQ = alignPQ)
  KL_QM <- tail_weighted_kl_divergence(Q, M, index = index, index.spacing = index.spacing,  parameter =  parameter,method = method, inverse=inverse, alignPQ = alignPQ)

  # Average the KL divergences
  js_divergence <- (KL_PM + KL_QM) / 2
  js_divergence <- round(js_divergence,4)

  return(js_divergence)
}






#' A tailweighted Version of 1 dimensional Wasserstein distance betwwen two probability vectors
#'
#' @param Q probability vector 1
#' @param P probability vector 2
#' @param inverse changes from right tail to left tail if TRUE
#' @param parameter list with lambda -> shape parameter (0 = constant weighting) & x0 -> curve offset (= inflexion point )
#' @param method weighting function along index. Available options are: c("constant", "linear, "exponential", "sigmoid", "gompertz")
#' @param index a positive numeric vector containing probability spacing e.g., depth
#' @param index.spacing whether index intervals are equally distant i.e., c(1,2,3,4....n), if "equal" than index is c(1,n)
#' @param baseline.weight  minimal weight between 0-1
#'
#' @return wasserstein metric
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
#' tail_weighted_wasserstein_distance(P,Q,
#'   inverse=F,method="constant",parameter = list(lambda = 0.2,x0=3))
tail_weighted_wasserstein_distance = function(Q,P,inverse=F,parameter = list(lambda = 0.2,x0=10),
                                              method= "step",baseline.weight = 0,
                                              index = 1:min(c(length(Q),length(P))), index.spacing = "equal"){
  # tail weighting
  weights= tail_weight_function(index = index, parameter, index.spacing = index.spacing, method = method,baseline.weight, inverse =FALSE)

  P = (P*weights) / sum(P*weights)
  Q = (Q*weights) / sum(Q*weights)

  # distance computation
  transport::wasserstein1d(a=P,b=Q)
}



#### root index
#' Root Weight Depth Index
#'
#' @param w weight, typically a vector containing depths
#' @param roots root cover
#'
#' @return a value between 0-1
#' @export
#'
#' @examples
#' w = seq(5,25,5)
#' roots = c(0,10,7,3,1)
#' RWDI(w,roots)
RWDI = function(w,roots){
  rwdi = sum((w+0.000000000000001)*roots)/sum(roots)
  return(rwdi)
}


#' Root Penetration Index
#'
#' @param roots root cover
#' @param w weight, typically a vector containing depths
#'
#' @return a value between 0-1
#' @export
#'
#' @examples
#' w = seq(5,25,5)
#' roots = c(0,10,7,3,1)
#' RPI(w,roots)
RPI = function(roots,w){

  rpi = 1-2*(sum(roots/sum(roots) * ((w+0.0000000000000000001)/(sum(w)))))
  return(rpi)
}

