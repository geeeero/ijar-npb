test <- function (at.times, survival.signature, test.data, nLower = 2,
                  nUpper = 2, yLower = 0.5, yUpper = 0.5)
{
  OPTl <- OPTu <- 0
  THYl <- THYu <- 0
  K <- ncol(survival.signature) - 1
  if (any(at.times < 0)) {
    stop("all at.times must be non-negative")
  }
  if (!any(prob.col <- (names(survival.signature) == "Probability"))) {
    stop("survival signature must contain a variable named 'Probability'")
  }
  if (!is.list(test.data) || length(test.data) != K) {
    stop("test.data must be a named list containing the same number of components as specified by the survival signature")
  }
  if (!all(sort(names(test.data)) == sort(names(survival.signature)[!prob.col]))) {
    stop("component names in survival.signature and test.data must match exactly")
  }
  if (is.vector(nLower) && is.vector(nUpper) && is.vector(yLower) &&
      is.vector(yUpper)) {
    if (!all.equal(length(nLower), length(nUpper), length(yLower),
                   length(yUpper))) {
      stop("nLower, nUpper, yLower and yUpper prior parameter vectors must be the same length")
    }
    if (length(nLower) != 1 && length(nLower) != length(at.times)) {
      stop("prior parameter vectors, nLower, nUpper, yLower and yUpper, must be either length 1 (for time homogeneous prior) or of the same length as at.times argument (where the prior parameters at time at.times[i] are now nLower[i], nUpper[i], yLower[i] and yUpper[i])")
    }
    nLower <- as.data.frame(matrix(rep(nLower, K), nrow = length(at.times),
                                   ncol = K, byrow = FALSE))
    names(nLower) <- names(survival.signature[!prob.col])
    nUpper <- as.data.frame(matrix(rep(nUpper, K), nrow = length(at.times),
                                   ncol = K, byrow = FALSE))
    names(nUpper) <- names(survival.signature[!prob.col])
    yLower <- as.data.frame(matrix(rep(yLower, K), nrow = length(at.times),
                                   ncol = K, byrow = FALSE))
    names(yLower) <- names(survival.signature[!prob.col])
    yUpper <- as.data.frame(matrix(rep(yUpper, K), nrow = length(at.times),
                                   ncol = K, byrow = FALSE))
    names(yUpper) <- names(survival.signature[!prob.col])
  }
  else if (is.data.frame(nLower) && is.data.frame(nUpper) &&
           is.data.frame(yLower) && is.data.frame(yUpper)) {
    if (ncol(nLower) != K || ncol(nUpper) != K || ncol(yLower) !=
        K || ncol(yUpper) != K) {
      stop("nLower, nUpper, yLower and yUpper must have priors for the same number of components as specified by the survival signature")
    }
    if (!all.equal(nrow(nLower), nrow(nUpper), nrow(yLower),
                   nrow(yUpper))) {
      stop("nLower, nUpper, yLower and yUpper must have matching size (they differ in number of rows)")
    }
    if (!all(sort(names(nLower)) == sort(names(survival.signature)[!prob.col])) ||
        !all(sort(names(nUpper)) == sort(names(survival.signature)[!prob.col])) ||
        !all(sort(names(yLower)) == sort(names(survival.signature)[!prob.col])) ||
        !all(sort(names(yUpper)) == sort(names(survival.signature)[!prob.col]))) {
      stop("component names in survival.signature and nLower, nUpper, yLower and yUpper prior lists must match exactly")
    }
    if (nrow(nLower) != 1 && nrow(nLower) != length(at.times)) {
      stop("prior parameter vectors, nLower, nUpper, yLower and yUpper, must be either length 1 (for time homogeneous prior) or of the same length as at.times argument (where the prior parameters at time at.times[i] are now nLower[i,j], nUpper[i,j], yLower[i,j] and yUpper[i,j]).")
    }
    nLower <- nLower[, names(survival.signature)[!prob.col]]
    nUpper <- nUpper[, names(survival.signature)[!prob.col]]
    yLower <- yLower[, names(survival.signature)[!prob.col]]
    yUpper <- yUpper[, names(survival.signature)[!prob.col]]
    if (nrow(nLower) == 1) {
      nLower <- as.data.frame(matrix(nLower, nrow = length(at.times),
                                     ncol = K, byrow = TRUE))
      names(nLower) <- names(survival.signature)[!prob.col]
      nUpper <- as.data.frame(matrix(nUpper, nrow = length(at.times),
                                     ncol = K, byrow = TRUE))
      names(nUpper) <- names(survival.signature)[!prob.col]
      yLower <- as.data.frame(matrix(yLower, nrow = length(at.times),
                                     ncol = K, byrow = TRUE))
      names(yLower) <- names(survival.signature)[!prob.col]
      yUpper <- as.data.frame(matrix(yUpper, nrow = length(at.times),
                                     ncol = K, byrow = TRUE))
      names(yUpper) <- names(survival.signature)[!prob.col]
    }
  }
  else {
    stop("nLower, nUpper, yLower and yUpper arguments must be either a vector or data frame and must match in type")
  }
  cores <- detectCores()
  cores <- ifelse(is.na(cores), 1, cores)
  pLower <- simplify2array(lapply(order(at.times), function(i,
                                                            at.times, nLower, nUpper, yLower, yUpper, sig, prob.col,
                                                            test.data, m, N, K) {
    t <- at.times[i]
    nLower <- unlist(nLower[i, ])
    nUpper <- unlist(nUpper[i, ])
    yLower <- unlist(yLower[i, ])
    s <- sapply(test.data, function(t_i, t) {
      sum(t_i > t)
    }, t = t)
    n <- split(cbind(nLower, nUpper), 1:K)
    if (all(N + m - 1 > 0)) {
      yl <- s/(N + m - 1)
      yu <- (s + m - 1)/(N + m - 1)
      for (k in 1:K) {
        if (yLower[k] < yl[k]) {
          n[[k]] <- nUpper[k]
        }
        else if (yLower[k] > yu[k]) {
          n[[k]] <- nLower[k]
        }
        else if ((lgamma(m[k] + nUpper[k] * (1 - yLower[k]) +
                         N[k] - s[k]) + lgamma(nLower[k] * (1 - yLower[k]) +
                                               N[k] - s[k]) + lgamma(nUpper[k] + N[k]) + lgamma(m[k] +
                                                                                                nLower[k] + N[k]) - lgamma(m[k] + nLower[k] *
                                                                                                                           (1 - yLower[k]) + N[k] - s[k]) - lgamma(nUpper[k] *
                                                                                                                                                                   (1 - yLower[k]) + N[k] - s[k]) - lgamma(nLower[k] +
                                                                                                                                                                                                           N[k]) - lgamma(m[k] + nUpper[k] + N[k]) >=
                  0) && (lgamma(m[k] + nUpper[k] * yLower[k] +
                                s[k]) + lgamma(nLower[k] * yLower[k] + s[k]) +
                         lgamma(nUpper[k] + N[k]) + lgamma(m[k] + nLower[k] +
                                                           N[k]) - lgamma(m[k] + nLower[k] * yLower[k] +
                                                                          s[k]) + lgamma(nUpper[k] * yLower[k] + s[k]) +
                         lgamma(nLower[k] + N[k]) + lgamma(m[k] + nUpper[k] +
                                                           N[k]) <= 0)) {
          n[[k]] <- nUpper[k]
        }
        else if ((lgamma(m[k] + nUpper[k] * (1 - yLower[k]) +
                         N[k] - s[k]) + lgamma(nLower[k] * (1 - yLower[k]) +
                                               N[k] - s[k]) + lgamma(nUpper[k] + N[k]) + lgamma(m[k] +
                                                                                                nLower[k] + N[k]) - lgamma(m[k] + nLower[k] *
                                                                                                                           (1 - yLower[k]) + N[k] - s[k]) - lgamma(nUpper[k] *
                                                                                                                                                                   (1 - yLower[k]) + N[k] - s[k]) - lgamma(nLower[k] +
                                                                                                                                                                                                           N[k]) - lgamma(m[k] + nUpper[k] + N[k]) <=
                  0) && (lgamma(m[k] + nUpper[k] * yLower[k] +
                                s[k]) + lgamma(nLower[k] * yLower[k] + s[k]) +
                         lgamma(nUpper[k] + N[k]) + lgamma(m[k] + nLower[k] +
                                                           N[k]) - lgamma(m[k] + nLower[k] * yLower[k] +
                                                                          s[k]) - lgamma(nUpper[k] * yLower[k] + s[k]) -
                         lgamma(nLower[k] + N[k]) - lgamma(m[k] + nUpper[k] +
                                                           N[k]) >= 0)) {
          n[[k]] <- nLower[k]
        }
      }
    }
    if (max(sapply(n, length)) > 1) {
      OPTl <<- OPTl+sum(sapply(n, length)>1)
      THYl <<- THYl+sum(sapply(n, length)==1)
      nGrid <- expand.grid(split(cbind(nLower, nUpper),
                                 1:K))
      res <- apply(nGrid, 1, function(n, sig, prob.col,
                                      s, m, N, K, yLower) {
        sum(apply(sig, 1, function(sigvec, prob.col,
                                   s, m, N, n, yLower) {
          l <- sigvec[!prob.col]
          sig <- sigvec[prob.col]
          sig * prod(choose(m, l) * beta(l + n * yLower +
                                           s, m - l + n * (1 - yLower) + N - s)/beta(n *
                                                                                       yLower + s, n * (1 - yLower) + N - s))
        }, prob.col = prob.col, s = s, m = m, N = N,
        n = n, yLower = yLower))
      }, sig = sig, prob.col = prob.col, s = s, m = m,
      N = N, K = K, yLower = yLower)
      res <- min(res)
    }
    else {
      THYl <<- THYl+K
      res <- apply(matrix(unlist(n, use.names = FALSE),
                          nrow = 1), 1, function(n, sig, prob.col, s, m,
                                                 N, K, yLower) {
                            sum(apply(sig, 1, function(sigvec, prob.col,
                                                       s, m, N, n, yLower) {
                              l <- sigvec[!prob.col]
                              sig <- sigvec[prob.col]
                              sig * prod(choose(m, l) * beta(l + n * yLower +
                                                               s, m - l + n * (1 - yLower) + N - s)/beta(n *
                                                                                                           yLower + s, n * (1 - yLower) + N - s))
                            }, prob.col = prob.col, s = s, m = m, N = N,
                            n = n, yLower = yLower))
                          }, sig = sig, prob.col = prob.col, s = s, m = m,
                   N = N, K = K, yLower = yLower)
    }
    res
  }, at.times = at.times, nLower = nLower, nUpper = nUpper,
  yLower = yLower, yUpper = yUpper, sig = survival.signature,
  prob.col = prob.col, test.data = test.data, m = apply(survival.signature[,
                                                                           -length(survival.signature), drop = FALSE], 2, max),
  N = sapply(test.data, length), K = K))[rank(at.times)]
  pUpper <- simplify2array(lapply(order(at.times), function(i,
                                                            at.times, nLower, nUpper, yLower, yUpper, sig, prob.col,
                                                            test.data, m, N, K) {
    t <- at.times[i]
    nLower <- unlist(nLower[i, ])
    nUpper <- unlist(nUpper[i, ])
    yUpper <- unlist(yUpper[i, ])
    s <- sapply(test.data, function(t_i, t) {
      sum(t_i > t)
    }, t = t)
    n <- split(cbind(nLower, nUpper), 1:K)
    if (all(N + m - 1 > 0)) {
      yl <- s/(N + m - 1)
      yu <- (s + m - 1)/(N + m - 1)
      for (k in 1:K) {
        if (yUpper[k] < yl[k]) {
          n[[k]] <- nLower[k]
        }
        else if (yUpper[k] > yu[k]) {
          n[[k]] <- nUpper[k]
        }
        else if ((lgamma(m[k] + nUpper[k] * (1 - yUpper[k]) +
                         N[k] - s[k]) + lgamma(nLower[k] * (1 - yUpper[k]) +
                                               N[k] - s[k]) + lgamma(nUpper[k] + N[k]) + lgamma(m[k] +
                                                                                                nLower[k] + N[k]) - lgamma(m[k] + nLower[k] *
                                                                                                                           (1 - yUpper[k]) + N[k] - s[k]) - lgamma(nUpper[k] *
                                                                                                                                                                   (1 - yUpper[k]) + N[k] - s[k]) - lgamma(nLower[k] +
                                                                                                                                                                                                           N[k]) - lgamma(m[k] + nUpper[k] + N[k]) >=
                  0) && (lgamma(m[k] + nUpper[k] * yUpper[k] +
                                s[k]) + lgamma(nLower[k] * yUpper[k] + s[k]) +
                         lgamma(nUpper[k] + N[k]) + lgamma(m[k] + nLower[k] +
                                                           N[k]) - lgamma(m[k] + nLower[k] * yUpper[k] +
                                                                          s[k]) + lgamma(nUpper[k] * yUpper[k] + s[k]) +
                         lgamma(nLower[k] + N[k]) + lgamma(m[k] + nUpper[k] +
                                                           N[k]) <= 0)) {
          n[[k]] <- nLower[k]
        }
        else if ((lgamma(m[k] + nUpper[k] * (1 - yUpper[k]) +
                         N[k] - s[k]) + lgamma(nLower[k] * (1 - yUpper[k]) +
                                               N[k] - s[k]) + lgamma(nUpper[k] + N[k]) + lgamma(m[k] +
                                                                                                nLower[k] + N[k]) - lgamma(m[k] + nLower[k] *
                                                                                                                           (1 - yUpper[k]) + N[k] - s[k]) - lgamma(nUpper[k] *
                                                                                                                                                                   (1 - yUpper[k]) + N[k] - s[k]) - lgamma(nLower[k] +
                                                                                                                                                                                                           N[k]) - lgamma(m[k] + nUpper[k] + N[k]) <=
                  0) && (lgamma(m[k] + nUpper[k] * yUpper[k] +
                                s[k]) + lgamma(nLower[k] * yUpper[k] + s[k]) +
                         lgamma(nUpper[k] + N[k]) + lgamma(m[k] + nLower[k] +
                                                           N[k]) - lgamma(m[k] + nLower[k] * yUpper[k] +
                                                                          s[k]) - lgamma(nUpper[k] * yUpper[k] + s[k]) -
                         lgamma(nLower[k] + N[k]) - lgamma(m[k] + nUpper[k] +
                                                           N[k]) >= 0)) {
          n[[k]] <- nUpper[k]
        }
      }
    }
    if (max(sapply(n, length)) > 1) {
      OPTu <<- OPTu+sum(sapply(n, length)>1)
      THYu <<- THYu+sum(sapply(n, length)==1)
      nGrid <- expand.grid(split(cbind(nLower, nUpper),
                                 1:K))
      res <- apply(nGrid, 1, function(n, sig, prob.col,
                                      s, m, N, K, yUpper) {
        sum(apply(sig, 1, function(sigvec, prob.col,
                                   s, m, N, n, yUpper) {
          l <- sigvec[!prob.col]
          sig <- sigvec[prob.col]
          sig * prod(choose(m, l) * beta(l + n * yUpper +
                                           s, m - l + n * (1 - yUpper) + N - s)/beta(n *
                                                                                       yUpper + s, n * (1 - yUpper) + N - s))
        }, prob.col = prob.col, s = s, m = m, N = N,
        n = n, yUpper = yUpper))
      }, sig = sig, prob.col = prob.col, s = s, m = m,
      N = N, K = K, yUpper = yUpper)
      res <- max(res)
    }
    else {
      THYu <<- THYu+K
      res <- apply(matrix(unlist(n, use.names = FALSE),
                          nrow = 1), 1, function(n, sig, prob.col, s, m,
                                                 N, K, yUpper) {
                            sum(apply(sig, 1, function(sigvec, prob.col,
                                                       s, m, N, n, yUpper) {
                              l <- sigvec[!prob.col]
                              sig <- sigvec[prob.col]
                              sig * prod(choose(m, l) * beta(l + n * yUpper +
                                                               s, m - l + n * (1 - yUpper) + N - s)/beta(n *
                                                                                                           yUpper + s, n * (1 - yUpper) + N - s))
                            }, prob.col = prob.col, s = s, m = m, N = N,
                            n = n, yUpper = yUpper))
                          }, sig = sig, prob.col = prob.col, s = s, m = m,
                   N = N, K = K, yUpper = yUpper)
    }
    res
  }, at.times = at.times, nLower = nLower, nUpper = nUpper,
  yLower = yLower, yUpper = yUpper, sig = survival.signature,
  prob.col = prob.col, test.data = test.data, m = apply(survival.signature[,
                                                                           -length(survival.signature), drop = FALSE], 2, max),

  N = sapply(test.data, length), K = K))[rank(at.times)]
  print("OPT")
  print(OPTl)
  print(OPTu)
  print("THY")
  print(THYl)
  print(THYu)
  list(lower = pLower, upper = pUpper)
}
