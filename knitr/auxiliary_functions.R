

library("Bmisc")
library("correlateR")

# Log determinant function
logdet <- function(x, ...) {
  z <- determinant(x, logarithm = TRUE, ...)
  stopifnot(z$sign == 1)
  return(z$modulus)
}


# Log likelihood
loglik <- function(Psi, nu, S, ns) {  # Log-likelihood
  stopifnot(length(nu) == 1)

  k <- length(S)
  p <- nrow(S[[1]])

  cs <- (nu + ns)/2
  logdetPsi <- logdet(Psi)
  logdetPsiPlusS <- sapply(S, function(s) logdet(Psi + s))

  const <- sum(((ns * p)/2) * log(2))
  t1 <- k*nu/2 * logdetPsi
  t2 <- sum(lgammap(cs, p = p))
  t3 <- -sum(cs * logdetPsiPlusS)
  t4 <- -k*lgammap(nu/2, p = p)

  return(const + t1 + t2 + t3 + t4)
}

dloglik <- function(Psi, nu, S, ns) {  # Derivative
  k <- length(S)
  p <- nrow(S[[1]])
  t1 <-  k/2 * logdet(Psi)
  t2 <-  1/2 * sum(sapply(ns, function(ni) digammap((nu + ni)/2, p = p)))
  t3 <- -1/2 * sum(sapply(S, function(s) logdet(Psi + s)))
  t4 <- -k/2 * digammap(nu/2, p = p)

  return(t1 + t2 + t3 + t4)
}

get.nu <- function(Psi, nu, S, ns, interval) {

  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    loglik(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)

  # # Find extrema as root
  # diff_loglik_nu <- function(nu) { # Derivative as a function of nu
  #   dloglik(Psi, nu, S, ns)
  # }
  # res2 <- uniroot(f = diff_dloglik_nu, interval = interval)$root
  #
  # if (abs((res-res2)/res) >= 1e-4) {
  #   stop("methods do not agree on maxima")
  # }
  # return((res + res2)/2)

}

EMstep <- function(Psi, nu, S, ns) {
  k <- length(S)
  p <- nrow(S[[1]])
  co <- 1/(k*nu)
  t <- lapply(seq_along(S), function(i) (co*(ns[i] + nu))*solve(Psi + S[[i]]))
  Psi_new <- solve(Reduce("+", t))  # Sum the matrices in t
  return(Psi_new)
}

Momentstep <- function(nu, S, ns) {
  k <- length(ns)
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) S[[i]]/ns[i]))/k
  p <- nrow(Psi)
  fac <- nu - p - 1
  return(fac*Psi)
}

MLEstep <- function(nu, S, ns) {
  n.tot <- sum(ns)
  fac <- nu + ns
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) fac[i]*S[[i]]))/n.tot
  return(Psi)
}

fit.grem.EM <- function(Psi.init, nu.init, S, ns,
                        max.ite = 1000, eps = 1e-3,
                        verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 1e-5, 1e6)

  Psi.old <- Psi.init
  nu.old <- nu.init

  for (i in seq_len(max.ite)) {
    ll.old <- loglik(Psi.old, nu.old, S, ns)

    Psi.new <- EMstep(Psi.old, nu.old, S, ns)
    nu.new  <- get.nu(Psi.new, nu.old, S, ns, interval)
    ll.new  <- loglik(Psi.new, nu.new, S, ns)

    stopifnot(ll.new > ll.old)
    if (ll.new - ll.old < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old <- nu.new
    }
    if (verbose) {
      cat("ite =", i, ":", "ll.new - ll.old =", ll.new - ll.old, "\n");
      flush.console()
    }
  }

  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")

  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}


fit.grem.MLE <- function(nu.init, S, ns,
                         max.ite = 1000, eps = 1e-3,
                         verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 1e-5, 1e6)

  nu.old <- nu.init
  Psi.old <- MLEstep(nu = nu.old, S = S, ns = ns)

  for (i in seq_len(max.ite)) {
    ll.old <- loglik(Psi.old, nu.old, S, ns)

    nu.new  <- get.nu(Psi.old, nu.old, S, ns, interval)
    Psi.new <- MLEstep(nu.new, S, ns)

    ll.new  <- loglik(Psi.new, nu.new, S, ns)

    stopifnot(ll.new > ll.old)
    if (ll.new - ll.old < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old  <- nu.new
    }

    if (verbose) {
      cat("ite =", i, ":", "ll.new - ll.old =", ll.new - ll.old, "\n");
      flush.console()
    }

  }

  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")

  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}


fit.grem.Moment <- function(nu.init, S, ns,
                            max.ite = 1000, eps = 1e-3,
                            verbose = FALSE) {
  p <- nrow(S[[1]])
  interval <- c(p - 1 + 1e-5, 1e6)

  nu.old <- nu.init
  Psi.old <- Momentstep(nu = nu.old, S = S, ns = ns)

  for (i in seq_len(max.ite)) {
    nu.new  <- get.nu(Psi.old, nu.old, S, ns, interval)
    fac     <- (nu.new - p - 1)/(nu.old - p - 1)
    Psi.new <- Psi.old * fac

    if (verbose) {
      cat("ite =", i, ":", "nu.new - nu.old =", nu.new - nu.old,
          "fac =", fac, "nu =", nu.new, "\n");
      flush.console()
    }

    if (abs(nu.new - nu.old) < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old  <- nu.new
    }

  }
  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")
  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}



#
# WDA
#

# dgrem3 <- function(x, mu, Psi, nu, log = FALSE) {
#   stopifnot(length(x) == length(mu))
#   p <- nrow(Psi)
#
#   n1 <- nu/2*logdet(Psi)
#   n2 <- lgammap((nu + 1)/2, p = p)
#   d1 <- -p/2*log(pi)
#   d2 <- (nu + 1)/2 * logdet(Psi + tcrossprod(x - mu))
#   d3 <- lgammap(nu/2, p = p)
#
#   ans <- n1 + n2 - d1 - d2 - d3
#
#   if (!log) {
#     ans <- exp(ans)
#   }
#   attributes(ans) <- NULL
#   return(ans)
# }
#
#
# dgrem2 <- function(x, mu, Psi, nu, log = FALSE) {
#   stopifnot(length(x) == length(mu))
#   p <- nrow(Psi)
#   t1 <- lgammap((nu + 1)/2, p = p)
#   t2 <- -lgammap(nu/2, p = p)
#   t3 <- p/2*log(pi)
#   t4 <- -1/2*logdet(Psi)
#   t5 <- -(nu + 1)/2 * log(1 + sum(solve(Psi) * tcrossprod(x - mu)))
#
#   ans <- t1 + t2 + t3 + t4 + t5
#
#   if (!log) {
#     ans <- exp(ans)
#   }
#   attributes(ans) <- NULL
#   return(ans)
# }

dgrem <- function(x, mu, Psi, nu, log = FALSE) {
  p <- nrow(Psi)
  if (is.null(dim(x))) {
    dim(x) <- c(1, p)
  }
  Q <- function(x, A) {
    rowSums(tcrossprod(x, A) * x)
  }
  stopifnot(ncol(x) == length(mu))
  t1 <- lgammap((nu + 1)/2, p = p)
  t2 <- -lgammap(nu/2, p = p)
  t3 <- p/2*log(pi)
  t4 <- -1/2*logdet(Psi)

  x.center <- t(t(x) - mu)
  t5 <- -(nu + 1)/2 * log(1 + Q(x.center, solve(Psi)))

  ans <- t1 + t2 + t3 + t4 + t5

  if (!log) {
    ans <- exp(ans)
  }
  attributes(ans) <- NULL
  return(ans)
}

# drem(x = 1:4, mu = 1:4, Psi, nu)
# drem2(x = 1:4 + 0.1, mu = 1:4, Psi, nu)
# drem3(x = 1:4 + 0.1, mu = 1:4, Psi, nu)
# drem3(x = rbind(1:4, 1:4, 1:4 + 0.1), mu = 1:4, Psi, nu)

#
# TESTING
#

library(MCMCpack)
library(GMCM)

# Simulate data
test.grem <- function(k = 4, n = 10, ns = rep(n, k), p = 10,
                      nu = 15, Psi) {
  stopifnot(nu > p + 1)
  rwishart <- function(n, Sigma) {
    X <- rmvnormal(n = n, mu = rep(0, nrow(Sigma)), sigma = Sigma)
    return(crossprod(X))
  }

  if (missing(Psi)) {  # Compound symmetry matrix
    rho <- 0.5
    std <- 1
    Psi <- matrix(rho*std^2, p, p) + diag(rep((1 - rho)*std^2, p))
  }

  # Simulate data
  Sigmas <- replicate(length(ns), riwish(nu, Psi), simplify = FALSE)
  S <- lapply(seq_along(ns), function(i) rwishart(ns[i], Sigmas[[i]]))
  SS <- lapply(seq_along(S), function(i) S[[i]]/ns[i])

  # Fit model
  t1 <- system.time({
    res1 <- fit.grem.EM(Psi.init = diag(p), nu.init = p + 1, S, ns, eps = 1e-2)
  })
  t2 <- system.time({
    res2 <- fit.grem.MLE(nu.init = p + 1, S = S, ns = ns, eps = 1e-2)
  })
  t3 <- system.time({
    res3 <- fit.grem.Moment(nu.init = p + 1e7, S = S, ns = ns, eps = 1e-2)
  })

  expected.covariance <- Psi/(nu - p - 1)
  mean.covariance     <- Reduce("+", SS)/length(ns)
  grem.em.covariance  <- res1$Psi/(res1$nu - p - 1)
  grem.mle.covariance <- res2$Psi/(res2$nu - p - 1)
  grem.mom.covariance <- res3$Psi/(res3$nu - p - 1)

  return(list(expected = expected.covariance,
              mean = mean.covariance,
              grem.em  = grem.em.covariance,
              grem.mle = grem.mle.covariance,
              grem.mom = grem.mom.covariance,
              res.em = res1,
              res.mle = res2,
              res.mom = res3,
              t1 = t1, t2 = t2, t3 = t3,
              Sigmas = Sigmas, S = S, k = k,
              ns = ns, nu = nu, Psi = Psi))
}


SSEs <- function(x) {

  Svar <- function(Sigma, n) {
    n*(Sigma^2 + tcrossprod(diag(Sigma)))
  }

  get.SSE <- function(O, E, n) {
    denom <- Svar(E, n)
    diff <- (E - O)^2/denom
    return(sum(diag(diff)) + sum(get.lower.tri(diff)))
  }


  n <- x$ns[1]
  sse.grem.em  <- get.SSE(x$grem.em,  x$expected, n = n)
  sse.grem.mle <- get.SSE(x$grem.mle, x$expected, n = n)
  sse.grem.mom <- get.SSE(x$grem.mom, x$expected, n = n)
  sse.mean     <- get.SSE(x$mean, x$expected, n = n)

  stopifnot(all(x$ns == x$ns[1]))
  return(c(n  = x$ns[1],
           k  = x$k,
           nu = x$nu,
           p  = nrow(x$Psi),
           SSE.grem.em  = sse.grem.em,
           SSE.grem.mle = sse.grem.mle,
           SSE.grem.mom = sse.grem.mom,
           SSE.mean     = sse.mean))

}




if (FALSE) {


  #
  # LDA/QDA/NDA
  #

  nda.fit <- function(classes, vars) {
    p <- ncol(vars) # Dimensionality
    counts <- table(classes) # Count observations in each class
    split.data <- split(vars, f = classes) # Split vars object by classes
    means <- t(sapply(split.data, colMeans))  # Compute means in each class
    # Compute the scatter matrix in each dataset
    S <- lapply(split.data, function(x) cov(as.matrix(x), method = "ML")*nrow(x))

    # Find "common" covariance
    res <- fit.grem(Psi.init = diag(p), nu.init = p, S = S, counts, eps = 1e-2)
    sigma <- res$Psi/(res$nu - p - 1)

    return(list(counts = counts,
                means = means,
                sigma = sigma,
                Psi = res$Psi,
                nu = res$nu))
  }

  nda.predict <- function(nda.fit, newdata) {
    K <- length(nda.fit$counts)
    probs <- as.numeric(nda.fit$counts/sum(nda.fit$counts))


    f <- function(k) {
      probs[k] * dgrem(x = newdata, mu = nda.fit$means[k, ],
                       Psi = nda.fit$Psi, nu = nda.fit$nu, )
    }

    scaled_dens <- sapply(seq_len(K), f)

    post <- scaled_dens/rowSums(scaled_dens)
    colnames(post) <- names(nda.fit$counts)
    pred.class <- apply(post, 1, which.max)

    return(list(class = pred.class,
                prob = post))

  }

  rm(list = ls(pattern = "^conf.|.pred$|.train$"))

  library("Bmisc")
  library("mlbench")

  data(Satellite)
  data(BreastCancer)
  data(Glass)
  data(iris)

  #data <- as.data.frame(Satellite)
  #data <- as.data.frame(BreastCancer)
  #names(Glass) <- gsub("Type", "classes", names(Glass))
  #data <- Glass



  data <- iris[,c(1,2,5)]
  plot(data, col = iris[,5], cex = 1.5, pch = 16)
  names(data) <- gsub("Species", "classes", names(data))
  #data <- as.data.frame(mlbench.smiley(n = 100))
  #data <- as.data.frame(mlbench.hypercube(n = 10000, d = 3, sd = 0.5))
  #data <- as.data.frame(mlbench.twonorm(n = 20000, d = 2))

  n <- 200
  acc <- nmis <-
    structure(rep(NA, 3*n), dim = c(n, 3),
              dimnames = list(NULL, c("LDA", "QDA", "NDA")))


  for (i in seq_len(n)) {

    get <- sort(sample.int(nrow(data), size = nrow(data)/2))

    train <- data[get, ]
    valid <- data[-get, ]

    K <- nlevels(data$classes)

    lda.train <- lda(classes ~ ., train, prior = rep(1/K, K))
    qda.train <- qda(classes ~ ., train, prior = rep(1/K, K))
    nda.train <- nda.fit(train$classes, train[, !grepl("classes", names(train))])

    lda.pred <- predict(lda.train, newdata = valid)
    qda.pred <- predict(qda.train, newdata = valid)
    nda.pred <- nda.predict(nda.train,
                            newdata = valid[, !grepl("classes", names(train))])

    conf.lda <- table(True = valid$classes, Pred = lda.pred$class)
    conf.qda <- table(True = valid$classes, Pred = qda.pred$class)
    conf.nda <- table(True = valid$classes, Pred = nda.pred$class)

    nMisclassified <- function(x) {
      sum(get.upper.tri(x)) + sum(get.lower.tri(x))
    }
    accuracy <- function(x) {
      sum(diag(x))/sum(x)
    }

    nmis[i, ] <- c(nMisclassified(conf.lda),
                   nMisclassified(conf.qda),
                   nMisclassified(conf.nda))
    acc[i, ] <- c(accuracy(conf.lda), accuracy(conf.qda), accuracy(conf.nda))

    cat(i, "\n"); flush.console()
  }

  x <- seq(4.4,8.4, by = 0.01)
  y <- seq(2,4.5, by = 0.01)
  newdata <- expand.grid(x,y)
  colnames(newdata) <- c("Sepal.Length", "Sepal.Width")

  lda.reg <- predict(lda.train, newdata = newdata)
  qda.reg <- predict(qda.train, newdata = newdata)
  nda.reg <- nda.predict(nda.train, newdata = newdata)

}
