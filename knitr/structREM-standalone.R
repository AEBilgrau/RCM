setwd("~/Documents/PhD/paper-structREM/knitr/")

## ---- initialize ----
rm(list = ls())
recompute <- FALSE
library("Hmisc")
library("Bmisc")
library("correlateR")
library("GMCM")
library("MASS")
library("dplyr")
library("MCMCpack")
load("saved.RData")


## ---- auxiliary_functions ----

# Compute the log determinant easily
logdet <- function(x, ...) {
  z <- determinant(x, logarithm = TRUE, ...)
  stopifnot(z$sign == 1)
  return(z$modulus)
}

# Log-likelihood of GREM model
loglik <- function(Psi, nu, S, ns) {
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

# Derivative of the log-likelihood of GREM model
dloglik <- function(Psi, nu, S, ns) {  #
  k <- length(S)
  p <- nrow(S[[1]])
  t1 <-  k/2 * logdet(Psi)
  t2 <-  1/2 * sum(sapply(ns, function(ni) digammap((nu + ni)/2, p = p)))
  t3 <- -1/2 * sum(sapply(S, function(s) logdet(Psi + s)))
  t4 <- -k/2 * digammap(nu/2, p = p)
  return(t1 + t2 + t3 + t4)
}

# Optimize loglik wrt nu
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

# Compute new Psi from Psi, nu, S, ns using the EM step
EMstep <- function(Psi, nu, S, ns) {
  k <- length(S)
  p <- nrow(S[[1]])
  co <- 1/(k*nu)
  t <- lapply(seq_along(S), function(i) (co*(ns[i] + nu))*solve(Psi + S[[i]]))
  Psi_new <- solve(Reduce("+", t))  # Sum the matrices in t
  return(Psi_new)
}

# Compute new Psi from nu, S, ns using moment estimate
Momentstep <- function(nu, S, ns) {
  k <- length(ns)
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) S[[i]]/ns[i]))/k
  p <- nrow(Psi)
  fac <- nu - p - 1
  return(fac*Psi)
}

# Compute new Psi from nu, S, ns using approximate MLE
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


# MLE alg
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

# Estimation using moment
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
# HDA
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

# Denisty of the GREM model
dgrem <- function(x, mu, Psi, nu, logarithm = FALSE) {
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

  if (!logarithm) {
    ans <- exp(ans)
  }
  attributes(ans) <- NULL
  return(ans)
}

# drem(x = 1:4, mu = 1:4, Psi, nu)
# drem2(x = 1:4 + 0.1, mu = 1:4, Psi, nu)
# drem3(x = 1:4 + 0.1, mu = 1:4, Psi, nu)
# drem3(x = rbind(1:4, 1:4, 1:4 + 0.1), mu = 1:4, Psi, nu)


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

to.df <- function(sim) {
  return(data.frame(sim$z, classes = factor(sim$K)))
}

misclassificationRisk <- function(x) {
  s <- sum(x)
  return((s - sum(diag(x)))/s)
}

accuracy <- function(x) {
  return(sum(diag(x))/sum(x))
}

#
# HDA
#

hda.fit <- function(classes, vars, ...) {
  p <- ncol(vars) # Dimensionality
  counts <- table(classes) # Count observations in each class
  split.data <- split(vars, f = classes) # Split vars object by classes
  means <- t(sapply(split.data, colMeans))  # Compute means in each class
  # Compute the scatter matrix in each dataset
  S <- lapply(split.data, function(x)
    cov(as.matrix(x, nrow = nrow(x)), method = "ML")*nrow(x))

  #Psi.init <- diag(p)
  Psi.init <- Reduce("+", S)/(nrow(vars) - length(S))

  # Find "common" covariance
  res <- fit.grem.EM(Psi.init = Psi.init, nu.init = p, S = S, counts,
                     eps = 1e-2, ...)
  sigma <- res$Psi/(res$nu - p - 1)

  return(list(counts = counts, means = means, sigma = sigma,
              Psi = res$Psi, nu = res$nu))
}

hda.predict <- function(hda.fit, newdata) {
  K <- length(hda.fit$counts)
  probs <- as.numeric(hda.fit$counts/sum(hda.fit$counts))

  f <- function(k) {
    probs[k] * dgrem(x = newdata, mu = hda.fit$means[k, ],
                     Psi = hda.fit$Psi, nu = hda.fit$nu)
  }
  scaled_dens <- sapply(seq_len(K), f)
  #post <- scaled_dens/rowSums(scaled_dens)
  #colnames(post) <- names(hda.fit$counts)
  #pred.class <- apply(post, 1, which.max)
  pred.class <- apply(scaled_dens, 1, which.max)
  return(list(class = pred.class))#, prob = post))
}
#
# #
# # QDA
# #
#
# # kappa2 <- function(A) {
# #   ev <- eigen(A)$val
# #   abs(max(ev)/min(ev))
# # }
# #
# # correctForStableInversion2 <- function(A, verbose = TRUE) {
# #   p <- ncol(A)
# #   while (rcond(A) < .Machine$double.eps) {
# #     if (verbose) {cat("+"); flush.console()}
# #     A <- A + 10*diag(.Machine$double.eps, p)
# #   }
# #   return(A)
# # }
#
# correctForStableInversion <- function(A) {
#   if (rcond(A) < .Machine$double.eps) {
#     eps <- .Machine$double.eps
#     spec <- eigen(A)
#     lambda.min <- min(spec$values)
#     lambda.max <- max(spec$values)
#     a <- (1e2*eps*lambda.max - lambda.min)/(1e2*eps + 1)
#     spec$values <- spec$values + a
#     A <- spec$vectors %*% (t(spec$vectors) * spec$values)
#   }
#   return(A)
# }
#
#
# qda.fit <- function(classes, vars) {
#   counts <- table(classes)
#   split.vars <- split(vars, f = classes)
#   means <- t(sapply(split.vars, colMeans))
#   sigmas <- list()
#   for (i in seq_along(counts)) {
#     tmp <- as.matrix(split.vars[[i]], nrow = nrow(split.vars[[i]]))
#     sigmas[[i]] <- stats::cov(tmp)
#   }
#   return(list(counts = counts, means = means, sigmas = sigmas))
# }
#
#
# qda.predict <- function(myqda, newdata) {
#   theta <- myqda
#   theta$means <- lapply(seq_along(theta$sigmas), function(i) theta$means[i,])
#   theta$counts <- as.vector(theta$counts)/sum(theta$counts)
#   theta <- c(list(m = length(theta$counts), d = nrow(theta$sigmas[[1]])), theta)
#   vars <- as.matrix(newdata)
#   names(theta) <- c("m", "d", "pie", "mu", "sigma")
#   theta$sigma <- lapply(theta$sigma, correctForStableInversion)
#   kap <- GMCM:::EStep(x = as.matrix(vars), theta = theta)
#   return(list(class = apply(kap, 1, which.max), posterior = kap))
# }
#
# #
# # LDA
# #
#
# # lda.fit2 <- function(classes, vars) {
# #   p <- ncol(vars)
# #   counts <- table(classes)
# #   K <- length(counts)
# #   split.vars <- split(vars, f = classes)
# #   means <- t(sapply(split.vars, colMeans))
# #   sigma <- matrix(0, p, p)
# #   for (i in seq_len(K)) {
# #     tmp <- as.matrix(split.vars[[i]], nrow = nrow(split.vars[[i]]))
# #     sigma <- sigma + (nrow(tmp) - 1)*stats::cov(tmp)
# #   }
# #   sigma <- sigma/(nrow(vars) - K)
# #   sigmas <- replicate(K, sigma, simplify = FALSE)
# #   return(list(counts = counts, means = means, sigmas = sigmas))
# # }
# qda2lda <- function(myqda) {
#   K <- length(myqda$counts)
#   p <- nrow(myqda$sigmas[[1]])
#   sigma <- matrix(0, p, p)
#   for (i in seq_len(K)) {
#     sigma <- sigma + (myqda$counts[i] - 1)*myqda$sigmas[[i]]
#   }
#   sigma <- sigma/(sum(myqda$counts) - K)
#   myqda$sigmas <- replicate(K, sigma, simplify = FALSE)
#   return(myqda)
# }
#
# lda.fit <- function(classes, vars) {
#   qda <- qda.fit(classes, vars)
#   return(qda2lda(qda))
# }
#
# lda.predict <- function(mylda, newdata) {
#   theta <- mylda
#   K <- length(theta$counts)
#   theta$means <- lapply(seq_len(K), function(i) theta$means[i,])
#   theta$counts <- as.vector(theta$counts)/sum(theta$counts)
#   theta <- c(list(m = length(theta$counts), d = nrow(theta$sigmas[[1]])), theta)
#   vars <- as.matrix(newdata)
#   names(theta) <- c("m", "d", "pie", "mu", "sigma")
#   theta$sigma <- lapply(theta$sigma, correctForStableInversion)
#   stopifnot(is.theta(theta))
#   kap <- GMCM:::EStep(x = as.matrix(vars), theta = theta)
#   return(list(class = apply(kap, 1, which.max), posterior = kap))
# }




## ---- numerical_experiment ----
par.ne <- list(k = 3,
               nu = 15,
               p = 10,
               n.sims = 1000,
               n.obs = 4 + seq_len(10))
#rm(res)
if (!exists("res") | recompute) {
  set.seed(644031)
  st <- proc.time()
  res <- list()
  it <- 1
  for (i in seq_along(par.ne$n.obs)) {
    for (j in seq_len(par.ne$n.sims)) {
      res[[it]] <- test.grem(k = par.ne$k, n = par.ne$n.obs[i],
                             p = par.ne$p, nu = par.ne$nu)
      it <- it  + 1
      cat("it =", it, "of", length(par.ne$n.obs)*par.ne$n.sims, "done after",
          (proc.time()-st)[3] %/% 60, "mins.\n")
      flush.console()
    }
  }
  resave(res, file = "saved.RData")
}
## ---- end ----

## ---- numerical_experiment_plot ----
df <- as.data.frame(t(sapply(res, SSEs)))
df <- aggregate(cbind(SSE.grem.em, SSE.grem.mle, SSE.grem.mom, SSE.mean) ~
                  n + nu + k + p, mean, data = df)
plot(df$n, df$SSE.grem.em, col = "blue", type = "b", axes = FALSE,
     xlab = expression(n = n[i]),
     ylim = range(df[,-(1:4)]),
     ylab = "average SSE", pch = 15, lty = 1,
     main = paste0("k = ", df$k, ", nu = ", df$nu, ", p = ", df$p)[1])
lines(df$n, df$SSE.mean, col = "red", type = "b", pch = 16, lty = 2)
lines(df$n, df$SSE.grem.mle, col = "orange", type = "b", pch = 17, lty = 3)
lines(df$n, df$SSE.grem.mom, col = "green", type = "b", pch = 18, lty = 4)
legend("topright", legend = c("GREM (EM)", "pool",
                              "GREM (MLE)", "GREM (Moment)"),
       lty = 1:4, pch = c(15, 16, 17, 18), lwd = 2,
       col = c("blue", "red", "orange", "green"))
axis(1)
axis(2)
grid()
## ---- end ----

#
# HDA vs LDA vs QDA
#

# Testing
# theta <- rtheta()
# theta$pie <- c(1,1,1)/3
#
# train <- to.df(SimulateGMMData(n = 100, theta = theta))
# valid <- to.df(SimulateGMMData(n = 1000, theta = theta))
#
# myxda.train <- lda.fit(train$classes, vars = select(train, -classes))
# myxda.pred  <- lda.predict(myxda.train, select(valid, -classes))
#
# xda.train <- lda(classes ~ ., data = train)
# xda.pred  <- predict(xda.train, newdata = valid)
#
# table(xda.pred$class, myxda.pred$class)
#
# plot(xda.pred$posterior[,1], myxda.pred$posterior[,1], cex = .2)
# for (i in 2:3)
#   points(xda.pred$posterior[,i], myxda.pred$posterior[,i], cex = .2, col = i)
# abline(0, 1)


## ---- discriminant_analysis ----
e <- function(i, p) { # ith standard basis vector of length p
  vec <- rep(0, p)
  vec[i] <- 1
  return(vec)
}
set.seed(10)
par.xda <- list(K = 3,
                n.obs = 40,
                n.obs.valid = 100,
                n.runs = 500,
                p.dims = c(5, 10, 20, 35))
#rm(misclassification.risks)
if (!exists("misclassification.risks") | recompute) {

  inner <- structure(vector("list", 4),
                     names = c("eq.sph", "neq.sph",
                               "eq.ell", "neq.ell"))
  misclassification.risks <-
    replicate(length(par.xda$p.dims), inner, simplify = FALSE)
  names(misclassification.risks) <- paste0("p", par.xda$p.dims)

  st <- proc.time()
  for (p.index in seq_along(par.xda$p.dims)) {

    p <- par.xda$p.dims[p.index]

    for (s in seq_len(4)) {
      cat("p =", p, "and s =", s, "started at",
          (proc.time()-st)[3] %/% 60, "mins.\n")
      flush.console()

      method <- switch(s,
                       "EqualSpherical",
                       "UnequalSpherical",
                       "EqualEllipsoidal",
                       "UnequalEllipsoidal")
      theta <- rtheta(m = par.xda$K, d = p, method = method)
      theta$pie <- c(1, 1, 1)/3
      theta$mu <- list(rep(0, p), 3*e(1, p), 4*e(1, p))

      if (method == "EqualEllipsoidal") { # Low-variance subspace
        e <- eigen(theta$sigma[[1]])$vectors
        l <- eigen(theta$sigma[[1]])$values
        w2 <- ( seq(sqrt(0.5), sqrt(4), length.out = d)^2)
        w3 <- (-seq(sqrt(0.5), sqrt(3), length.out = d)^2)
        mu1 <- rep(0, d)
        mu2 <- c(e %*% (sqrt(l)*w2))
        mu3 <- c(e %*% (sqrt(l)*w3))
        theta$mu <- list(mu1, mu2, mu3)
#         theta$mu <-
#           list(rep(0, p),
#                tcrossprod(s, t(seq(2, 0, length.out = p)^4)),
#                tcrossprod(s, t(seq(0, 2, length.out = p)^4)))
      }

      mis.risk <-
        structure(rep(NA, 3*par.xda$n.runs), dim = c(par.xda$n.runs, 3),
                  dimnames = list(NULL, c("LDA", "QDA", "HDA")))

      for (i in seq_len(par.xda$n.runs)) {
        cat("i =", i, " "); flush.console()

        repeat {
          train <- to.df(SimulateGMMData(par.xda$n.obs,       theta = theta))
          valid <- to.df(SimulateGMMData(par.xda$n.obs.valid, theta = theta))
          # Make sure that we have a least two observations in each group
          if (all(table(train$classes) > 2) && all(table(valid$classes) > 2)) {
            break
          }
        }

        qda.train <- qda.fit(train$classes, select(train, -classes))
        lda.train <- qda2lda(qda.train)
        hda.train <- hda.fit(train$classes, select(train, -classes))

        lda.pred <- lda.predict(lda.train, newdata = select(valid, -classes))
        qda.pred <- qda.predict(qda.train, newdata = select(valid, -classes))
        hda.pred <- hda.predict(hda.train, newdata = select(valid, -classes))

        conf.lda <- table(True = valid$classes, Pred = lda.pred$class)
        conf.qda <- table(True = valid$classes, Pred = qda.pred$class)
        conf.hda <- table(True = valid$classes, Pred = hda.pred$class)

        mis.risk[i, ] <- c(misclassificationRisk(conf.lda),
                           misclassificationRisk(conf.qda),
                           misclassificationRisk(conf.hda))
        cat("\n")
      }
      misclassification.risks[[p.index]][[s]] <- mis.risk

    }
  }

  resave(misclassification.risks, file = "saved.RData")
}

## ---- hda_table_res ----
# Formatting results in misclassification.risks
cm  <- rapply(misclassification.risks, colMeans)
csd <- rapply(misclassification.risks, colSds)
tmp <- t(simplify2array(strsplit(names(cm), "\\.")))

# Format
misclass.df <-
  data.frame(t(simplify2array(strsplit(names(cm), "\\.")))) %>%
  dplyr::select(p = X1, equal = X2, spherical = X3, method = X4) %>%
  mutate(p = gsub("p", "", p)) %>%
  mutate(risk = sprintf("%.03f (%.3f)", cm, csd), mean = cm, sd = csd) %>%
  mutate(risk = gsub("0\\.", ".", risk)) %>%
  group_by(p, equal, spherical) %>%
  mutate(is.min = mean==min(mean), is.max = mean==max(mean))


misclass.wide <-
  reshape(misclass.df, idvar = c("method", "spherical", "equal"),
          drop = c("mean", "sd", "is.min", "is.max"),
          timevar = "p", v.names = "risk", direction = "wide")
names(misclass.wide) <-
  gsub("risk\\.([0-9]+)", "$p = \\1$", names(misclass.wide))

misclass.wide.is.min <-
  reshape(misclass.df, idvar = c("method", "spherical", "equal"),
          drop = c("mean", "sd", "risk", "is.max"),
          timevar = "p", v.names = "is.min", direction = "wide")

misclass.wide.is.max <-
  reshape(misclass.df, idvar = c("method", "spherical", "equal"),
          drop = c("mean", "sd", "risk", "is.min"),
          timevar = "p", v.names = "is.max", direction = "wide")

tmp <- misclass.wide[, -(1:3)]
tmp.min <- as.matrix(misclass.wide.is.min[, -(1:3)])
tmp.max <- as.matrix(misclass.wide.is.max[, -(1:3)])
cmd <- ifelse(tmp.min, "green", ifelse(tmp.max, "red", ""))

df.rownames <-
  sprintf("%s, %s",
          ifelse(misclass.wide$equal == "eq", "Equal", "Unequal"),
          ifelse(misclass.wide$spherical == "sph", "spherical", "ellipsoidal"))


# Create LaTeX table
latex(tmp, file = "",
      title = "$\\vSigma_1, ..., \\vSigma_3$",
      rowname = gsub("NDA", "HDA", misclass.wide$method),
      n.rgroup = table(df.rownames)[unique(df.rownames)],
      rgroup = unique(df.rownames),
      n.cgroup = ncol(misclass.wide) - 3,
      cgroup = "Mean misclassification risk (sd)",
      cellTexCmds = cmd,
      label = "HDA_tab",
      caption.loc = "bottom",
      caption = paste("The results for the different scenarios. The minimum",
                      "and maximum risk are highlighted in green and red,",
                      "respectively."))
rm(cm, csd, tmp, df.rownames)
## ---- end ----


## ---- one_dimensional_loglik ----
par(mar = c(2,2,0,0) + 0.2)
dl <- function(phi, k = 1, nu = 1, ni = 1, xi = 1)  {
  k*nu/2*1/phi - (nu + ni)/2 * 1/(phi + xi^2)
}
phi <- seq(0.5, 10, by = 0.01)
plot(phi, dl(phi), type = "l", col = "red", lwd = 2)
abline(h = 0, col = "grey", lty = 2, lwd = 2)
## ---- end ----
