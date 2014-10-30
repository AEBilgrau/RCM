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

# Multicore support
library("foreach")
library("doMC") #library("doParallel") # Use this package on windows
registerDoMC(detectCores())

## ---- auxiliary_functions ----

#
# HDA
#

# Density of the GREM model
dgrem <- correlateR::dgrem

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
  S  <- lapply(seq_along(ns), function(i) rwishart(ns[i], Sigmas[[i]]))
  SS <- lapply(seq_along(S), function(i) S[[i]]/ns[i])

  # Fit model
  t1 <- system.time({
    res1 <- fit.grem(Psi.init = diag(p), nu.init = p + 1, S, ns, eps = 1e-2)
  })
  t2 <- system.time({
    res2 <- fit.grem.MLE(nu.init = p + 1, S = S, ns = ns, eps = 1e-2)
  })
  t3 <- system.time({
    res3 <- fit.grem.moment(nu.init = p + 1e7, S = S, ns = ns, eps = 1e-2)
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
  sse.mean     <- get.SSE(x$mean,     x$expected, n = n)

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
  Psi.init <- Reduce("+", S)/(nrow(vars) - length(S))

  # Find "common" covariance
  res <- fit.grem(Psi.init = Psi.init, nu.init = p, S = S, ns = counts, ...)
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
  post <- scaled_dens/rowSums(scaled_dens)
  colnames(post) <- names(hda.fit$counts)
  pred.class <- apply(scaled_dens, 1, which.max)
  return(list(class = pred.class))#, prob = post))
}



## ---- numerical_experiment ----
par.ne <- list(k = 3,
               nu = 15,
               p = 10,
               n.sims = 1000,
               n.obs = seq(4, 10, by = 1))
#rm(res)
if (!exists("res") | recompute) {
  set.seed(64403101)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.ne$n.obs)) {
    tmp <- foreach(j = seq_len(par.ne$n.sims)) %dopar% {
      test.grem(k = par.ne$k, n = par.ne$n.obs[i], p = par.ne$p, nu = par.ne$nu)
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.ne$n.obs), "done after",
        (proc.time()-st)[3] %/% 60, "mins.\n")
  }
  rm(tmp)
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

# par.ne$nu
# str(res[[1]])
# res[[1]]$res.mle$nu
# x <- res[[1]]
#
# get.nu <- function(x) {
#   c(n  = x$ns[1], k  = x$k, nu = x$nu, p  = nrow(x$Psi),
#     em = x$res.em$nu, mle = x$res.mle$nu, mom = x$res.mom$nu)
# }
# get.psi.diag <- function(x) {
#   cbind(n  = x$ns[1], k  = x$k, nu = x$nu, p  = nrow(x$Psi),
#         em = mean(diag(x$res.em$Psi))
#         mle = mean(diag(x$res.mle$nu$Psi)),
#         mom = diag(x$res.mom$nu)
# }
#

# nu.res <- as.data.frame(t(sapply(res, get.nu)))
#
# boxplot(em ~ n, data = nu.res, col = "grey", outline = FALSE)
# points(jitter(nu.res$n) - nu.res$k, nu.res$em, cex = 0.2, pch = 16)
# abline(h = par.ne$nu)
#
#
#
#
#
# abline(v = par.ne$nu)
# (sapply(res, function(x) x[["res.mle"]][["nu"]]))
# res$

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

par.xda <- list(K = 3,
                n.obs = 40,
                n.obs.valid = 100,
                n.runs = 100, #500,
                p.dims = c(5, 10, 20, 35))

#set.seed(15)
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
        evec <- eigen(theta$sigma[[1]])$vectors
        eval <- eigen(theta$sigma[[1]])$values
        w2 <- ( seq(sqrt(0.5), sqrt(4), length.out = p)^2) # reverse t
        w3 <- (-seq(sqrt(0.5), sqrt(3), length.out = p)^2) # reverse
        mu1 <- rep(0, p)
        mu2 <- c(evec %*% (sqrt(eval)*w2))
        mu3 <- c(evec %*% (sqrt(eval)*w3))
        theta$mu <- list(mu1, mu2, mu3)
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
        hda.train <- hda.fit(train$classes, select(train, -classes), eps = 1e-1)

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
misclass.df <- as.data.frame(misclass.df)


misclass.wide <-
  reshape(misclass.df,
          idvar = c("method", "spherical", "equal"),
          drop = c("mean", "sd", "is.min", "is.max"),
          timevar = "p", v.names = "risk", direction = "wide")
names(misclass.wide) <-
  gsub("risk\\.([0-9]+)", "$p = \\1$", names(misclass.wide))

misclass.wide.is.min <-
  reshape(as.data.frame(misclass.df), idvar = c("method", "spherical", "equal"),
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
l <- function(psi, k = 1, nu = 1, ni = 1, xi = 1) {
  k*nu/2*log(psi) - (nu + ni)/2*log(psi + xi^2)
}
dl <- function(psi, k = 1, nu = 1, ni = 1, xi = 1)  {
  k*nu/2*1/psi - (nu + ni)/2 * 1/(psi + xi^2)
}
par(mfrow = 1:2, mar = c(2, 2, 2, 0) + 0.2)
psi <- seq(0.5, 10, by = 0.01)
plot(psi, l(psi), type = "l", col = "red", lwd = 2, main = "loglik")
abline(v = 1, col = "grey", lty = 2, lwd = 2)
plot(psi, dl(psi), type = "l", col = "red", lwd = 2, main = "dloglik")
abline(h = 0, v = 1, col = "grey", lty = 2, lwd = 2)
## ---- end ----


## ---- log_gamma_ratio ----
logGammaRatio <- function(x, a) {
  lgamma(x + a) - lgamma(x) #= log(gamma(x + a)/gamma(x))
}
xs <- seq(0.01, 2, by = 0.01)
par(mfrow = c(1,2), mar = c(2, 2, 0, 0)+ 0.5)
plot(xs, logGammaRatio(xs, a = 2),    type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2)
plot(xs, logGammaRatio(xs, a = 1e-3), type = "l", xlab = "", ylab = "",
     col = "red", lwd = 2)
## ---- end ----

#
# DLBCL analysis
#

## ---- dlbcl ----
library("WGCNA")
library("correlateR")
library("Bmisc")
library("affy")
load("saved.RData")
load("studies.RData")
dlbcl.par <- list(top.n = 100)
if (!exists("dlbcl.grem") | !exists("dlbcl.pool") |
      !exists("var.pool") | recompute | TRUE) {
  load("gep.ensg.RData")
  vars      <- sapply(gep, function(x) rowSds(exprs(x))*(ncol(x) - 1))
  var.pool  <- rowSums(vars)/(sum(sapply(gep, ncol)) - length(gep))
  use.genes <- names(sort(var.pool, decreasing = TRUE)[seq_len(dlbcl.par$top.n)])
  gep <- lapply(gep, function(x) exprs(x)[use.genes, ])
  ns  <- sapply(gep, ncol)
  S   <- lapply(gep, function(x) scatter(t(x)))
  dlbcl.grem <- fit.grem(S = S, ns = ns, verbose = TRUE)
  dlbcl.pool <- Reduce("+", S)/sum(ns - 1)
  dimnames(dlbcl.grem$Psi) <- dimnames(dlbcl.pool)
  resave(dlbcl.grem, dlbcl.pool, use.genes, var.pool, file = "saved.RData")
}

dlbcl.exp <- with(dlbcl.grem, Psi2Sigma(Psi, nu))

## ---- dlbcl_plot_1 ----
# hist((get.lower.tri(dlbcl.adjMat)), main = "cor", breaks = 50,
#      col = "grey", prob = TRUE)
dlbcl.adjMat <- abs(cov2cor(dlbcl.exp))
hist(-log(get.lower.tri(dlbcl.adjMat)), main = "-log(cor)", breaks = 50,
     col = "grey", prob = TRUE)


## ---- dlbcl_plot_2 ----
dlbcl.tom <- TOMdist(adjMat = dlbcl.adjMat)
dimnames(dlbcl.tom) <- dimnames(dlbcl.adjMat) # Keep names
dlbcl.clust  <- flashClust(as.dist(dlbcl.tom), method = "average")
dlbcl.modules <- labels2colors(cutreeDynamicTree(dlbcl.clust, minModuleSize=9))
names(dlbcl.modules) <- dlbcl.clust$labels


o <- dlbcl.clust$order
layout(rbind(c(0,0,5,6), c(0,0,2,6), c(4,1,3,6)),
       widths = c(4, 1, 15, 15), heights = c(4, 1, 15))
TOMplot(dissim = dlbcl.adjMat, dendro = dlbcl.clust,
        Colors = dlbcl.modules, setLayout = FALSE)

layout.custom <- function(graph,...) {
  l <- layout.circle(graph)
  layout.fruchterman.reingold(graph, niter = 10000,
                              area = vcount(graph)/2,
                              maxdelta = 10*vcount(graph),
                              repulserad = vcount(graph),
                              weights = 10*E(graph)$weight,
                              start = l,
                              ...)
}
get.size <- function(x) {
  s <- rowSums(x)
  return((s - min(s))/max(s))
}
tmp <- soft(dlbcl.adjMat[o,o], .175)
gr <- plotModuleGraph(tmp,
                      labels = "",
                      diff.exprs = get.size(tmp)*3 + 3,
                      layout = layout.custom,
                      mark.shape = .5,
                      ecol = "black",
                      vcol = dlbcl.modules[o])
scale.layout <- function(x) {
  xx <- apply(x, 2, function(x) (x - min(x))/max(x - min(x)))
  return(2*xx - 1)
}
points(scale.layout(gr$layout), pch = 16, cex = 0.7)

# dlbcl.g <- graph.adjacency(tmp, mode = "undirected",
#                            weighted = TRUE, diag = FALSE)
# library("ape")
# plotHierarchicalEdgeBundles(as.phylo(dlbcl.clust),
#                             dlbcl.g,
#                             beta = 0.4,
#                             use.mrca = FALSE,
#                             simplify = FALSE,
#                             type = "fan",
#                             debug = FALSE,
#                             args.lines = list(col = "#FF000021", lwd = 2))
# aa <- plot.phylo(as.phylo(dlbcl.clust), type = "fan")

# ## ---- dlbcl_go_analysis ----
# library("biomaRt")
# if (!exists("gene.info") | recompute) {
#   mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#   attributes <- c("hgnc_symbol", "chromosome_name", "start_position",
#                   "end_position", "strand", "band", "ensembl_gene_id", "go_id",
#                   "name_1006", "definition_1006", "go_linkage_type",
#                   "namespace_1003")
#   gene.info <-
#     getBM(attributes = attributes,
#           filters = "arrayexpress", #"ens_hs_gene",
#           values = gsub("_at$", "", names(sort(var.pool, decreasing = TRUE))),
#           mart = mart)
#
#   resave(gene.info, file = "saved.RData")
# }
#
# # Mappings between ENSG and GO id
# gid2go <- split(gene.info$go_id, gene.info$ensembl_gene_id)
#
#
# dlbcl.modules
#
#
#
# 1-phyper(q  = 69, m = 1870, n = 10000 - 1870, k = 290)

