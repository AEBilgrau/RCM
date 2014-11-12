setwd("~/Documents/PhD/paper-RCM/knitr/")

## ---- initialize_script ----
rm(list = ls())
recompute <- FALSE
library("Hmisc")
library("Bmisc")
library("correlateR")
library("GMCM")
library("MASS")
library("WGCNA")
library("affy")
library("biomaRt")
library("ape")
library("dplyr")
library("topGO")
library("igraph")
library("MCMCpack")
if (file.exists("saved.RData")) load("saved.RData")

# Multicore support
library("foreach")
library("doMC") #library("doParallel") # Use this package on windows
registerDoMC(detectCores())

## ---- auxiliary_functions ----

#
# HDA
#

# Density of the RCM model
drcm <- correlateR::drcm

# Simulate data
test.rcm <- function(k = 4, n = 10, ns = rep(n, k), p = 10,
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
    res1 <- fit.rcm(Psi.init = diag(p), nu.init = p + 1, S, ns, eps = 1e-2)
  })
  t2 <- system.time({
    res2 <- correlateR:::fit.rcm.MLE(nu.init = p + 1, S = S,
                                     ns = ns, eps = 1e-2)
  })
  t3 <- system.time({
    res3 <- correlateR:::fit.rcm.moment(nu.init = p + 1e7, S = S,
                                        ns = ns, eps = 1e-2)
  })

  expected.covariance <- Psi/(nu - p - 1)
  mean.covariance     <- Reduce("+", SS)/length(ns)
  rcm.em.covariance  <- res1$Psi/(res1$nu - p - 1)
  rcm.mle.covariance <- res2$Psi/(res2$nu - p - 1)
  rcm.mom.covariance <- res3$Psi/(res3$nu - p - 1)

  return(list(expected = expected.covariance,
              mean = mean.covariance,
              rcm.em  = rcm.em.covariance,
              rcm.mle = rcm.mle.covariance,
              rcm.mom = rcm.mom.covariance,
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
  sse.rcm.em  <- get.SSE(x$rcm.em,  x$expected, n = n)
  sse.rcm.mle <- get.SSE(x$rcm.mle, x$expected, n = n)
  sse.rcm.mom <- get.SSE(x$rcm.mom, x$expected, n = n)
  sse.mean    <- get.SSE(x$mean,    x$expected, n = n)

  stopifnot(all(x$ns == x$ns[1]))
  return(c(n  = x$ns[1],
           k  = x$k,
           nu = x$nu,
           p  = nrow(x$Psi),
           SSE.rcm.em  = sse.rcm.em,
           SSE.rcm.mle = sse.rcm.mle,
           SSE.rcm.mom = sse.rcm.mom,
           SSE.mean    = sse.mean))
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
  res <- fit.rcm(Psi.init = Psi.init, nu.init = p, S = S, ns = counts, ...)
  sigma <- res$Psi/(res$nu - p - 1)

  return(list(counts = counts, means = means, sigma = sigma,
              Psi = res$Psi, nu = res$nu))
}

hda.predict <- function(hda.fit, newdata) {
  K <- length(hda.fit$counts)
  probs <- as.numeric(hda.fit$counts/sum(hda.fit$counts))

  f <- function(k) {
    probs[k] * drcm(x = newdata, mu = hda.fit$means[k, ],
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
if (!exists("res") | recompute) {
  set.seed(64403101)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.ne$n.obs)) {
    tmp <- foreach(j = seq_len(par.ne$n.sims)) %dopar% {
      test.rcm(k = par.ne$k, n = par.ne$n.obs[i], p = par.ne$p, nu = par.ne$nu)
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.ne$n.obs), "done after",
        (proc.time()-st)[3] %/% 60, "mins.\n")
  }
  rm(tmp)
  resave(res, file = "saved.RData")
}

## ---- numerical_experiment_plot ----
df <- as.data.frame(t(sapply(res, SSEs)))
df <- aggregate(cbind(SSE.rcm.em, SSE.rcm.mle, SSE.rcm.mom, SSE.mean) ~
                  n + nu + k + p, mean, data = df)
plot(df$n, df$SSE.rcm.em, col = "blue", type = "b", axes = FALSE,
     xlab = expression(n = n[i]),
     ylim = range(df[,-(1:4)]),
     ylab = "average SSE", pch = 15, lty = 1,
     main = paste0("k = ", df$k, ", nu = ", df$nu, ", p = ", df$p)[1])
lines(df$n, df$SSE.mean, col = "red", type = "b", pch = 16, lty = 2)
lines(df$n, df$SSE.rcm.mle, col = "orange", type = "b", pch = 17, lty = 3)
lines(df$n, df$SSE.rcm.mom, col = "green", type = "b", pch = 18, lty = 4)
legend("topright", legend = c("RCM (EM)", "pool",
                              "RCM (MLE)", "RCM (Moment)"),
       lty = 1:4, pch = c(15, 16, 17, 18), lwd = 2,
       col = c("blue", "red", "orange", "green"))
axis(1)
axis(2)
grid()
## ---- end ----


#
# Discriminant analysis
#

## ---- discriminant_analysis ----
e <- function(i, p) { # ith standard basis vector of length p
  vec <- rep(0, p)
  vec[i] <- 1
  return(vec)
}

par.xda <- list(K = 3,
                n.obs = 40,
                n.obs.valid = 100,
                n.runs = 2500,
                p.dims = c(5, 10, 20, 35))
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
      misclassification.risks[[p.index]][[s]]  <-
        tmp <- foreach(i = seq_len(par.xda$n.runs), .combine = rbind,
                       .packages = c("dplyr", "correlateR")) %dopar% {
                         cat("i =", i, "\n"); flush.console()

                         repeat {
                           train <- to.df(SimulateGMMData(par.xda$n.obs,       theta = theta))
                           valid <- to.df(SimulateGMMData(par.xda$n.obs.valid, theta = theta))
                           # Make sure that we have a least two observations in each group
                           if (all(table(train$classes) > 2) && all(table(valid$classes) > 2)){
                             break
                           }
                         }

                         qda.train <- qda.fit(train$classes, dplyr::select(train, -classes))
                         lda.train <- qda2lda(qda.train)
                         hda.train <-
                           hda.fit(train$classes, dplyr::select(train, -classes), eps=1e-2)

                         lda.pred <-
                           lda.predict(lda.train, newdata = dplyr::select(valid, -classes))
                         qda.pred <-
                           qda.predict(qda.train, newdata = dplyr::select(valid, -classes))
                         hda.pred <-
                           hda.predict(hda.train, newdata = dplyr::select(valid, -classes))

                         conf.lda <- table(True = valid$classes, Pred = lda.pred$class)
                         conf.qda <- table(True = valid$classes, Pred = qda.pred$class)
                         conf.hda <- table(True = valid$classes, Pred = hda.pred$class)

                         return(c(misclassificationRisk(conf.lda),
                                  misclassificationRisk(conf.qda),
                                  misclassificationRisk(conf.hda)))
                       }
      colnames(misclassification.risks[[p.index]][[s]]) <- c("LDA","QDA","HDA")

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
dev.off()

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
dev.off()

#
# DLBCL analysis
#

## ---- dlbcl_analysis ----
load("studies.RData")
load("gep.ensg.RData")
studies <- studies[studies$Study != "Celllines", ]
dlbcl.dims <- sapply(rev(gep)[-1], dim)

dlbcl.par <- list(top.n = 300,
                  linkage = "ward",
                  go.alpha.level = 0.01,
                  ontology = "MF",
                  minModuleSize = 20,
                  threshold = 0.3)

vars      <- sapply(gep, function(x) rowSds(exprs(x))*(ncol(x) - 1))
var.pool  <- rowSums(vars)/(sum(sapply(gep, ncol)) - length(gep))
use.genes <- names(sort(var.pool, decreasing = TRUE)[seq_len(dlbcl.par$top.n)])
gep <- lapply(gep, function(x) exprs(x)[use.genes, ])

if (!exists("dlbcl.rcm") | !exists("var.pool") | recompute) {
  dlbcl.ns  <- sapply(gep, ncol)
  dlbcl.S   <- lapply(gep, function(x) correlateR::scatter(t(x)))
  dlbcl.rcm <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE)
  dimnames(dlbcl.rcm$Psi) <- dimnames(dlbcl.S[[1]])
  resave(dlbcl.rcm, file = "saved.RData")
}
# Expected covariance matrix
#s <- sample(1:nrow(dlbcl.rcm$Psi), 45)
dlbcl.exp <- with(dlbcl.rcm, Psi2Sigma(Psi, nu))#[s,s]

## ---- dlbcl_plot_1 ----
dlbcl.cor <- cov2cor(dlbcl.exp)
dlbcl.adjMat <- abs(dlbcl.cor)

par(mfrow = c(1,2), mar = c(4,4,0,0) + .2)
hist(get.lower.tri(dlbcl.cor), col = "grey", breaks = 50, main = "",
     xlab = "correlation",  prob = TRUE)
hist(-log(get.lower.tri(dlbcl.adjMat)), breaks = 50, col = "grey", main = "",
     xlab = "-log(abs(correlation))",  prob = TRUE)



## ---- dlbcl_mappings ----
all.genes <- gsub("_at$", "", names(sort(var.pool, decreasing = TRUE)))
if (!exists("gene.info") || recompute) {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("hgnc_symbol", "chromosome_name", "start_position",
                  "end_position", "strand", "band", "ensembl_gene_id", "go_id",
                  "name_1006", "definition_1006", "go_linkage_type",
                  "namespace_1003")
  gene.info <- getBM(attributes = attributes, filters = "arrayexpress",
                     values = all.genes, mart = mart)
  resave(gene.info, file = "saved.RData")
}

# Mappings between ENSG and GO id
gid2go <- split(gene.info$go_id, gene.info$ensembl_gene_id)
gid2hugo <- with(dplyr::select(gene.info, hgnc_symbol, ensembl_gene_id) %>%
                   distinct(hgnc_symbol, ensembl_gene_id),
                 structure(hgnc_symbol, names = ensembl_gene_id))

map2hugo <- function(ensg) {
  ensg <- gsub("_at$", "", ensg)
  ans <- gid2hugo[ensg]
  isna <- is.na(ans)
  ans[isna] <- ensg[isna]
  names(ans) <- ensg
  return(ans)
}


## ---- dlbcl_clustering ----
dlbcl.clust <- flashClust(as.dist(1 - dlbcl.adjMat), method = dlbcl.par$linkage)
# Cluster

dlbcl.modules <- labels2colors(cutree(dlbcl.clust, k = 5))
names(dlbcl.modules) <- dlbcl.clust$labels

dlbcl.g <- graph.adjacency(dlbcl.cor, mode = "undirected",
                           weighted = TRUE, diag = FALSE)
w <- E(dlbcl.g)$weight
E(dlbcl.g)$color <- alp(ifelse(w < 0, "steelblue","tomato"), abs(w))
V(dlbcl.g)$name <- map2hugo(V(dlbcl.g)$name)

# Phylo
phylo <- as.phylo(dlbcl.clust)
phylo$tip.label <- map2hugo(phylo$tip.label)

## ---- end ----


## ---- dlbcl_plot_2 ----
dlbcl_plot_2 <- "figure/dlbcl_plot_2-1.png"
png(dlbcl_plot_2, width = 1500, height = 2200, res = 150)
if (!file.exists("figure/dlbcl_plot_2-1.png") || recompute) {

  # LAYOUT
  layout(rbind(c(0,0,5,6), c(0,0,2,6), c(4,1,3,6),7),
         widths = c(4, 1, 15, 15), heights = c(4, 1, 10, 20))

  # PANEL A
  TOMplot(dissim = abs(dlbcl.cor),
          dendro = dlbcl.clust,
          Colors = dlbcl.modules,
          setLayout = FALSE)

  # PANEL B
  layout.custom <- function(graph,...) {
    l <- layout.circle(graph)
    layout.fruchterman.reingold(graph, niter = 10000,
                                area = vcount(graph)/2,
                                maxdelta = vcount(graph),
                                repulserad = vcount(graph),
                                weights = E(graph)$weight,
                                start = l,
                                ...)
  }
  get.size <- function(x) {
    s <- rowSums(x)
    return((s - min(s))/max(s))
  }

  #o <- dlbcl.clust$order
  topn <- function(x, n = 6) {
    nms <- names(x)
    top <- names(tail(sort(x), n = n))
    ifelse(nms %in% top, nms, "")
  }

  thresholded <- soft(dlbcl.cor, dlbcl.par$threshold)
  gr <- plotModuleGraph(abs(thresholded),
                        labels = "", #map2hugo(topn(rowSums(abs(dlbcl.cor)))),
                        diff.exprs = 3*get.size(abs(dlbcl.cor)) + 3,
                        layout = layout.custom,
                        mark.shape = .5,
                        ecol = "black",
                        vcol = dlbcl.modules)
  scaleToLayout <- function(x) {
    return(2*apply(x, 2, function(x) (x - min(x))/max(x - min(x))) - 1)
  }
  points(scaleToLayout(gr$layout), pch = 16, cex = 0.7)

  # PANEL C
  plotHierarchicalEdgeBundles(phylo, dlbcl.g, beta = 0.95,
                              cex = 0.7, type = "fan",
                              tip.color = dlbcl.modules,
                              e.cols = alp(E(dlbcl.g)$color, 0.6))
}
dev.off()
## ---- end ----





## ---- dlbcl_go_analysis ----
# GO analysis
if (!exists("dlbcl.module.genes") || !exists("dlbcl.go.analysis") ||
      recompute) {
  dlbcl.module.genes <- list()
  dlbcl.go.analysis <- list()

  for (col in unique(dlbcl.modules)) {
    dlbcl.module.genes[[col]] <- mod.genes <-
      gsub("_at$", "", names(dlbcl.modules[dlbcl.modules == col]))
    geneList <- factor(as.integer(all.genes %in% mod.genes))
    names(geneList) <- all.genes
    GOdata <- new("topGOdata", ontology = dlbcl.par$ontology,
                  allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = gid2go)
    result.fisher <-
      runTest(GOdata, algorithm = "classic", statistic = "fisher")
    mod.res <- GenTable(GOdata, classicFisher = result.fisher, topNodes = 100)
    dlbcl.go.analysis[[col]] <-
      mod.res %>% filter(classicFisher <= dlbcl.par$go.alpha.level)
  }
  dlbcl.module.genes <- lapply(dlbcl.module.genes, map2hugo)
  resave(dlbcl.module.genes, dlbcl.go.analysis, file = "saved.RData")
}





## ---- dlbcl_mod_tab ----
mod.genes <- lapply(dlbcl.module.genes, function(x) paste0(names(x), "_at"))
# Order by rowSums
dlbcl.exp.sub <- lapply(mod.genes, function(ensg) dlbcl.exp[ensg, ensg])
dlbcl.exp.sub <- lapply(dlbcl.exp.sub, function(x) {
  o <- order(rowSums(x), decreasing = TRUE)
  x[o,o]
})
tmp <- lapply(dlbcl.exp.sub, rownames)
m <- max(table(dlbcl.modules))  # Get largest module

# Construct table
dlbcl.mod.tab.genes <-
  sapply(tmp, function(x) unname(c(map2hugo(x), rep(NA, m - length(x)))))

# First letter capitalized
cgroup <- gsub("(^|[[:space:]])([[:alpha:]])",
               "\\1\\U\\2", names(dlbcl.module.genes), perl = TRUE)
colnames(dlbcl.mod.tab.genes) <- # Number of gens in each module
  paste0("n = ", colSums(!is.na(dlbcl.mod.tab.genes)))

is.ensg <- structure(grepl("^ENSG", dlbcl.mod.tab.genes),
                     dim = dim(dlbcl.mod.tab.genes))
nn <- 40
latex(dlbcl.mod.tab.genes[seq_len(nn), ],
      cgroup = cgroup,
      size = "tiny",
      caption = paste("The identified modules and their sizes and selected",
                      "member genes. For each module, the genes are sorted",
                      "by their intra-module connectivity from highest to",
                      "lowest. Only a maximum of", nn, "genes is shown."),
      cellTexCmds = ifelse(is.ensg, "tiny", "")[seq_len(nn), ],
      caption.loc = "bottom",
      label = "tab:dlbcl_mod_tab",
      landscape = TRUE,
      file = "")


## ---- GO_tabs ----
go.table <- do.call(rbind, dlbcl.go.analysis)
go.col <- gsub("(^|[[:space:]])([[:alpha:]])",
               "\\1\\U\\2", gsub("\\.[0-9]+$", "", rownames(go.table)),
               perl = TRUE)
colnames(go.table)[3:6] <- c("$N$", "$O$", "$E$", "$P$-val")
rownames(go.table) <- NULL
latex(go.table[, -c(1, 3)],
      rowname = go.table[, 1],
      rgroup = unique(go.col),
      n.rgroup = table(go.col)[unique(go.col)],
      title = "GO ID",
      size = "tiny",
      caption = paste0("The significant GO terms in the modules at ",
                       "$\\alpha$-level ", dlbcl.par$go.alpha.level, "."),
      #caption.loc = "bottom",
      label = paste0("tab:GO_tabs"),
      longtable = TRUE,
      lines.page = 80,
      file = "")
## ---- end ----

## ---- survival_analysis ----

plot.cis <- function(coxph,...) {
  x95 <- summary(cph.fit, conf.int = c(0.95))
  x99 <- summary(cph.fit, conf.int = c(0.99))
  ci95 <- x95$conf.int
  ci99 <- x99$conf.int
  rng <- range(ci99)

  plot(1, type = "n", xlim = rng, ylim = c(0, nrow(ci99)+1),
       xlab = "", ylab = "", axes = FALSE, log = "x", ...)
  h <- 0.2
  col <- gsub("^ME", "", rownames(ci95))
  for (ci in  list(ci95, ci99)) {
    rect(ci[,3], 1:nrow(ci) - h, ci[,4], 1:nrow(ci) + h,
         col = alp(col, 0.5), border = NA)

    segments(x0 = ci[,1], y0 = 1:nrow(ci) - 1.5*h, y1 = 1:nrow(ci) + 1.5*h,
             lwd = 2)
  }
  rect(ci99[,3], 1:nrow(ci99) - h, ci99[,4], 1:nrow(ci99) + h)
  rect(ci95[,3], 1:nrow(ci95) - h, ci95[,4], 1:nrow(ci95) + h)
  abline(v = 1, lty = 2, col = "darkgrey")
  axis(3, at = axTicks(3), label = formatC(axTicks(3)))
  axis(2, at = 1:nrow(ci), label = col,
       las = 2, tick = FALSE, pos = axTicks(3)[1])
}

load("metadata.RData")
library("WGCNA")
library("survival")
library("rms")

par(mar = c(.1,5,5,0.1), oma = c(2,0,0,0))
layout(cbind(1:3,4:6), heights = c(1,2,2))

for (j in 1:2) {
  meta <- switch(j, metadataLLMPPCHOP, metadataLLMPPRCHOP)
  rownames(meta) <- as.character(meta$GEO.ID)
  expr <- switch(j,
                 (gep$GEPLLMPPCHOP.ensg)[names(dlbcl.modules), ],
                 (gep$GEPLLMPPRCHOP.ensg)[names(dlbcl.modules), ])
  meta <- meta[colnames(expr), ] # Reorder

  # Check order
  stopifnot(rownames(meta) == colnames(expr))

  res <- moduleEigengenes(t(expr), dlbcl.modules)
  eg <- res$eigengenes
  col <- gsub("^ME", "", colnames(eg))

  cph.fit <- coxph(meta$OS ~ ., data = eg, x = TRUE, y = TRUE)

  plot.cis(cph.fit, main = switch(j, "LLMPP CHOP", "LLMPP R-CHOP"))

  eg.hl <- data.frame(eg >= apply(eg, 2, median))
  for (i in seq_len(ncol(eg.hl))) {
    if (col[i] %in% c("turquoise", "yellow"))
      plot(survfit(meta$OS ~ factor(eg.hl[,i])), conf.int = TRUE, main = col[i],
           col = c(col[i], "black"), lwd = 2)
  }
}
## ----


