# This R-script contains the code to produce simulated data and run analyses
# to get results for the Paper. Please note that running the simulations,
# especially estimating p-values, takes a long time. With parallel 
# processing in 3 threads on an i7-5600U laptop, estimated runtime is ~1 week
# Rasmus Brøndum, June 2017


## Load libraries
library("gplots")
library("ggplot2")
library("ggbiplot")
library("correlateR")
library("dendextend")
library("Matrix")
library("matrixStats")
library("flashClust")
library("fpc")
library("cluster")
library("MASS")
library("parallel")
library("lsr")
library("Hmisc")
library("Bmisc")
library("GMCM")
library("affy")
library("biomaRt")
library("ape")
library("topGO")
library("MCMCpack")
library("dplyr")
library("tawny")
library("foreach")
library("WGCNA")
library("igraph")
library("gProfileR")

#########################
### Multicore support ###
#########################
library("foreach")
if (Sys.info()[1] == "Windows") {
  library("doParallel") # Use this package on windows
  registerDoParallel(detectCores())
} else {
  library("doMC")
  registerDoMC(detectCores())
}

#########################
### Initialize script ###
#########################
rm(list = ls())
set.seed(1)
recompute <- FALSE

if (file.exists("saved.RData"))
  loaded <- load("saved.RData")


##############################
### User defined functions ###
##############################
get.ICC <- function(x) {
  with(x, ICC(nu, nrow(Psi)))
}

## Functions for heatmap
my.dist <- function(X){
  as.dist(1 - abs(X))
}
my.hclust <- function(my.dist){
  flashClust(my.dist, method= "ward")
}

### Test RCM model on simulated data
test.rcm <- function(k = 4, n = 10, ns = rep(n, k),
                     p = 10, nu = 15, Psi, e = 1e-2, ...) {
  stopifnot(nu > p - 1)
  
  if (missing(Psi)) {
    rho <- 0.5  # Compound symmetry matrix
    std <- 1
    Psi <- matrix(rho*std^2, p, p) + diag(rep((1 - rho)*std^2, p))
  }
  
  # Create data
  warningHandler <- function(warning) {
    if (any(grepl("singular Wishart distribution", warning))) {
      invokeRestart("muffleWarning")
    }
  }
  S <- withCallingHandlers(createRCMData(ns = ns, psi = Psi, nu = nu),
                           warning = warningHandler)
  
  # Check if initial values are given
  args <- list(...)
  if (is.null(args[['nu.init']])) nu.init <- sum(ns) + ncol(S[[1]]) + 1
  if (is.null(args[['Psi.init']])) Psi.init <- nu.init*correlateR:::pool(S, ns)
  
  t_em   <- system.time(res.em   <-
                          factory(fit.rcm)(S, ns,
                                           Psi.init = Psi.init,
                                           nu.init = nu.init,
                                           method = "EM",   eps=e,
                                           ...))
  t_pool <- system.time(res.pool <-
                          factory(fit.rcm)(S, ns,
                                           Psi.init = Psi.init,
                                           nu.init = nu.init,
                                           method = "pool", eps=e,
                                           ...))
  t_mle  <- system.time(res.mle  <-
                          factory(fit.rcm)(S, ns,
                                           Psi.init = Psi.init,
                                           nu.init = nu.init,
                                           method = "appr", eps=e,
                                           ...))
  time <- c(em = t_em[3], pool = t_pool[3], mle = t_mle[3])
  
  return(list(S = S, ns = ns, nu = nu, Psi = Psi,
              rcm.em = res.em, rcm.pool = res.pool, rcm.mle = res.mle,
              time = time))
}

SSEs <- function(x) {
  Svar <- function(Sigma, n) {
    n*(Sigma^2 + tcrossprod(diag(Sigma)))
  }
  getSSE <- function(O, E, n) {
    denom <- Svar(E, n)
    diff <- (E - O)^2/denom
    return(sum(diag(diff)) + sum(get.lower.tri(diff)))
  }
  getSigma <- function(x) {
    with(x, Psi2Sigma(Psi, nu))
  }
  n <- x$ns[1]
  expected  <- Psi2Sigma(Psi = x$Psi, nu = x$nu)
  sse.rcm.em   <- getSSE(getSigma(x$rcm.em[[1]]),   expected, n = n)
  sse.rcm.mle  <- getSSE(getSigma(x$rcm.mle[[1]]),  expected, n = n)
  sse.rcm.pool <- getSSE(getSigma(x$rcm.pool[[1]]), expected, n = n)
  
  get <- c("nu", "iterations", "loglik")
  stopifnot(all(x$ns == x$ns[1]))
  return(c(n  = x$ns[1],
           k  = length(x$ns),
           nu = x$nu,
           p  = nrow(x$Psi),
           SSE.rcm.em   = sse.rcm.em,
           SSE.rcm.mle  = sse.rcm.mle,
           SSE.rcm.pool = sse.rcm.pool,
           nu.rcm.em = x$rcm.em[[1]]$nu,
           nu.rcm.mle = x$rcm.mle[[1]]$nu,
           nu.rcm.pool = x$rcm.pool[[1]]$nu,
           em = unlist(x$rcm.em[get]),
           mle = unlist(x$rcm.mle[get]),
           pool = unlist(x$rcm.pool[get]),
           time = x$time))
}

### Simulate data from a true sigma, and fit different model nSims times
test.rcm.simData <- function(nSims, sigmaSim, nuSim, n, k, linkage.method="ward"){
  
  ## Number of studies and obs in each study
  nk = rep(n,k) 
  
  ## Initialize list for results
  results <- list()
  
  ## Clustering of "True" data
  hcluSim  <- flashClust(as.dist(1-abs(cov2cor(sigmaSim))), method=linkage.method)
  
  for(sim in 1:nSims){
    cat(date(), "- Iteration ", sim, "of ", nSims, "with n =", n, "and nu =", nuSim, "\n")
    ## Simulate data
    simData   <- createRCMData(ns = nk, sigma = sigmaSim, nu = nuSim)
    
    ## Fit RCM model
    rcm       <- fit.rcm(simData,nk, method="EM", max.ite = 2500)
    sigma.rcm <- Psi2Sigma(rcm$Psi, rcm$nu)
    row.names(sigma.rcm) = colnames(sigma.rcm) = row.names(sigmaSim)
    hclu.rcm  <- flashClust(my.dist(cov2cor(sigma.rcm)), method=linkage.method)
    
    ## Fit pool model
    pool       <- fit.rcm(simData,nk, method="pool", max.ite=2500)
    sigma.pool <- Psi2Sigma(pool$Psi, pool$nu)
    row.names(sigma.pool) = colnames(sigma.pool) = row.names(sigmaSim)
    hclu.pool  <- flashClust(my.dist(cov2cor(sigma.pool)), method=linkage.method)
    
    ## Fit apprMLE model
    mle       <- fit.rcm(simData,nk, method="approxMLE", max.ite=2500)
    sigma.mle <- Psi2Sigma(mle$Psi, mle$nu)
    row.names(sigma.mle) = colnames(sigma.mle) = row.names(sigmaSim)
    hclu.mle  <- flashClust(my.dist(cov2cor(sigma.mle)), method=linkage.method)
    
    results[[sim]] <- list("nu"=nuSim,
                           "n" =n,
                           "k" =k,
                           "rcm.fit"=rcm,
                           "pool.fit"=pool,
                           "mle.fit"=mle,
                           "sigma.rcm"=sigma.rcm,
                           "sigma.pool"=sigma.pool,
                           "sigma.mle"=sigma.mle,
                           "hclu.rcm"=hclu.rcm,
                           "hclu.pool"=hclu.pool,
                           "hclu.mle"=hclu.mle)
    
  }
  return(results)
}

## Calculate different measures of similarity between true and estimated
## sigma matrices given a results object from the function above
test.sigmas <- function(real.hclu, real.sigma, results_list){
  
  # Convert true clusters to dendrogram
  real.dendro <- as.dendrogram(real.hclu)
  
  # Convert estimated clustering to dendrogram
  rcm.dendro  <- as.dendrogram(results_list$hclu.rcm)
  mle.dendro  <- as.dendrogram(results_list$hclu.mle)
  pool.dendro <- as.dendrogram(results_list$hclu.pool)
  
  test.results <- c("nu"=results_list$nu,
                    "n"=results_list$n,
                    "rcm.copheno"=cor_cophenetic(real.dendro, rcm.dendro),
                    "mle.copheno"=cor_cophenetic(real.dendro, mle.dendro),
                    "pool.copheno"=cor_cophenetic(real.dendro, pool.dendro),
                    "rcm.kl"=divergence.kl(real.sigma, results_list$sigma.rcm),
                    "mle.kl"=divergence.kl(real.sigma, results_list$sigma.mle),
                    "pool.kl"=divergence.kl(real.sigma, results_list$sigma.pool))
  
  
  return(test.results)
}

## Make a neat table with mean and confidence interval for the mean
## given a a dataframe from test.sigmas()
makeTable <- function(results){
  means <- round(aggregate(results[,3:ncol(results)], list("n"=results[,2], "nu"=results[,1]), mean),2)
  ci <- round(aggregate(results[,3:ncol(results)], list("n"=results[,2], "nu"=results[,1]), ciMean),2)
  ci.low <- sapply(ci[,3:ncol(results)], function(x) x[, 1, drop = FALSE])
  ci.high <- sapply(ci[,3:ncol(results)], function(x) x[, 2, drop = FALSE])
  
  mean_ci <- matrix(paste(as.matrix(means[,3:ncol(results)]),
                          " (", ci.low, ";", ci.high, ")", sep=""), nrow=nrow(means))
  mean_ci <- data.frame(means[,1:2], mean_ci)
  names(mean_ci) <- names(means)
  return(mean_ci)
}

## Calculate P-value given a list of fitted models from permuted data
get.TestPValue <- function(the.list, the.object) {
  n <- sum(sapply(the.list, "[[", "nu") < the.object$nu) + 1
  N <- length(the.list) + 1
  return(n / N)
}


# Simulate data from model and fit models to permuted data
test.rcm.h0_par <- function(k = 4, n = 10, ns = rep(n, k), nsims = 10, 
                            p = 10, nu = 15, Psi, nCores=3, nPerm=500, e = 1e-2, ...) {
  stopifnot(nu > p - 1)
  
  if (missing(Psi)) {
    rho <- 0.5  # Compound symmetry matrix
    std <- 1
    Psi <- matrix(rho*std^2, p, p) + diag(rep((1 - rho)*std^2, p))
  }
  
  results <- list()
  
  for(nsim in 1:nsims){
    cat(paste("Simulation", nsim, "of", nsims, "with nu=", nu, ",n=", n, "\n"))
    # Sample random covariance matrices and data
    if (nu == Inf) {
      sigmas <- replicate(k, Psi)
    } else {
      sigmas <- rinvwishart(n = k, psi = Psi, nu = nu)
    }
    samples <- list()
    for(i in 1:k){
      samples[[i]] <- mvrnorm(n = n, mu = rep(0,p), Sigma=sigmas[,,i])
    }
    
    # Fit the model
    samples.ns  <- sapply(samples, nrow)
    samples.S   <- lapply(samples, correlateR::scatter)
    nu.i <- sum(samples.ns) + ncol(samples.S[[1]]) + 1
    psi <- nu.i*correlateR:::pool(samples.S, samples.ns)
    
    samples.rcm <- fit.rcm(S = samples.S, ns = samples.ns, verbose = FALSE,
                           Psi.init = psi, nu.init = nu.i, eps = 0.01,
                           max.ite = 2500)
    
    
    # Do the test for homogeneitiy
    dat <- do.call(rbind, samples)
    dat <- data.frame(dat)
    class.lab <- rep(1:k, each=n)
    
    homo.test <- function(tmp){
      library(correlateR)
      s.class.lab <- sample(class.lab)
      dat.sub <- split(dat, s.class.lab)
      dat.ns  <- sapply(dat.sub, nrow)
      dat.S   <- lapply(dat.sub, correlateR::scatter)
      nu <- sum(dat.ns) + ncol(dat.S[[1]]) + 1
      psi <- nu*correlateR:::pool(dat.S, dat.ns)
      
      rand.rcm.i <- fit.rcm(S = dat.S, ns = dat.ns, verbose = FALSE,
                            Psi.init = psi, nu.init = nu, eps = 0.01,
                            max.ite = 2500)
      
      return(rand.rcm.i)
    }
    
    cl <- makeCluster(nCores)
    ## make this reproducible
    clusterSetRNGStream(cl, 123)
    clusterExport(cl, varlist = c("dat", "class.lab"))
    homogeneity.rcm <- parLapply(cl, vector("list", nPerm), homo.test)
    stopCluster(cl)
    
    perm.nu <- sapply(homogeneity.rcm, "[[", "nu")
    p.value <- get.TestPValue(homogeneity.rcm, samples.rcm)
    results[[nsim]] <- list("k"=k,
                            "n"=n,
                            "nu"=nu,
                            "p"=p,
                            "sigmas"=sigmas, 
                            "samples"=samples,
                            "rcm"=samples.rcm,
                            "perm.nu",
                            "p.value"=p.value)
  }
  return(results)
}


## aux functions for survival_analysis
get.cis <- function(coxph) {
  x95 <- summary(coxph, conf.int = c(0.95))
  x99 <- summary(coxph, conf.int = c(0.99))
  
  return(data.frame("x95" = x95$conf.int, "x99" = x99$conf.int))
}


plot.cis <- function(dat, ...) {
  
  plot(1, type = "n",
       xlim = range(dat),
       ylim = c(0, nrow(dat) + 1),
       xlab = "", ylab = "", axes = FALSE, log = "x",
       xaxs = "r", ...)
  
  h <- 0.2
  col <- gsub("^ME", "", rownames(dat))
  axis(1, at = axTicks(3), label = formatC(axTicks(3)))
  axis(2, at = seq_along(col), label = cleanName(col), las = 2, tick = FALSE,
       pos = min(dat))
  
  
  n.seq <- 1:nrow(dat)
  rect(dat[,5], n.seq - h, dat[,6], n.seq + h, col = alp(col, 0.5), border = NA)
  rect(dat[,2], n.seq - h, dat[,3], n.seq + h)
  rect(dat[,5], n.seq - h, dat[,6], n.seq + h)
  segments(x0 = dat[,1], y0 = n.seq- 1.5*h, y1 = n.seq + 1.5*h, lwd = 2)
  segments(x0 = 1, y0 = 0, y1 = 6, col = "darkgrey", lty = 2)
  
}

plot.cis_bw <- function(dat, ...) { # Plot confidence interval i B/W
  
  plot(1, type = "n",
       xlim = range(dat),
       ylim = c(0, nrow(dat) + 1),
       xlab = "", ylab = "", axes = FALSE, log = "x",
       xaxs = "r", ...)
  
  h <- 0.2
  col <- gsub("^ME", "", rownames(dat))
  axis(1, at = axTicks(3), label = formatC(axTicks(3)))
  axis(2, at = seq_along(col), label = cleanName(col), las = 2, tick = FALSE,
       pos = min(dat))
  
  
  n.seq <- 1:nrow(dat)
  rect(dat[,5], n.seq - h, dat[,6], n.seq + h, col = alp("lightgrey", 0.5), border = NA)
  rect(dat[,2], n.seq - h, dat[,3], n.seq + h)
  rect(dat[,5], n.seq - h, dat[,6], n.seq + h)
  segments(x0 = dat[,1], y0 = n.seq- 1.5*h, y1 = n.seq + 1.5*h, lwd = 2)
  segments(x0 = 1, y0 = 0, y1 = 6, col = "darkgrey", lty = 2)
  
}

cleanName <- function(x) {
  return(gsub("[0-9]+|dark|medium|light", "", x))
}

## ---- auxiliary_functions ----
# Borrowed from Martin Morgan
# http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
factory <- function(fun) {
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(fun(...), error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(res, warn=warn, err=err)
  }
}




##################################
### Analysis of simulated data ###
##################################

## -- Define "True" Sigma matrix -- ##
sigma.diag.block <- matrix(0.5, 10, 10) + diag(rep(0.5, 10))
sigmaTrue <- as.matrix(bdiag(sigma.diag.block,
                             sigma.diag.block,
                             sigma.diag.block,
                             sigma.diag.block))

## Set 0.3 correlation between block 1&2 and 3&4, and 0.1 for the rest
sigmaTrue[11:20,1:10] = sigmaTrue[1:10,11:20] = 0.3
sigmaTrue[31:40,21:30] = sigmaTrue[21:30,31:40] = 0.3
sigmaTrue[which(sigmaTrue==0)] <- 0.1


## -- Create Sigma matrix from IDRC data -- ##
load("~/GitHub/RCM/studies.RData")
load("~/GitHub/RCM/gep.ensg.RData")
studies <- studies[studies$Study != "Celllines", ]
gep <- gep[grep("Celllines", names(gep), invert = TRUE)]

## Find top 100 genes for IDRC based on variance and subset
idrc.var <- apply(exprs(gep$GEPIDRC.ensg),1,var)
idrc.var.top <- names(sort(idrc.var, decreasing = TRUE))[1:100]
idrc.top <- exprs(gep$GEPIDRC.ensg)[idrc.var.top,]

## Build new "True" sigma matrix based on scatter matrix from IDRC
idrc.sigma.true <- cov2cor(correlateR::scatter(t(idrc.top)))

## -- Plot heatmaps for the two Sigma matrices ##
groupColors = c(rep("red", 10),
                rep("blue", 10),
                rep("green", 10),
                rep("purple", 10))

hmColors <- colorRampPalette(c("blue","white","red"))

jpeg("FigureS1a.jpg", height=7, width=7, units = "in", res = 200)
heatmap.2(sigmaTrue, # Figure S1a
          dendrogram = "column",
          trace = "none",
          dist = my.dist,
          hclustfun = my.hclust,
          labCol = NA,
          labRow = NA,
          breaks = seq(-1,1, length.out = 100),
          col= hmColors,
          ColSideColors = groupColors,
          main="Simulated Sigma matrix")
dev.off()

jpeg("FigureS1b.jpg", height=7, width=7, units = "in", res = 200)
heatmap.2(idrc.sigma.true, # Figure S1b
          distfun = my.dist,
          hclustfun = my.hclust,
          dendrogram = "column",
          trace = "none",
          labCol = NA,
          labRow = NA,
          col= hmColors,
          breaks = seq(-1,1, length.out = 100),
          main=paste("Sigma matrix from IDRC top100"))
dev.off()

## Make clustes from Sigma matrices
hclu.real <- flashClust(my.dist(sigmaTrue), method="ward")
hclu.idrc.top <- flashClust(my.dist(idrc.sigma.true), method="ward")

## -- Simulate data from Sigma matrices and fit models -- ##

## Sim data
if(!exists("results.clustering") || recompute){
  set.seed(42)
  results.clustering <- list()
  counter=1
  for(n in c(20,30,50,100,500,1000)){
    for(nuSim in c(50,100,1000,10000)){ 
      results.clustering[[counter]] <- test.rcm.simData(100, sigmaTrue, nuSim=nuSim, n=n, k=3)
      counter <- counter+1
    }
  }
  resave(results.clustering, file = "saved.RData")
}

if(!exists("results.test.sigmas") || recompute){
  results.test.sigmas <- lapply(results.clustering, sapply, 
                                function(x) test.sigmas(real.hclu=hclu.real, real.sigma=sigmaTrue, x))
  results.test.sigmas <- t(do.call(cbind, results.test.sigmas))
  resave(results.test.sigmas, file = "saved.RData")
}

## Table with only EM and Pool
results.test.sigmas.table <- makeTable(results.test.sigmas)
results.test.sigmas.table <- results.test.sigmas.table[,-c(4,6)] #Remove MLE results for copheno
names(results.test.sigmas.table) <- c("$n_i$", "$\\nu$", "EM","Pool", "EM", "Pool")

caption <- 'Mean cophenetic correlation and 
            Kullback-Leibler divergence with $95\\%$ confidence,
            for estimated  vs true network for 
            different values of $\\nu$ and $n_i$ using the EM or Pool method'

table1 <- latex(results.test.sigmas.table, file = "table1.tex",
                title = "Clustering results",
                caption = caption,
                size = "tiny",
                label ="tab:results.clustering",
                rowname=NULL,
                landscape = FALSE,
                cgroup = c("", "Cophenetic Correlation", "Kullback-Leibler divergence"),
                n.cgroup = c(2,2,2))

## Table with full results
results.test.sigmas.table <- makeTable(results.test.sigmas)
names(results.test.sigmas.table) <- c("$n_i$", "$\\nu$", "EM","MLE","Pool", "EM","MLE", "Pool")

caption <- 'Mean cophenetic correlation and 
Kullback-Leibler divergence with $95\\%$ confidence,
for estimated  vs true network for 
different values of $\\nu$ and $n_i$ using the EM, MLE or Pool method'

table1 <- latex(results.test.sigmas.table, file = "tableS2.tex",
                title = "Clustering results",
                caption = caption,
                size = "tiny",
                label ="tab:results.clustering.full",
                rowname=NULL,
                landscape = TRUE,
                cgroup = c("", "Cophenetic Correlation", "Kullback-Leibler divergence"),
                n.cgroup = c(2,3,3))


## IDRC data
if(!exists("results.clustering.idrc") || recompute){
  set.seed(42)
  results.clustering.idrc <- list()
  counter=1
  for(n in c(40,50,70,90, 150,500,1000)){
    for(nuSim in c(150,200,1000,10000)){ 
      results.clustering.idrc[[counter]] <- test.rcm.simData(100, idrc.sigma.true, nuSim=nuSim, n=n, k=3)
      counter <- counter+1
    }
  }
  resave(results.clustering.idrc, file = "saved.RData")
}

if(!exists("results.idrc.test.sigmas") || recompute){
  results.idrc.test.sigmas <- lapply(results.clustering.idrc, sapply, 
                                     function(x) test.sigmas(real.hclu=hclu.idrc.top, real.sigma=idrc.sigma.true, x))
  results.idrc.test.sigmas <- t(do.call(cbind, results.idrc.test.sigmas))
  resave(results.idrc.test.sigmas, file = "saved.RData")
}

results.idrc.test.sigmas.table <- makeTable(results.idrc.test.sigmas)
results.idrc.test.sigmas.table <- results.idrc.test.sigmas.table #Remove MLE results
names(results.idrc.test.sigmas.table) <- c("$n_i$", "$\\nu$", "EM","MLE","Pool", "EM","MLE", "Pool")


caption <- 'Simulation results based on IDRC data. Mean cophenetic correlation and 
            Kullback-Leibler divergence  with $95\\%$ confidence,
            for estimated vs true network for 
            different values of $\\nu$ and $n_i$ using the EM, MLE or Pool method.'

tableS1 <- latex(results.idrc.test.sigmas.table,
                file = "tableS1.tex",
                title = "Clustering results",
                caption = caption,
                size = "tiny",
                landscape = TRUE,
                label = "tab:results.clustering.idrc",
                rowname=NULL,
                cgroup = c("", "Cophenetic Correlation", "Kullback-Leibler divergence"),
                n.cgroup = c(2,3,3))

## Example Tanglegrams
exIndex <- 13
simIndex <- 2

## get Values for example data
exNu <- results.clustering[[exIndex]][[simIndex]]$nu
exN  <- results.clustering[[exIndex]][[simIndex]]$n
copheno.rcm  = subset(data.frame(results.test.sigmas),nu==exNu & n ==exN)$rcm.copheno[simIndex]
copheno.pool = subset(data.frame(results.test.sigmas),nu==exNu & n ==exN)$pool.copheno[simIndex]

jpeg("FigureS2A.jpg", height=7/2, width=7, units = "in", res = 200)
tanglegram(hclu.real, 
           results.clustering[[exIndex]][[simIndex]]$hclu.rcm, 
           main_left ="True", 
           main_right = "EM",
           sort=T,
           cex_main = 1,
           color_lines=groupColors, 
           main=paste("Simulation: nu=", exNu, ", n=", exN,
                      ", Cophenetic correlation =", round(copheno.rcm,2)))
dev.off()

jpeg("FigureS2B.jpg", height=7/2, width=7, units = "in", res = 200)
tanglegram(hclu.real, 
           results.clustering[[exIndex]][[simIndex]]$hclu.pool, 
           main_left ="True", 
           main_right = "Pool",
           sort=T,
           cex_main = 1,
           color_lines=groupColors, 
           main=paste("Simulation: nu=", exNu, ", n=", exN,
                      ", Cophenetic correlation =", round(copheno.pool,2)))
dev.off()

########################
### Computation time ###
########################

### Set parameters for testing runtime
par.t <- list(k = 3,
              nu = 300,
              p = c(100, 150, 200, 250),
              n.sims = 10,
              n.obs = 210)


### Test runtime
if(!exists("df.t") || recompute){
  set.seed(722040)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.t$p)) {
    tmp <- foreach(j = seq_len(par.t$n.sims)) %do% {
      with(par.t, test.rcm(k = k, n = n.obs, p = p[i], nu = nu))
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.t$p), "done after",
        (proc.time()-st)[3] %/% 60, "mins.\n")
  }
  df.t <- as.data.frame(t(sapply(res, SSEs)))
  resave(df.t, file = "saved.RData")
}

### Extract runtime results
tm.elapsed <-
  aggregate(cbind(time.em.elapsed,  time.pool.elapsed, time.mle.elapsed) ~
              n + nu + k + p, FUN = mean,
            data = df.t)

### Plot runtimes
cols.fig1 <- c("black", "darkgray", "lightgray")
jpeg("Figure1.jpg", height=7, width=10, units = "in", res = 200)
plot(tm.elapsed$p, tm.elapsed$time.em.elapsed,
     type = "b",
     col = cols.fig1[1],
     xlim =c(100, 250),
     ylab = "Computation time (s)",
     xlab = "p",
     axes = FALSE,
     pch = 15)
axis(1)
axis(2)
lines(tm.elapsed$p, tm.elapsed$time.pool.elapsed, type = "b",
      col = cols.fig1[2], pch = 16)
lines(tm.elapsed$p, tm.elapsed$time.mle.elapsed, type = "b",
      col = cols.fig1[3], pch = 17)

legend("topleft", legend = c("EM", "pool", "Approx. MLE"),
       lty = 1:4, pch = c(15, 16, 17), lwd = 2, bty = "n",
       col = cols.fig1, inset = 0.05)

lgnd <- c(substitute(k == x,    list(x = unique(tm.elapsed$k))),
          substitute(nu == x,   list(x = unique(tm.elapsed$nu))),
          substitute(n[i] == x, list(x = unique(tm.elapsed$n))))
legend("left", inset = 0.01, bty = "n", horiz = FALSE,
       legend = as.expression(lgnd))
dev.off()

#####################
### Test P-values ###
#####################
if(!exists("results.test.p.values") || recompute){
  par_l <- list(k=3, p=20, n=c(10,20,50), nu=c(30,300,Inf))
  par_h <- list(k=3, p=50, n=c(20,50), nu=c(115,1150,Inf))
  
  set.seed(42)
  counter <- 1
  results.test.p.values <- list()
  
  for (i in seq_along(par_l$nu)){
    for(j in seq_along(par_l$n)){
      results.test.p.values[[counter]] <- with(par_l, test.rcm.h0_par(k = k,
                                                                      n = n[j],
                                                                      p = p,
                                                                      nu = nu[i],
                                                                      nsims=100,
                                                                      nCores=3))
      counter <- counter+1
    }
  }
  
  for (i in seq_along(par_h$nu)){
    for(j in seq_along(par_h$n)){
      results.test.p.values[[counter]] <- with(par_h, test.rcm.h0_par(k = k,
                                                                      n = n[j],
                                                                      p = p,
                                                                      nu = nu[i],
                                                                      nsims=100,
                                                                      nCores=3))
      counter <- counter+1
    }
  }
  resave(results.test.p.values, file = "saved.RData")
}

results.test.p.values2=list()
for(i in 1:length(results.test.p.values)){
  p.value = sapply(results.test.p.values[[i]],"[[","p.value")
  nu = sapply(results.test.p.values[[i]],"[[","nu")
  p = sapply(results.test.p.values[[i]],"[[","p")
  n = sapply(results.test.p.values[[i]],"[[","n")
  nuhat = unlist(sapply(results.test.p.values[[i]],"[[","rcm")["nu",])
  results.test.p.values2[[i]]<-data.frame(p,n,nu,nuhat,p.value)
  rm(p.value,nu,n,p,nuhat)
}

results.test.p.values.plot <- do.call(rbind, results.test.p.values2)
results.test.p.values.plot$scenario <- with(results.test.p.values.plot, paste0("p=",p, ", nu=", nu, ", n_i=", n))

p.value.plot <- ggplot(results.test.p.values.plot, aes(x=scenario, y=p.value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),text = element_text(size=15))

jpeg("FigureS3.jpg", height=7, width=10, units = "in", res = 200)
  print(p.value.plot)
dev.off()



###############################
### ----  DLBCL DATA ---- #####
###############################

### Set parameters for DLBCL analysis
dlbcl.dims <- sapply(gep, dim)
dlbcl.par <- list(top.n = 300,
                  linkage = "ward",
                  go.alpha.level = 0.01,
                  ontology = "MF",
                  minModuleSize = 20,
                  n.modules = 5,
                  threshold = NA)

### Subset to top 300 genes
vars      <- sapply(gep, function(x) rowSds(exprs(x))*(ncol(x) - 1))
var.pool  <- rowSums(vars)/(sum(sapply(gep, ncol)) - length(gep))
use.genes <- names(sort(var.pool, decreasing = TRUE)[seq_len(dlbcl.par$top.n)])
gep.sub <- lapply(gep, function(x) exprs(x)[use.genes, ])

### PCA plot using all datasets
gep.all <- do.call(cbind, gep.sub)
gep.all.pca <- prcomp(t(gep.all))
gep.all.pca.names <- rep(names(gep.sub), sapply(gep.sub,ncol))
gep.all.pca.names <- gsub("GEP", "", gep.all.pca.names)
gep.all.pca.names <- gsub(".ensg", "", gep.all.pca.names)

### Figure 4
g <- ggbiplot(gep.all.pca, obs.scale = 1, var.scale = 1,
              choices=1:2,
              ellipse = TRUE,
              circle = TRUE, 
              var.axes = FALSE,
              groups=gep.all.pca.names) +
              labs(color="Dataset") #+
#              theme(aspect.ratio = 1)
              
jpeg("FigureS4.jpg", height=5, width=7, units = "in", res = 200)
print(g)
dev.off()

### Fit RCM model three times with different initial
### values to investigate convergence
if (!exists("dlbcl.rcm") || recompute) {
  dlbcl.ns  <- sapply(gep.sub, ncol)
  dlbcl.S   <- lapply(gep.sub, function(x) correlateR::scatter(t(x)))
  
  nu <- sum(dlbcl.ns) + ncol(dlbcl.S[[1]]) + 1
  psi <- nu*correlateR:::pool(dlbcl.S, dlbcl.ns)
  
  dlbcl.time <- system.time({
    dlbcl.trace <- capture.output({
      dlbcl.rcm <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                           Psi.init = psi, nu.init = nu, eps = 0.01,
                           max.ite = 1500)
    })
  })
  dimnames(dlbcl.rcm$Psi) <- dimnames(dlbcl.S[[1]])
  dlbcl.rcm$time <- dlbcl.time
  dlbcl.rcm$trace <- dlbcl.trace
  resave(dlbcl.rcm, file = "saved.RData")
  
  # Alternative fit 1
  dlbcl.time.alt <- system.time({
    dlbcl.trace.alt <- capture.output({
      dlbcl.rcm.alt <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                               Psi.init = psi/2, nu.init = ceil(nu/2),
                               eps = 0.01, max.ite = 1500)
    })
  })
  dimnames(dlbcl.rcm.alt$Psi) <- dimnames(dlbcl.S[[1]])
  dlbcl.rcm.alt$time <- dlbcl.time.alt
  dlbcl.rcm.alt$trace <- dlbcl.trace.alt
  resave(dlbcl.rcm.alt, file = "saved.RData")
  
  # Alternative fit 2
  tmp <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                 Psi.init = diag(diag(psi)), nu.init = 1000,
                 eps = 1, max.ite = 2, method = "pooled") # Very crude start
  dlbcl.time.alt2 <- system.time({
    dlbcl.trace.alt2 <- capture.output({
      dlbcl.rcm.alt2 <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                                Psi.init = tmp$Psi,
                                nu.init = tmp$nu,
                                eps = 0.01, max.ite = 5000)
    })
  })
  dimnames(dlbcl.rcm.alt2$Psi) <- dimnames(dlbcl.S[[1]])
  dlbcl.rcm.alt2$time <- dlbcl.time.alt2
  dlbcl.rcm.alt2$trace <- dlbcl.trace.alt2
  resave(dlbcl.rcm.alt2, file = "saved.RData")
  
  ### Fit pool model
  dlbcl.pool <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                        Psi.init = psi, nu.init = nu, eps = 0.01,
                        max.ite = 1500, method="pool")
    resave(dlbcl.pool, file = "saved.RData")
}

# Check equality
diff.psi <- max(dlbcl.rcm.alt$Psi - dlbcl.rcm$Psi)/max(dlbcl.rcm.alt$Psi, dlbcl.rcm$Psi)
diff.nu <- (dlbcl.rcm.alt$nu - dlbcl.rcm$nu)/dlbcl.rcm$nu
diff.psi2 <- max(dlbcl.rcm.alt$Psi - dlbcl.rcm$Psi)/max(dlbcl.rcm.alt2$Psi, dlbcl.rcm$Psi)
diff.nu2 <- (dlbcl.rcm.alt2$nu - dlbcl.rcm$nu)/dlbcl.rcm$nu
stopifnot(diff.psi < 1e-3)
stopifnot(diff.nu < 1e-3)
stopifnot(diff.psi2 < 1e-3)
stopifnot(diff.nu2 < 1e-3)

# Make convergence plot
get.ll.trace <- function(x) {
  as.numeric(gsub(" L = ","",sapply(strsplit(x, "\\|"),"[",2)))
}
ll.tr <-get.ll.trace(dlbcl.rcm$trace)
ll.tr.a <- get.ll.trace(dlbcl.rcm.alt$trace)
ll.tr.a2 <- get.ll.trace(dlbcl.rcm.alt2$trace)


jpeg("Figure2.jpg", width = 7, height = 7/2, units = "in", res = 200)
{
  fig1.cols <- c("black", "darkgrey", "lightgrey")
  par(mar = c(4, 4.5, 2, 0) + 0.1, mgp = c(2.5, 1, 0))
  plot(ll.tr, type = "l", ylab = "", xlab = "Iteration", lwd = 2,
       xlim = c(0,2000), axes = FALSE, col = fig1.cols[1])
  axis(1)
  axis(2, las = 2)
  title(ylab = "log-likelihood", line = 3.6)
  grid()
  lines(ll.tr.a, lty = 2, col = fig1.cols[2], lwd = 2)
  lines(ll.tr.a2, lty = 3, col = fig1.cols[3], lwd = 2)
  legend("bottomright", col = fig1.cols, lty = 1:3,
         bty = "n", y.intersp = 1.5,
         legend = paste0("Fit ", 1:3, " (",
                         round(c(dlbcl.rcm$time[3],
                                 dlbcl.rcm.alt$time[3],
                                 dlbcl.rcm.alt2$time[3])/60, 0), " min)"),
         inset = 0.1, lwd = 2)
  axis(3, at = length(ll.tr),    col = fig1.cols[1])
  axis(3, at = length(ll.tr.a),  col = fig1.cols[2])
  axis(3, at = length(ll.tr.a2), col = fig1.cols[3])
}
dev.off()

### Get ICC for DLBCL data
with(dlbcl.rcm, ICC(nu, nrow(Psi)))

### Do the test for homogeneitiy
set.seed(100)
if (!exists("homogeneity.rcm") || recompute) {
  dat <- do.call(rbind, lapply(gep.sub, function(x) t(x)))
  dat <- data.frame(dat)
  class.lab <- rep(names(gep), sapply(gep.sub, ncol))
  
  homo.test <- function(tmp){
    library(correlateR)
    s.class.lab <- sample(class.lab)
    dat.sub <- split(dat, s.class.lab)
    dat.ns  <- sapply(dat.sub, nrow)
    dat.S   <- lapply(dat.sub, correlateR::scatter)
    nu <- sum(dat.ns) + ncol(dat.S[[1]]) + 1
    psi <- nu*correlateR:::pool(dat.S, dat.ns)
    
    rand.rcm.i <- fit.rcm(S = dat.S, ns = dat.ns, verbose = FALSE,
                          Psi.init = psi, nu.init = nu, eps = 0.01,
                          max.ite = 1500)
    
    return(rand.rcm.i)
  }
  
  cl <- makeCluster(3)
  ## make this reproducible
  clusterSetRNGStream(cl, 123)
  clusterExport(cl, varlist = c("dat", "class.lab"))
  homogeneity.rcm <- parLapply(cl, vector("list", 500), homo.test)
  stopCluster(cl)
  
  resave(homogeneity.rcm, file = "saved.RData")
}

dlbcl.p.value <- get.TestPValue(homogeneity.rcm, dlbcl.rcm)


## Conversion between ENSG and HGNC
if(!exists("gene.info") || recompute){
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


### Retrieve Sigma and Correlation matrices and calculate Hclusts for RCM
dlbcl.rcm.sigma  <- with(dlbcl.rcm, Psi2Sigma(Psi,nu))
dlbcl.rcm.cor <- cov2cor(dlbcl.rcm.sigma)
dlbcl.rcm.hclu   <- flashClust(my.dist(dlbcl.rcm.cor), method=dlbcl.par$linkage)


### Define colors for clusters
num2col <- c("gray32",
             "darkolivegreen3",
             "mediumorchid3",
             "lightskyblue3",
             "coral3")

### Make clusters from RCM model
dlbcl.cut <- cutree(dlbcl.rcm.hclu, k = dlbcl.par$n.modules)
dlbcl.modules <- num2col[dlbcl.cut]
names(dlbcl.modules) <- dlbcl.rcm.hclu$labels
tangle.colors <- dlbcl.modules[order.dendrogram(as.dendrogram(dlbcl.rcm.hclu))]


## Heatmap for RCM
jpeg("Figure3A.jpg", width = 7, height = 7, units="in", res = 200)
heatmap.2(cov2cor(dlbcl.rcm.sigma),
          distfun = my.dist,
          col = hmColors,
          breaks = seq(-1,1, length.out = 100),
          hclustfun = my.hclust,
          dendrogram = "column",
          trace = "none",
          ColSideColors = dlbcl.modules,
          labRow = NA,
          labCol = NA,
          main="EM")
dev.off()


## iGraph for RCM
# PANEL GRAPH layout
layout.custom <- function(graph,...) {
  l <- layout.circle(graph)
  layout_with_fr(graph, niter = 5000,
                 weights = abs(E(graph)$weight),
                 coords = l)
}


scaleToLayout <- function(x) {
  return(2*apply(x, 2, function(x) (x - min(x))/max(x - min(x))) - 1)
}

# Colors for igraph
n.col <- 256
keybreaks <- seq(min(dlbcl.rcm.cor), max(dlbcl.rcm.cor), length.out = n.col)

keycols <- rep("White", n.col - 1)
tmp <- rle(diff(sign(keybreaks)))
neg <- tmp$lengths[1]
pos <- tmp$lengths[3]
keycols[1:neg] <- colorRampPalette(c("Blue", "White"))(neg)
keycols[(neg+2):(neg+pos+1)] <- colorRampPalette(c("White", "Red"))(pos)


# Make graph
dlbcl.g <- graph.adjacency(abs(dlbcl.rcm.cor), mode = "undirected",
                           weighted = TRUE, diag = FALSE)
V(dlbcl.g)$name  <- map2hugo(V(dlbcl.g)$name)
V(dlbcl.g)$color <- dlbcl.modules
E(dlbcl.g)$color <- keycols[cut(E(dlbcl.g)$weight, keybreaks)]

w <- E(dlbcl.g)$weight
dlbcl.par$threshold <- quantile(abs(w), prob = 0.90)

## Define labels for color legend
num2names <- c("Antigen & receptor binding",
               "Fatty acid binding & peptidase activity",
               "Inflamation & immuneresponse. Lipid metabolisme.",
               "Extracellular matrix structure & Growthfactor.",
               "Cytoskeleton organization")
names(num2names) <- num2col


# Plot graph
jpeg("Figure3B.jpg", width = 7, height = 7, units="in", res = 200)
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,4))

tab <- table(dlbcl.modules)
go.func <- num2names[match(names(tab), names(num2names))]
stopifnot(names(tab) == names(go.func))

par(xaxs = "i", yaxs = "i", mar = c(0, 0, 0, 0) + 0.1)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(0.7, 1, bty = "n", col = "black", pt.bg = names(tab),
       cex = 1.2, pch = 22, pt.cex = 3,
       title = "Modules", legend = capitalize(cleanName(names(tab))),
       xjust = 0.5, yjust = 0.5)
legend(.8, 1, bty = "n", col = "black", cex = 1.2, xjust = 0.5, yjust = 0.5,
       title = "Size", legend =  sprintf("(%d)", tab))
legend(1.07, 1, bty = "n", col = "black", cex = 1.2, xjust = 0.5, yjust = 0.5,
       title = "Suggested GO function", legend =  rep("", 5))
legend(1.12, 1, bty = "n", col = "black", cex = 1.0,
       title = "", xjust = 0.5, yjust = 0.5,
       legend =  go.func, y.intersp = 1.25)

# Plot Graph
tmp.g <- delete.edges(dlbcl.g, which(abs(w) < dlbcl.par$threshold))
V(tmp.g)$size <- 3
E(tmp.g)$width <- 1.3
V(tmp.g)$frame.color <- V(tmp.g)$color
l <- layout.custom(tmp.g)
plot(tmp.g, layout = l, vertex.label = "")
points(scaleToLayout(l), pch = 16, cex = 0.3)
dev.off()


### Survival analysis with EM modules
### Survival analysis for clusters
load("~/GitHub/RCM/metadata.RData")
the.module <- num2col[2] # = "darkolivegreen3"

eigenvars <- list()

figure4 <- "Figure4.jpg"
jpeg(figure4, height = 1.2*7, width = 2*7, units = "in", res = 200)
{
  par(mar = c(4, 4, 1, 0) + 0.1, oma = c(0,0,0,0), xpd = TRUE, cex = 1.2,
      mgp = c(3, 1, 0)*0.75)
  layout(cbind(1:2,3:4), heights = c(1,2))
  
  for (j in 1:2) {
    
    meta <- switch(j, metadataLLMPPCHOP, metadataLLMPPRCHOP)
    rownames(meta) <- as.character(meta$GEO.ID)
    expr <- switch(j,
                   (gep.sub$GEPLLMPPCHOP.ensg)[names(dlbcl.modules), ],
                   (gep.sub$GEPLLMPPRCHOP.ensg)[names(dlbcl.modules), ])
    meta <- meta[colnames(expr), ] # Reorder
    
    # Check order
    stopifnot(rownames(meta) == colnames(expr))
    
    res <- moduleEigengenes(t(expr), dlbcl.modules)
    eigenvars[[j]] <- res$varExplained
    eg <- res$eigengenes
    col <- gsub("^ME", "", colnames(eg))
    eg <- as.data.frame(lapply(eg, function(x) x/sd(x))) # Standardize
    
    #   # Multivariate
    cph.fits <- coxph(meta$OS ~ ., data = eg, x = TRUE, y = TRUE)
    dats <- get.cis(cph.fits)
    dats <- dats[, -c(2,6)]
    plot.cis_bw(dats) #use B/W
    #plot.cis(dats)
    
    for (i in seq_len(ncol(eg))) {
      if (col[i] == the.module) {
        eg.i <- eg[, paste0("ME", col[i])]
        eg.high <- eg.i >= mean(eg.i)
        sfit.h <- survfit(meta$OS ~ 1, subset = eg.high)
        sfit.l <- survfit(meta$OS ~ 1, subset = !eg.high)
        
        plot(sfit.h,
             conf.int = TRUE,
             main = "",
             col = "grey", #col[i], B/W
             lwd = 2,
             lty = c(1,2,2),
             axes = FALSE,
             xlab = "years",
             ylab = "Survival proportion")
        lines(sfit.l,
              lwd = 2,
              lty = c(1,3,3))
        axis(1)
        axis(2)
        legend("bottomleft", bty = "n", lwd = 2, lty = 2:3,
               legend = c("High eigengene", "Low eigengene"),
               #col = c(col[i], "black"), horiz = FALSE)
               col = c("grey", "black"), horiz = FALSE) # B/W
        legend("topright", legend = paste(cleanName(col[i]), "eigengene"),
               #bty = "n", text.col = col[i]) #B/W
               bty = "n", text.col = "grey") 
      }
    }
    title(switch(j, "GSE10846 CHOP", "GSE10846 R-CHOP"), xpd = TRUE)
  }
}
dev.off()


#####################################
### Top genes in DLBCL RCM modules ##
#####################################


# Genes in modules
mod.genes <- list()
for (col in unique(dlbcl.modules)) {
  mod.genes[[col]] <- names(dlbcl.modules[dlbcl.modules == col])
}
  
stopifnot(all(names(mod.genes) == names(num2names)))

# Order by rowSums
dlbcl.cor.sub <- lapply(mod.genes, function(ensg) dlbcl.rcm.cor[ensg, ensg])
dlbcl.cor.sub <- lapply(dlbcl.cor.sub, function(x) {
  o <- order(rowSums(abs(x)), decreasing = TRUE)
  return(x[o,o])
})

tmp <- lapply(dlbcl.cor.sub, rownames)
n.tmp <- unlist(lapply(tmp, length))
m <- max(table(dlbcl.modules))  # Get largest module

# Construct table
convert_to_hugo <- function(x){
   x <- map2hugo(x)
   x <- x[grep("ENSG",x, invert=TRUE)]
  return(unname(c(x, rep(NA, m - length(x)))))
}

dlbcl.mod.tab.genes <- sapply(tmp, convert_to_hugo)

# First letter capitalized
cgroup <- capitalize(cleanName(names(mod.genes)))
# cgroup <- paste0(cgroup, num2names)

colnames(dlbcl.mod.tab.genes) <- # Number of gens in each module
  paste0("n = ", n.tmp)


nn <- 40
seq_tmp <- seq_len(min(nn, nrow(dlbcl.mod.tab.genes)))
mod.tab.table <- latex(dlbcl.mod.tab.genes[seq_tmp, ],
                       cgroup = cgroup,
                       n.cgroup = rep(1,5),
                       size = "tiny",
                       caption = paste("The identified modules, their sizes, and",
                                       "member genes. The genes are sorted decreasingly",
                                       "by their intra-module connectivity (sum of the incident",
                                       "edge weights). Only the top", nn, "genes are shown."),
                       caption.loc = "bottom",
                       label = "tab:top.genes.em",
                       landscape = FALSE,
                       file = "table.top.genes.em.tex")




########################################
### Enrichment analysis for clusters ###
########################################

go.genes <- gsub("_at", "", use.genes) # Consider only top 300 genes

dlbcl.rcm.gprofile <- list()
for(module in unique(dlbcl.modules)){
  module.tmp <- sub("_at", "", names(dlbcl.modules[dlbcl.modules==module]))
  dlbcl.rcm.gprofile[[module]] <- gprofiler(module.tmp, custom_bg = go.genes)#, correction_method = "fdr")  
}

pickcolumns <- function(x){
  x <- x[,c("term.id",
            "domain",
            "term.name",
            "p.value",
            "query.size",
            "overlap.size",
            "recall",
            "precision")]
  return(x)
}

dlbcl.rcm.gprofile.table <- do.call(rbind,lapply(dlbcl.rcm.gprofile, pickcolumns))
dlbcl.rcm.gprofile.table$term.id <- paste0("$",dlbcl.rcm.gprofile.table$term.id,"$")
dlbcl.rcm.gprofile.table$term.name <- substr(dlbcl.rcm.gprofile.table$term.name,1,30)
names(dlbcl.rcm.gprofile.table) <-c("Term ID", "Domain", "Term", "P", "N", "O", "Recall","Precision")



gp.col.rcm <- gsub("(^|[[:space:]])([[:alpha:]])",
                   "\\1\\U\\2", gsub("\\.[0-9]+$", "", row.names(dlbcl.rcm.gprofile.table)),
                   perl = TRUE)



go.table.rcm.gp <- latex(dlbcl.rcm.gprofile.table[,-c(1,7,8)],
                         rgroup = capitalize(cleanName(unique(tolower(gp.col.rcm)))),
                         rowname = dlbcl.rcm.gprofile.table[,1],
                         label  = "tab:enrichment.em",
                         n.rgroup = table(gp.col.rcm)[unique(gp.col.rcm)],
                         title = "Term ID",
                         size = "tiny",
                         caption = paste0("The significant terms for the gene enrichment ,
                                          analysis of the DLBCL EM method modules. Number of genes ,
                                          in each term (N), and the overlap to module (O)."),
                         lines.page = 80,
                         longtable = TRUE,
                         center = "centering",
                         file = "table.enrichment.em.tex")


##########################################
### Repeat analysis with pool clusters ###
##########################################


#### Retrieve Sigma and Correlation matrices and calculate Hclusts for Pool
dlbcl.pool.sigma <- with(dlbcl.pool, Psi2Sigma(Psi,nu))
rownames(dlbcl.pool.sigma) <- colnames(dlbcl.pool.sigma) <- colnames(dlbcl.rcm.sigma)
dlbcl.pool.cor <- cov2cor(dlbcl.pool.sigma)
dlbcl.pool.hclu  <- flashClust(my.dist(dlbcl.pool.cor), method=dlbcl.par$linkage)

### Colors for different ordering of clusters in pool model
num2col_pool <- c("gray32",
                  "darkolivegreen3",
                  "lightskyblue3",
                  "mediumorchid3",
                  "coral3")

### Make clusters from POOL model
dlbcl.pool.cut <- cutree(dlbcl.pool.hclu, k = dlbcl.par$n.modules)
dlbcl.pool.modules <- num2col[dlbcl.pool.cut]
names(dlbcl.pool.modules) <- dlbcl.pool.hclu$labels


table(dlbcl.modules, dlbcl.pool.modules)


## Heatmap for Pool
jpeg("FigureS5A.jpg", width = 7, height = 7, units="in", res = 200)
heatmap.2(cov2cor(dlbcl.pool.sigma),
          distfun = my.dist,
          col = hmColors,
          breaks = seq(-1,1, length.out = 100),
          hclustfun = my.hclust,
          dendrogram = "column",
          trace = "none",
          ColSideColors = dlbcl.pool.modules,
          labRow = NA,
          labCol = NA,
          main="Pool")
dev.off()


### iGraph for pool
# Colors for igraph
n.col <- 256
keybreaks <- seq(min(dlbcl.pool.cor), max(dlbcl.pool.cor), length.out = n.col)

keycols <- rep("White", n.col - 1)
tmp <- rle(diff(sign(keybreaks)))
neg <- tmp$lengths[1]
pos <- tmp$lengths[3]
keycols[1:neg] <- colorRampPalette(c("Blue", "White"))(neg)
keycols[(neg+2):(neg+pos+1)] <- colorRampPalette(c("White", "Red"))(pos)

dlbcl.pool.g <- graph.adjacency(abs(dlbcl.pool.cor), mode = "undirected",
                                weighted = TRUE, diag = FALSE)
V(dlbcl.pool.g)$name  <- map2hugo(V(dlbcl.pool.g)$name)
V(dlbcl.pool.g)$color <- dlbcl.pool.modules
E(dlbcl.pool.g)$color <- keycols[cut(E(dlbcl.pool.g)$weight, keybreaks)]


w.pool <- E(dlbcl.pool.g)$weight
dlbcl.par$threshold.pool <- quantile(abs(w.pool), prob = 0.90)

jpeg("FigureS5B.jpg", width = 7, height = 7, units="in", res = 200)
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
       widths=c(1,1), heights=c(1,4))

# PANEL MODULE OVERVIEW
tab <- table(dlbcl.pool.modules)
go.func <- num2names[match(names(tab), names(num2names))]
stopifnot(names(tab) == names(go.func))

par(xaxs = "i", yaxs = "i", mar = c(0, 0, 0, 0) + 0.1)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend(0.7, 1, bty = "n", col = "black", pt.bg = names(tab),
       cex = 1.2, pch = 22, pt.cex = 3,
       title = "Modules", legend = capitalize(cleanName(names(tab))),
       xjust = 0.5, yjust = 0.5)
legend(.8, 1, bty = "n", col = "black", cex = 1.2, xjust = 0.5, yjust = 0.5,
       title = "Size", legend =  sprintf("(%d)", tab))
legend(1.07, 1, bty = "n", col = "black", cex = 1.2, xjust = 0.5, yjust = 0.5,
       title = "Suggested GO function", legend =  rep("", 5))
legend(1.12, 1, bty = "n", col = "black", cex = 1.0,
       title = "", xjust = 0.5, yjust = 0.5,
       legend =  go.func, y.intersp = 1.25)

# Plot graph
tmp.g <- delete.edges(dlbcl.pool.g, which(abs(w.pool) < dlbcl.par$threshold.pool))
V(tmp.g)$size <- 3
E(tmp.g)$width <- 1.3
V(tmp.g)$frame.color <- V(tmp.g)$color
l <- layout.custom(tmp.g)
plot(tmp.g, layout = l, vertex.label = "")
points(scaleToLayout(l), pch = 16, cex = 0.3)
dev.off()


### Tanglegram of RCM vs Pool
dlbcl.rcm.dend <- as.dendrogram(dlbcl.rcm.hclu)
labels_colors(dlbcl.rcm.dend) <- dlbcl.modules[order.dendrogram(dlbcl.rcm.dend)]

dlbcl.pool.dend <- as.dendrogram(dlbcl.pool.hclu)
labels_colors(dlbcl.pool.dend) <- dlbcl.pool.modules[order.dendrogram(dlbcl.pool.dend)]

copheno_rcm_pool <- cor_cophenetic(dlbcl.rcm.dend, dlbcl.pool.dend)
jpeg("FigureS5C.jpg", width = 7, height = 7/2, units="in", res = 200)
tanglegram(dlbcl.rcm.dend, dlbcl.pool.dend,
           sort=F, 
           main_left="EM", 
           main_right="Pool", 
           lwd=2,
           cex_main = 1.2,
           color_lines=tangle.colors,
           lab.cex=0.5,
           edge.lwd=2,
           # k_labels = 3,
           main=paste("Cophenetic correlation =", round(copheno_rcm_pool,2))
)
dev.off()

eigenvarspool <- list()
figureS6 <- "FigureS6.jpg"
jpeg(figureS6, height = 1.2*7, width = 2*7, units = "in", res = 200)
{
  par(mar = c(4, 4, 1, 0) + 0.1, oma = c(0,0,0,0), xpd = TRUE, cex = 1.2,
      mgp = c(3, 1, 0)*0.75)
  layout(cbind(1:2,3:4), heights = c(1,2))
  
  for (j in 1:2) {
    
    meta <- switch(j, metadataLLMPPCHOP, metadataLLMPPRCHOP)
    rownames(meta) <- as.character(meta$GEO.ID)
    expr <- switch(j,
                   (gep.sub$GEPLLMPPCHOP.ensg)[names(dlbcl.pool.modules), ],
                   (gep.sub$GEPLLMPPRCHOP.ensg)[names(dlbcl.pool.modules), ])
    meta <- meta[colnames(expr), ] # Reorder
    
    # Check order
    stopifnot(rownames(meta) == colnames(expr))
    
    res <- moduleEigengenes(t(expr), dlbcl.pool.modules)
    eigenvarspool[[j]] <- res$varExplained
    eg <- res$eigengenes
    col <- gsub("^ME", "", colnames(eg))
    eg <- as.data.frame(lapply(eg, function(x) x/sd(x))) # Standardize
    
    #   # Multivariate
    cph.fits <- coxph(meta$OS ~ ., data = eg, x = TRUE, y = TRUE)
    dats <- get.cis(cph.fits)
    dats <- dats[, -c(2,6)]
    plot.cis(dats)
    
    for (i in seq_len(ncol(eg))) {
      if (col[i] == the.module) {
        eg.i <- eg[, paste0("ME", col[i])]
        eg.high <- eg.i >= mean(eg.i)
        sfit.h <- survfit(meta$OS ~ 1, subset = eg.high)
        sfit.l <- survfit(meta$OS ~ 1, subset = !eg.high)
        
        plot(sfit.h,
             conf.int = TRUE,
             main = "",
             col = col[i],
             lwd = 2,
             lty = c(1,2,2),
             axes = FALSE,
             xlab = "years",
             ylab = "Survival proportion")
        lines(sfit.l,
              lwd = 2,
              lty = c(1,3,3))
        axis(1)
        axis(2)
        legend("bottomleft", bty = "n", lwd = 2, lty = 2:3,
               legend = c("High eigengene", "Low eigengene"),
               col = c(col[i], "black"), horiz = FALSE)
        legend("topright", legend = paste(cleanName(col[i]), "eigengene"),
               bty = "n", text.col = col[i])
      }
    }
    title(switch(j, "GSE10846 CHOP", "GSE10846 R-CHOP"), xpd = TRUE)
  }
}
dev.off()



### Top genes in DLBCL Pool modules
# Genes in modules
mod.genes.pool <- list()
for (col in unique(dlbcl.pool.modules)) {
  mod.genes.pool[[col]] <- names(dlbcl.pool.modules[dlbcl.pool.modules == col])
}

stopifnot(all(names(mod.genes.pool) == names(num2names)))

# Order by rowSums
dlbcl.cor.sub <- lapply(mod.genes.pool, function(ensg) dlbcl.pool.cor[ensg, ensg])
dlbcl.cor.sub <- lapply(dlbcl.cor.sub, function(x) {
  o <- order(rowSums(abs(x)), decreasing = TRUE)
  return(x[o,o])
})

tmp <- lapply(dlbcl.cor.sub, rownames)
n.tmp <- unlist(lapply(tmp, length))
m <- max(table(dlbcl.pool.modules))  # Get largest module

# Construct table
convert_to_hugo <- function(x){
  x <- map2hugo(x)
  x <- x[grep("ENSG",x, invert=TRUE)]
  return(unname(c(x, rep(NA, m - length(x)))))
}

dlbcl.mod.tab.genes <- sapply(tmp, convert_to_hugo)

# First letter capitalized
cgroup <- capitalize(cleanName(names(mod.genes)))
# cgroup <- paste0(cgroup, num2names)

colnames(dlbcl.mod.tab.genes) <- # Number of gens in each module
  paste0("n = ", n.tmp)

nn <- 40
seq_tmp <- seq_len(min(nn, nrow(dlbcl.mod.tab.genes)))
mod.tab.table <- latex(dlbcl.mod.tab.genes[seq_tmp, ],
                       cgroup = cgroup,
                       size = "tiny",
                       caption = paste("The identified modules from the Pool method, their sizes, and",
                                       "member genes. The genes are sorted decreasingly",
                                       "by their intra-module connectivity (sum of the incident",
                                       "edge weights). Only the top", nn, "genes are shown."),
                       caption.loc = "bottom",
                       label = "tab:top.genes.pool",
                       landscape = FALSE,
                       file = "table.top.genes.pool.tex")

### Enrichment for pool modules
dlbcl.pool.gprofile <- list()
for(module in unique(dlbcl.pool.modules)){
  module.tmp <- sub("_at", "", names(dlbcl.pool.modules[dlbcl.pool.modules==module]))
  dlbcl.pool.gprofile[[module]] <- gprofiler(module.tmp, custom_bg = go.genes)#, correction_method = "fdr")  
}


dlbcl.pool.gprofile.table <- do.call(rbind,lapply(dlbcl.pool.gprofile, pickcolumns))
dlbcl.pool.gprofile.table$term.id <- paste0("$",dlbcl.pool.gprofile.table$term.id,"$")
dlbcl.pool.gprofile.table$term.name <- substr(dlbcl.pool.gprofile.table$term.name,1,30)
names(dlbcl.pool.gprofile.table) <-c("Term ID", "Domain", "Term", "P", "N", "O", "Recall","Precision")

gp.col.pool <- gsub("(^|[[:space:]])([[:alpha:]])",
                    "\\1\\U\\2", gsub("\\.[0-9]+$", "", row.names(dlbcl.pool.gprofile.table)),
                    perl = TRUE)


go.table.pool.gp <- latex(dlbcl.pool.gprofile.table[,-c(1,7,8)],
                          rgroup = capitalize(cleanName(unique(tolower(gp.col.pool)))),
                          rowname = dlbcl.pool.gprofile.table[,1],
                          label  = "tab:enrichment.pool",
                          n.rgroup = table(gp.col.pool)[unique(gp.col.pool)],
                          title = "Term ID",
                          size = "tiny",
                          caption = paste0("The significant terms for the gene enrichment ,
                                          analysis of the DLBCL Pool method modules. Number of genes ,
                                          in each term (N), and the overlap to module (O)."),
                          lines.page = 80,
                          longtable = TRUE,
                          center = "centering",
                          file = "table.enrichment.pool.tex")



