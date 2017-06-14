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
                    "rcm.kl"=divergence.kl(real.sigma, cov2cor(results_list$sigma.rcm)),
                    "mle.kl"=divergence.kl(real.sigma, cov2cor(results_list$sigma.mle)),
                    "pool.kl"=divergence.kl(real.sigma, cov2cor(results_list$sigma.pool)))
  
  
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

results.test.sigmas.table <- makeTable(results.test.sigmas)
results.test.sigmas.table <- results.test.sigmas.table[,-c(4,7)] #Remove MLE results
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
                cgroup = c("", "Cophenetic Correlation", "Kullback-Leibler divergence"),
                n.cgroup = c(2,2,2))


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
results.idrc.test.sigmas.table <- results.idrc.test.sigmas.table[,-c(4,7)] #Remove MLE results
names(results.idrc.test.sigmas.table) <- c("$n_i$", "$\\nu$", "EM","Pool", "EM", "Pool")


caption <- 'Simulation results based on IDRC data. Mean cophenetic correlation and 
            Kullback-Leibler divergence  with $95\\%$ confidence,
            for estimated vs true network for 
            different values of $\\nu$ and $n_i$ using the EM or Pool method'

tableS1 <- latex(results.idrc.test.sigmas.table,
                file = "tableS1.tex",
                title = "Clustering results",
                caption = caption,
                size = "tiny",
                label = "tab:results.clustering.idrc",
                rowname=NULL,
                cgroup = c("", "Cophenetic Correlation", "Kullback-Leibler divergence"),
                n.cgroup = c(2,2,2))

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
                  n.modules = 3,
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





### Retrieve Sigma matrices and calculate Hclusts
dlbcl.rcm.sigma  <- with(dlbcl.rcm, Psi2Sigma(Psi,nu))
dlbcl.pool.sigma <- with(dlbcl.pool, Psi2Sigma(Psi,nu))

rownames(dlbcl.pool.sigma) <- colnames(dlbcl.pool.sigma) <- colnames(dlbcl.rcm.sigma)

dlbcl.rcm.hclu   <- flashClust(my.dist(cov2cor(dlbcl.rcm.sigma)), method=dlbcl.par$linkage)
dlbcl.pool.hclu  <- flashClust(my.dist(cov2cor(dlbcl.pool.sigma)), method=dlbcl.par$linkage)

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

### Make clusters from POOL model
dlbcl.pool.cut <- cutree(dlbcl.pool.hclu, k = dlbcl.par$n.modules)
dlbcl.pool.modules <- num2col[dlbcl.pool.cut]
names(dlbcl.pool.modules) <- dlbcl.pool.hclu$labels

table(dlbcl.modules, dlbcl.pool.modules)

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

## Heatmap for Pool
jpeg("Figure3B.jpg", width = 7, height = 7, units="in", res = 200)
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

### Tanglegram of RCM vs Pool
dlbcl.rcm.dend <- as.dendrogram(dlbcl.rcm.hclu)
labels_colors(dlbcl.rcm.dend) <- dlbcl.modules[order.dendrogram(dlbcl.rcm.dend)]

dlbcl.pool.dend <- as.dendrogram(dlbcl.pool.hclu)
labels_colors(dlbcl.pool.dend) <- dlbcl.pool.modules[order.dendrogram(dlbcl.pool.dend)]

copheno_rcm_pool <- cor_cophenetic(dlbcl.rcm.dend, dlbcl.pool.dend)
jpeg("Figure3C.jpg", width = 7, height = 7/2, units="in", res = 200)
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

### Survival analysis for clusters
load("~/GitHub/RCM/metadata.RData")
the.module <- num2col[2] # = "darkolivegreen3"



plot_surv_fig <- function(meta, gep, sampleID){
  par(mar = c(4, 4, 1, 0) + 0.1, oma = c(0,0,0,0), xpd = TRUE, cex = 1.2,
      mgp = c(3, 1, 0)*0.75)
  layout(cbind(1:2,3:4), heights = c(1,2))
  
  rownames(meta) <- as.character(sampleID)
  
  for (j in 1:2) {
    
    modules <- switch(j, dlbcl.modules, dlbcl.pool.modules)
    expr <- gep[names(modules),]
    meta <- meta[colnames(expr), ] # Reorder
    
    # Check order
    stopifnot(rownames(meta) == colnames(expr))
    
    res <- moduleEigengenes(t(expr), modules)
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
    title(switch(j, "EM", "Pool"), xpd = TRUE)
  }
}

jpeg("Figure4.jpg", height = 1.2*7, width = 2*7, units = "in", res = 200)
  plot_surv_fig(metadataLLMPPCHOP, gep.sub$GEPLLMPPCHOP.ensg, metadataLLMPPCHOP$GEO.ID)
dev.off()

jpeg("FigureS5.jpg", height = 1.2*7, width = 2*7, units = "in", res = 200)
  plot_surv_fig(metadataLLMPPRCHOP, gep.sub$GEPLLMPPRCHOP.ensg, metadataLLMPPRCHOP$GEO.ID)
dev.off()

#plot_surv_fig(metadataIDRC, gep.sub$GEPIDRC.ensg, gsub(".CEL", "", metadataIDRC$Array.Data.File))
#plot_surv_fig(metadataCHEPRETRO, gep.sub$GEPCHEPRETRO.ensg, gsub(".CEL", "", metadataCHEPRETRO$file))


### TOPGO analysis
all.genes <- gsub("_at$", "", names(sort(var.pool, decreasing = TRUE)))

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

go.genes <- gsub("_at", "", use.genes) # Consider only top 300 genes

### Get HGNC symbol for genes within each group for RCM
if(!exists("dlbcl.go.analysis") || recompute){
  dlbcl.go.analysis <- list()
  dlbcl.module.genes <- list()
  for (col in unique(dlbcl.modules)) {
    dlbcl.module.genes[[col]] <- mod.genes <-
      gsub("_at$", "", names(dlbcl.modules[dlbcl.modules == col]))
    geneList <- factor(as.integer(go.genes %in% mod.genes))
    names(geneList) <- go.genes
    GOdata <- new("topGOdata", ontology = dlbcl.par$ontology,
                  allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = gid2go)
    result.fisher <-
      runTest(GOdata, algorithm = "classic", statistic = "fisher")
    mod.res <- GenTable(GOdata, classicFisher = result.fisher, topNodes = 100)
    dlbcl.go.analysis[[col]] <-
      mod.res %>% filter(classicFisher <= dlbcl.par$go.alpha.level)
  }
  resave(dlbcl.go.analysis, file = "saved.RData")
  dlbcl.module.genes <- lapply(dlbcl.module.genes, map2hugo)
  resave(dlbcl.module.genes, file = "saved.RData")
}


### Get HGNC symbol for genes within each group for POOL
if(!exists("dlbcl.pool.go.analysis") || recompute){
  dlbcl.pool.go.analysis <- list()
  dlbcl.pool.module.genes <- list()
  for (col in unique(dlbcl.pool.modules)) {
    dlbcl.pool.module.genes[[col]] <- mod.genes <-
      gsub("_at$", "", names(dlbcl.pool.modules[dlbcl.pool.modules == col]))
    geneList <- factor(as.integer(go.genes %in% mod.genes))
    names(geneList) <- go.genes
    GOdata <- new("topGOdata", ontology = dlbcl.par$ontology,
                  allGenes = geneList,
                  annot = annFUN.gene2GO, gene2GO = gid2go)
    result.fisher <-
      runTest(GOdata, algorithm = "classic", statistic = "fisher")
    mod.res <- GenTable(GOdata, classicFisher = result.fisher, topNodes = 100)
    dlbcl.pool.go.analysis[[col]] <-
      mod.res %>% filter(classicFisher <= dlbcl.par$go.alpha.level)
  }
  resave(dlbcl.pool.go.analysis, file = "saved.RData")
  dlbcl.pool.module.genes <- lapply(dlbcl.pool.module.genes, map2hugo)
  resave(dlbcl.pool.module.genes, file = "saved.RData")
}

### Make table with GO terms

# EM
tmp <- dlbcl.go.analysis
names(tmp) <- cleanName(names(tmp))
go.table <- do.call(rbind, tmp)
go.col <- gsub("(^|[[:space:]])([[:alpha:]])",
               "\\1\\U\\2", gsub("\\.[0-9]+$", "", row.names(go.table)),
               perl = TRUE)
colnames(go.table)[3:6] <- c("$N$", "$O$", "$E$", "$P$-val")
rownames(go.table) <- NULL


# Pool
tmp <- dlbcl.pool.go.analysis
names(tmp) <- cleanName(names(tmp))
go.table2 <- do.call(rbind, tmp)
go.col2 <- gsub("(^|[[:space:]])([[:alpha:]])",
               "\\1\\U\\2", gsub("\\.[0-9]+$", "", row.names(go.table2)),
               perl = TRUE)
colnames(go.table2)[3:6] <- c("$N$", "$O$", "$E$", "$P$-val")
rownames(go.table2) <- NULL

# Combined table
go.table3 <- latex(rbind(go.table[,-c(1,3)], go.table2[,-c(1,3)]),
                   rgroup = c(paste("EM", unique(go.col), sep=" - "), paste("Pool", unique(go.col2), sep=" - ")), 
                   rowname = c(go.table[, 1], go.table2[, 1]),
                   label  = "tab:results.go",
                   n.rgroup = c(table(go.col)[unique(go.col)],table(go.col2)[unique(go.col2)]),
                   title = "GO ID",
                   size = "tiny",
                   caption = paste0("The significant GO terms in the modules from the EM and Pool method at ",
                                    "$\\alpha$-level ", dlbcl.par$go.alpha.level, "."),
                   lines.page = 80,
                   center = "centering",
                   file = "table3.tex")


# Combined table Wide
#go.table3 <- merge(go.table[,c(1,2,6)], go.table2[,c(1,2,6)], by=0, all.y=T)
#go.col <- gsub("(^|[[:space:]])([[:alpha:]])",
#               "\\1\\U\\2", gsub("\\.[0-9]+$", "", go.table3$Row.names),
#               perl = TRUE)
#go.table3 <- go.table3[order(go.col, go.table3$classicFisher.y),]
#go.table3 <- go.table3[,-1]
#names(go.table3) <- c("GO ID", "Term", "P", "GO ID", "Term", "P")
#go.table3[is.na(go.table3)] <- " "
#rownames(go.table3) <- NULL

#go.table <- latex(go.table3[,-1],
#                  rowname = go.table3[,1],
#                  rgroup = unique(go.col),
#                  label  = "tab:results.go",
#                  cgroup = c("EM", "Pool"),
#                  n.cgroup = c(2,3),
#                  n.rgroup = table(go.col)[unique(go.col)],
#                  title = "GO ID",
#                  size = "tiny",
#                  caption = paste0("The significant GO terms in the EM- and Pool-modules at ",
#                                   "$\\alpha$-level ", dlbcl.par$go.alpha.level, "."),
#                  #longtable = TRUE,
 #                 lines.page = 80,
#                  center = "centering",
#                  file = "table3.tex")



### Correlate eigenGenes to Genes in modules for CHOP dataset
eg.cor <- list()
for (j in 1:2) {
  modules <- switch(j, dlbcl.modules, dlbcl.pool.modules)
  name <- switch(j, "RCM", "POOL")
  expr <- gep.sub$GEPLLMPPCHOP.ensg[names(modules),]
  
  res <- moduleEigengenes(t(expr), modules)
  eg <- res$eigengenes
  
  for(module in num2col[1:3]){
    modeuleGenes <- which(modules == module)
    expr2 <- cbind(eg[, paste0("ME", module)], t(expr[modeuleGenes,]))
    eg.cor.j <- cor(expr2)[-1,1]
    names(eg.cor.j) <- map2hugo(gsub("_at", "",names(eg.cor.j)))
    eg.cor[[name]][[cleanName(module)]] <- sort(eg.cor.j, decreasing = TRUE)
  }
}

## Make top Gene tables
rcm.eigen.cors <- array(0, dim=c(20,6))
pool.eigen.cors <- array(0, dim=c(20,6))
counter <- 1
for(module in cleanName(num2col[1:3])){
  rcm.eigen.cors[,counter] <- names(eg.cor[["RCM"]][[module]])[1:20]
  rcm.eigen.cors[,counter+1] <- round(eg.cor[["RCM"]][[module]][1:20],2)
  
  pool.eigen.cors[,counter] <- names(eg.cor[["POOL"]][[module]])[1:20]
  pool.eigen.cors[,counter+1] <- round(eg.cor[["POOL"]][[module]][1:20],2)
  
  counter <- counter +2
}

colnames(rcm.eigen.cors) <- c("Gene", "Cor","Gene", "Cor","Gene", "Cor")
colnames(pool.eigen.cors) <- c("Gene", "Cor","Gene", "Cor","Gene", "Cor")

## Comparison table for olivegreen module
olive.eg.cors <- array(0, dim=c(50,4))
olive.eg.cors[,1] <- names(eg.cor[["RCM"]][["olivegreen"]])[1:50]
olive.eg.cors[,2] <- round(eg.cor[["RCM"]][["olivegreen"]][1:50],2)
olive.eg.cors[,3] <- names(eg.cor[["POOL"]][["olivegreen"]])[1:50]
olive.eg.cors[,4] <- round(eg.cor[["POOL"]][["olivegreen"]][1:50],2)

colnames(olive.eg.cors) <- c("Gene", "Cor","Gene", "Cor")


## Make Latex tables
rcm.eg.table <- latex(rcm.eigen.cors,
                  label  = "tab:rcm.eg.cors",
                  cgroup = names(eg.cor$RCM),
                  n.cgroup = c(2,2,2),
                  size = "tiny",
                  caption = paste("The top 20 genes in each module for the EM method and",
                                   "correlation to the eigenGene for the CHOP dataset"),
                  longtable = TRUE,
                  lines.page = 80,
                  center = "centering",
                  file = "table.rcm.eg.tex")


pool.eg.table <- latex(pool.eigen.cors,
                      label  = "tab:pool.eg.cors",
                      cgroup = names(eg.cor$POOL),
                      n.cgroup = c(2,2,2),
                      size = "tiny",
                      caption = paste("The top 20 genes in each module for the Pool method and",
                                       "correlation to the eigenGene for the CHOP dataset"),
                      longtable = TRUE,
                      lines.page = 80,
                      center = "centering",
                      file = "table.pool.eg.tex")

olive.eg.table <- latex(olive.eg.cors,
                      label  = "tab:olive.eg.cors",
                      cgroup = c("EM", "Pool"),
                      n.cgroup = c(2,2),
                      size = "tiny",
                      caption = paste("The top 50 genes in the olive module for the EM vs Pool method and",
                                       "correlation to the eigenGene for the CHOP dataset"),
                      longtable = TRUE,
                      lines.page = 80,
                      center = "centering",
                      file = "table.olive.eg.tex")


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


################################################
### Repeat survival analysis with 5 clusters ###
################################################

### Define colors for clusters
num2col <- c("gray32",
             "darkolivegreen3",
             "mediumorchid3",
             "lightskyblue3",
             "coral3")


### Colors for different ordering of clusters in pool model
num2col_pool <- c("gray32",
                  "darkolivegreen3",
                  "lightskyblue3",
                  "mediumorchid3",
                  "coral3")



### Make clusters from RCM model
dlbcl.cut <- cutree(dlbcl.rcm.hclu, k = 5)
dlbcl.modules <- num2col[dlbcl.cut]
names(dlbcl.modules) <- dlbcl.rcm.hclu$labels
tangle.colors <- dlbcl.modules[order.dendrogram(as.dendrogram(dlbcl.rcm.hclu))]

### Make clusters from POOL model
dlbcl.pool.cut <- cutree(dlbcl.pool.hclu, k = 5)
dlbcl.pool.modules <- num2col_pool[dlbcl.pool.cut]
names(dlbcl.pool.modules) <- dlbcl.pool.hclu$labels

table(dlbcl.modules, dlbcl.pool.modules)

## Heatmap for RCM
jpeg("FigureS6A.jpg", width = 7, height = 7, units="in", res = 200)
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

## Heatmap for Pool
jpeg("FigureS6B.jpg", width = 7, height = 7, units="in", res = 200)
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

### Tanglegram of RCM vs Pool
dlbcl.rcm.dend <- as.dendrogram(dlbcl.rcm.hclu)
labels_colors(dlbcl.rcm.dend) <- dlbcl.modules[order.dendrogram(dlbcl.rcm.dend)]

dlbcl.pool.dend <- as.dendrogram(dlbcl.pool.hclu)
labels_colors(dlbcl.pool.dend) <- dlbcl.pool.modules[order.dendrogram(dlbcl.pool.dend)]

copheno_rcm_pool <- cor_cophenetic(dlbcl.rcm.dend, dlbcl.pool.dend)
jpeg("FigureS6C.jpg", width = 7, height = 7/2, units="in", res = 200)
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

jpeg("FigureS7.jpg", height = 1.2*7, width = 2*7, units = "in", res = 200)
plot_surv_fig(metadataLLMPPCHOP, gep.sub$GEPLLMPPCHOP.ensg, metadataLLMPPCHOP$GEO.ID)
dev.off()

jpeg("FigureS8.jpg", height = 1.2*7, width = 2*7, units = "in", res = 200)
plot_surv_fig(metadataLLMPPRCHOP, gep.sub$GEPLLMPPRCHOP.ensg, metadataLLMPPRCHOP$GEO.ID)
dev.off()
