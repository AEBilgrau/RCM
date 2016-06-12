# setwd("~/Documents/PhD/paper-RCM/")
# setwd("~/GitHub/RCM/")

## ---- initialize_script ----
rm(list = ls())
set.seed(1)
recompute <- FALSE
library("Hmisc")
library("Bmisc")
library("correlateR")
library("GMCM")
library("MASS")
library("flashClust")
library("WGCNA")
library("affy")
library("biomaRt")
library("ape")
library("topGO")
library("igraph")
library("MCMCpack")
library("dplyr")

if (file.exists("saved.RData")) loaded <- load("saved.RData")

# Multicore support
library("foreach")
library("doMC")
registerDoMC(detectCores())
# library("doParallel") # Use this package on windows
# registerDoParallel(detectCores())

num2col <- c("gray32",
             "darkolivegreen3",
             "mediumorchid3",
             "lightskyblue3",
             "coral3")

num2names <- c("Antigen & receptor binding",
               "Fatty acid binding & peptidase activity",
               "Inflamation & immuneresponse. Lipid metabolisme.",
               "Extracellular matrix structure & Growthfactor.",
               "Cytoskeleton organization")
names(num2names) <- num2col

cleanName <- function(x) {
  return(gsub("[0-9]+|dark|medium|light", "", x))
}
#c("#C55F4B", "#98BE53", "#A463B5", "#94B9B5", "#614051")
#names(num2col) <- c("Coral", "Pistachio", "Orchid", "Viridian", "Eggplant")

## ---- auxiliary_functions ----

#
# HDA
#

# Density of the RCM model
drcm <- correlateR::drcm

# Simulate data
test.rcm <- function(k = 4, n = 10, ns = rep(n, k), p = 10, nu = 15, Psi, ...) {
  #stopifnot(nu > p - 1)

  if (missing(Psi)) {
    rho <- 0.5  # Compound symmetry matrix
    std <- 1
    Psi <- matrix(rho*std^2, p, p) + diag(rep((1 - rho)*std^2, p))
  }

  S <- createRCMData(ns = ns, psi = Psi, nu = nu)
  e <- 1e-2
  t_em   <- system.time(res.em   <- fit.rcm(S, ns, method = "EM",   eps=e, ...))
  t_pool <- system.time(res.pool <- fit.rcm(S, ns, method = "pool", eps=e, ...))
  t_mle  <- system.time(res.mle  <- fit.rcm(S, ns, method = "appr", eps=e, ...))
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
  sse.rcm.em   <- getSSE(getSigma(x$rcm.em),   expected, n = n)
  sse.rcm.mle  <- getSSE(getSigma(x$rcm.mle),  expected, n = n)
  sse.rcm.pool <- getSSE(getSigma(x$rcm.pool), expected, n = n)

  get <- c("nu", "iterations", "loglik")
  stopifnot(all(x$ns == x$ns[1]))
  return(c(n  = x$ns[1],
           k  = length(x$ns),
           nu = x$nu,
           p  = nrow(x$Psi),
           SSE.rcm.em   = sse.rcm.em,
           SSE.rcm.mle  = sse.rcm.mle,
           SSE.rcm.pool = sse.rcm.pool,
           em = unlist(x$rcm.em[get]),
           mle = unlist(x$rcm.mle[get]),
           pool = unlist(x$rcm.pool[get]),
           time = x$time))
}

to.df <- function(sim) {
  return(data.frame(sim$z, classes = factor(sim$K)))
}

misclassificationRisk <- function(x) {
  s <- sum(x)
  return((s - sum(diag(x)))/s)
}

misclassificationRiskSE <- function(phat, n) {
  return(sqrt(phat*(1 - phat)/n))
}

accuracy <- function(x) {
  return(sum(diag(x))/sum(x))
}


## ---- numerical_experiment ----
# As a function of n_i
par.ne <- list(k = 3,
               nu = 30,
               p = 20,
               n.sims = 1000,
               n.obs = ceil(c(7, seq(10, 40, length.out = 7))))

if (!exists("df.numerical") || recompute) {
  set.seed(987654321)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.ne$n.obs)) {
    tmp <- foreach(j = seq_len(par.ne$n.sims)) %do% {
      with(par.ne, test.rcm(k = k, n = n.obs[i], p = p, nu = nu))
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.ne$n.obs), "done after",
        (proc.time() - st)[3] %/% 60, "mins.\n")
  }
  df.numerical <- as.data.frame(t(sapply(res, SSEs)))
  rm(tmp)
  resave(df.numerical, file = "saved.RData")
}

# As a function of p
par.ne2 <- list(k = 3,
                nu = 150,
                p = ceil(seq(5, 40, length.out = 8)),
                n.sims = 1000,
                n.obs = 40)

if (!exists("df.numerical2") || recompute) {
  set.seed(1234567890)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.ne2$p)) {
    tmp <- foreach(j = seq_len(par.ne2$n.sims)) %dopar% {
      with(par.ne2, test.rcm(k = k, n = n.obs, p = p[i], nu = nu))
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.ne2$p), "done after",
        (proc.time()-st)[3] %/% 60, "mins.\n")
  }
  df.numerical2 <- as.data.frame(t(sapply(res, SSEs)))
  rm(tmp)
  resave(df.numerical2, file = "saved.RData")
}

# Time pr fit
par.ne3 <- list(k = 3,
                nu = 300,
                p = c(100, 150, 200, 250),
                n.sims = 10,
                n.obs = 210)

if (!exists("df.numerical3") || recompute) {
  set.seed(3842856)
  st <- proc.time()
  res <- list()
  for (i in seq_along(par.ne3$p)) {
    tmp <- foreach(j = seq_len(par.ne3$n.sims)) %do% {
      with(par.ne3, test.rcm(k = k, n = n.obs, p = p[i], nu = nu))
    }
    res <- c(res, tmp)
    cat("loop =", i, "of", length(par.ne3$p), "done after",
        (proc.time()-st)[3] %/% 60, "mins.\n")
  }
  df.numerical3 <- as.data.frame(t(sapply(res, SSEs)))
  resave(df.numerical3, file = "saved.RData")
}

## ---- end ----


## ---- numerical_experiment_plot ----
ci <- function(x) 1.96*sd(x)/sqrt(length(x))
df <- aggregate(cbind(SSE.rcm.em, SSE.rcm.mle, SSE.rcm.pool) ~
                  n + nu + k + p, FUN = mean, data = df.numerical)
df.mad <- aggregate(cbind(SSE.rcm.em, SSE.rcm.mle, SSE.rcm.pool) ~
                      n + nu + k + p, FUN = ci,
                    data = df.numerical)
tm.elapsed <-
  aggregate(cbind(time.em.elapsed,  time.pool.elapsed, time.mle.elapsed) ~
              n + nu + k + p, FUN = mean,
            data = df.numerical3)

figure1 <- "figure1.jpg"
jpeg(figure1, height=7/2, width=7, units = "in", res = 200)
{
  par(mfrow = c(1,2))
  plot(df$n, df$SSE.rcm.em,
       col = num2col[3],
       type = "b",
       axes = FALSE,
       xlab = expression(n = n[i]),
       ylim = range(df[,5:7]),
       ylab = "mean SSE",
       pch = 15,
       lty = 1,
       main = "")

  legend_expressions <-
    sapply(1:3, function(i) {
      as.expression(substitute(x == y,list(x = as.name(c("k", "nu", "p")[i]),
                                    y = unique(c(df$k, df$nu, df$p))[i])))
    })
  legend("bottomleft", inset = 0.01, bty = "n", horiz = TRUE,
         legend = legend_expressions)
  axis(1)
  axis(2)
  grid()

  lines(df$n+0.2, df$SSE.rcm.pool, col = num2col[4], type = "b", pch=16, lty=2, lwd=2)
  lines(df$n+0.4, df$SSE.rcm.mle, col = num2col[5], type = "b", pch=17, lty=3, lwd=2)
  legend("topright", legend = c("EM", "pool", "Approx. MLE"),
         lty = 1:4, pch = c(15, 16, 17), lwd = 2, bty = "n",
         col = num2col[c(3,4,5)], inset = 0.05)

  arrows(df$n, df$SSE.rcm.em-df.mad$SSE.rcm.em,
         df$n, df$SSE.rcm.em+df.mad$SSE.rcm.em,
         length=0.05, angle=90, code=3, col = num2col[3])
  arrows(df$n+0.2, df$SSE.rcm.pool-df.mad$SSE.rcm.pool,
         df$n+0.2, df$SSE.rcm.pool+df.mad$SSE.rcm.pool,
         length=0.05, angle=90, code=3, col = num2col[4])
  arrows(df$n+0.4, df$SSE.rcm.mle-df.mad$SSE.rcm.mle,
         df$n+0.4, df$SSE.rcm.mle+df.mad$SSE.rcm.mle,
         length=0.05, angle=90, code=3, col = num2col[5])



  # Panel 2
  plot(tm.elapsed$p, tm.elapsed$time.em.elapsed,
       type = "b",
       col = num2col[3],
       xlim =c(100, 250),
       ylab = "Computation time (s)",
       xlab = "p",
       axes = FALSE,
       pch = 16)
  grid()
  axis(1)
  axis(2)
  lines(tm.elapsed$p, tm.elapsed$time.pool.elapsed, type = "b",
        col = num2col[4], pch = 16)
  lines(tm.elapsed$p, tm.elapsed$time.mle.elapsed, type = "b",
        col = num2col[5], pch = 16)

  legend("topleft", legend = c("EM", "pool", "Approx. MLE"),
         lty = 1:4, pch = c(15, 16, 17), lwd = 2, bty = "n",
         col = num2col[c(3,4,5)], inset = 0.05)
}
dev.off()


## ---- end ----




################################################################################
#
# DLBCL analysis
#
################################################################################

## ---- dlbcl_analysis ----
the.module <- num2col[3] # = "mediumorchid3"

load("studies.RData")
load("gep.ensg.RData")
studies <- studies[studies$Study != "Celllines", ]
gep <- gep[grep("Celllines", names(gep), invert = TRUE)]

dlbcl.dims <- sapply(gep, dim)
dlbcl.par <- list(top.n = 300,
                  linkage = "ward",
                  go.alpha.level = 0.01,
                  ontology = "MF",
                  minModuleSize = 20,
                  n.modules = 5,
                  threshold = NA)

vars      <- sapply(gep, function(x) rowSds(exprs(x))*(ncol(x) - 1))
var.pool  <- rowSums(vars)/(sum(sapply(gep, ncol)) - length(gep))
use.genes <- names(sort(var.pool, decreasing = TRUE)[seq_len(dlbcl.par$top.n)])
gep.sub <- lapply(gep, function(x) exprs(x)[use.genes, ])

if (!exists("dlbcl.rcm") || !exists("var.pool") || recompute) {
  dlbcl.ns  <- sapply(gep.sub, ncol)
  dlbcl.S   <- lapply(gep.sub, function(x) correlateR::scatter(t(x)))
  nu <- sum(dlbcl.ns) + ncol(dlbcl.S[[1]]) + 1
  psi <- nu*correlateR:::pool(dlbcl.S, dlbcl.ns)
  dlbcl.rcm <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                       Psi.init = psi, nu.init = nu, eps = 0.01,
                       max.ite = 1500)
  dimnames(dlbcl.rcm$Psi) <- dimnames(dlbcl.S[[1]])
  resave(dlbcl.rcm, file = "saved.RData")
}

# Expected covariance matrix
dlbcl.exp <- with(dlbcl.rcm, Psi2Sigma(Psi, nu))
dlbcl.cor <- cov2cor(dlbcl.exp)




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

n.col <- 256
keybreaks <- seq(min(dlbcl.cor), max(dlbcl.cor), length.out = n.col)

keycols <- rep("White", n.col - 1)
tmp <- rle(diff(sign(keybreaks)))
neg <- tmp$lengths[1]
pos <- tmp$lengths[3]
keycols[1:neg] <- colorRampPalette(c("Blue", "White"))(neg)
keycols[(neg+2):(neg+pos+1)] <- colorRampPalette(c("White", "Red"))(pos)

# Clustering and tree cut
dlbcl.clust <- flashClust(as.dist(1 - abs(dlbcl.cor)),
                          method = dlbcl.par$linkage)
dlbcl.cut <- cutree(dlbcl.clust, k = dlbcl.par$n.modules)
dlbcl.modules <- num2col[dlbcl.cut]
names(dlbcl.modules) <- dlbcl.clust$labels


# Construct graph
dlbcl.g <- graph.adjacency(dlbcl.cor, mode = "undirected",
                           weighted = TRUE, diag = FALSE)

V(dlbcl.g)$name  <- map2hugo(V(dlbcl.g)$name)
V(dlbcl.g)$color <- dlbcl.modules
E(dlbcl.g)$color <- keycols[cut(E(dlbcl.g)$weight, keybreaks)]

# Phylo / dendrogram
phylo <- as.phylo(dlbcl.clust)
phylo$tip.label <- map2hugo(phylo$tip.label)
## ---- end ----


## ---- dlbcl_plot ----
plotColourMat <- function(x, add = FALSE, ...) {
  if (missing(x)) x <- matrix(sample(colors())[1:(6*10)], 6, 10)
  nr <- nrow(x)
  nc <- ncol(x)
  if (!add) {
    plot(1, type = "n", xlab = "", ylab = "", axes = FALSE,
         xlim = c(.5, nc + .5), ylim = c(.5, nr + .5))
  }
  xx <- rep(seq_len(nc), each = nr)
  yy <- rep(seq(nr, 1), nc)
  col <- c(x)
  rect(xx - 0.5, yy - 0.5, xx + 0.5, yy + 0.5, col = col, border = NA, ...)
  return(invisible(list(x = xx, y = yy, col = col)))
}

plotColorkey <- function(breaks, col, add = FALSE, ...) {
  stopifnot(length(col) + 1 == length(breaks))
  nb <- length(breaks)
  if (!add) {
    plot(breaks, xlim = range(breaks), ylim = c(0,1), type = "n", axes = FALSE,
         xlab ="", ylab = "")
  }
  rect(breaks[-nb], 0, breaks[-1], 1, col = col, border = NA, ...)
}

w <- E(dlbcl.g)$weight
dlbcl.par$threshold <- quantile(abs(w), prob = 0.95)
dlbcl_plot <- "figure2.jpg"
if (!file.exists(dlbcl_plot) || recompute) {
  jpeg(dlbcl_plot, width = 1800, height = 2640, res = 200)

  # LAYOUT
  lmat <- rbind(c(8,8,5,9), c(8,8,2,0), c(4,1,3,6), 7)
  lwid <- c(4, 1, 15, 15)
  lhei <- c(4, 1, 10, 20)
  layout(lmat, widths = c(4, 1, 15, 15), heights = c(4, 1, 10, 20))
  op <- par()
  par(mar = c(0,0,0,0)+0.1, xaxs = "i", yaxs = "i")

  # PANEL CLASSIFICATIONS
  o <- dlbcl.clust$order
  plotColourMat(cbind(dlbcl.modules[o]))
  plotColourMat(rbind(dlbcl.modules[o]))

  # HEATMAP
  nc <- ncol(dlbcl.cor)
  image(dlbcl.cor[o,o][, nc:1], axes = FALSE,
        breaks = keybreaks, col = keycols)

  # DENDROGRAMS
  tmp.dend <- as.dendrogram(dlbcl.clust)
  plot(rev(tmp.dend), leaflab = "none", horiz = TRUE)
  abline(v = 5.5, col = "grey", lty = 2)
  plot(tmp.dend, leaflab = "none", axes = FALSE)
  abline(h = 5.5, col = "grey", lty = 2)

  # PANEL GRAPH
  layout.custom <- function(graph,...) {
    l <- layout.circle(graph)
    layout.fruchterman.reingold(graph, niter = 5000,
                                area = vcount(graph)/2,
                                repulserad = vcount(graph)/2,
                                weights = abs(E(graph)$weight),
                                start = l, ...)
  }

  scaleToLayout <- function(x) {
    return(2*apply(x, 2, function(x) (x - min(x))/max(x - min(x))) - 1)
  }

  tmp.g <- delete.edges(dlbcl.g, which(abs(w) < dlbcl.par$threshold))
  V(tmp.g)$size <- 3
  E(tmp.g)$width <- 1.3
  V(tmp.g)$frame.color <- V(tmp.g)$color
  l <- layout.custom(tmp.g)
#   plot(1)
  plot(tmp.g, layout = l, vertex.label = "")
  points(scaleToLayout(l), pch = 16, cex = 0.3)


  # PANEL HEB
  # http://github.com/AEBilgrau/HierarchicalEdgeBundles
  # devtools::install_github("AEBilgrau/HierarchicalEdgeBundles")
  library("HierarchicalEdgeBundles")
  E(tmp.g)$width <- 1
#   plot(1)
  plotHEB(tmp.g, phylo, beta = 0.85,
          cex = 0.7, type = "fan",
          tip.color = dlbcl.modules,
          e.cols = E(tmp.g)$color)


  # PANEL COLORKEY
  par(xaxs = "r", yaxs = "r", mar = c(5.5, 0, 5.5, 0) + 0.5)
  plotColorkey(keybreaks, keycols )
  rect(ybottom = par("usr")[3], xleft = -dlbcl.par$threshold,
       xright = dlbcl.par$threshold, ytop = par("usr")[4], lwd = 0.5)
  axis(1, labels = FALSE)
  axis(3)

   # PANEL MODULE OVERVIEW
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

  dev.off()
}
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

stopifnot(all(names(mod.genes) == names(num2names)))

# Order by rowSums
dlbcl.exp.sub <- lapply(mod.genes, function(ensg) dlbcl.exp[ensg, ensg])
dlbcl.exp.sub <- lapply(dlbcl.exp.sub, function(x) {
  o <- order(rowSums(x), decreasing = TRUE)
  return(x[o,o])
})
tmp <- lapply(dlbcl.exp.sub, rownames)
m <- max(table(dlbcl.modules))  # Get largest module

# Construct table
dlbcl.mod.tab.genes <-
  sapply(tmp, function(x) unname(c(map2hugo(x), rep(NA, m - length(x)))))

# First letter capitalized
cgroup <- capitalize(cleanName(names(dlbcl.module.genes)))
# cgroup <- paste0(cgroup, num2names)

colnames(dlbcl.mod.tab.genes) <- # Number of gens in each module
  paste0("n = ", colSums(!is.na(dlbcl.mod.tab.genes)))

is.ensg <- structure(grepl("^ENSG", dlbcl.mod.tab.genes),
                     dim = dim(dlbcl.mod.tab.genes))
nn <- 20
seq_tmp <- seq_len(min(nn, nrow(dlbcl.mod.tab.genes)))
latex(dlbcl.mod.tab.genes[seq_tmp, ],
      cgroup = cgroup,
      size = "tiny",
      caption = paste("The identified modules, their sizes, and",
                      "member genes. The genes are sorted decreasingly",
                      "by their intra-module connectivity (sum of the incident",
                      "edge weights). Only the top", nn, "genes are shown."),
      cellTexCmds = ifelse(is.ensg, "tiny", "")[seq_tmp, ],
      caption.loc = "bottom",
      label = "tab:dlbcl_mod_tab",
      landscape = FALSE,
      file = "")

colnames(dlbcl.mod.tab.genes) <-
  paste0(cgroup, " (", colnames(dlbcl.mod.tab.genes), ")")
write.table(dlbcl.mod.tab.genes, file = "SuppB_module_genes.txt",
            quote = FALSE, sep = "\t")
## ---- end ----


## ---- GO_tabs ----
tmp <- dlbcl.go.analysis
names(tmp) <- cleanName(names(tmp))
go.table <- do.call(rbind, tmp)
write.table(go.table, file = "SuppB_go_table.txt", quote = FALSE, sep = "\t")
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


load("metadata.RData")
library("WGCNA")
library("survival")

figure3 <- "figure3.jpg"
jpeg(figure3, height = 1.2*7, width = 2*7, units = "in", res = 200)
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
    eg <- res$eigengenes
    col <- gsub("^ME", "", colnames(eg))
    eg <- as.data.frame(lapply(eg, function(x) x/sd(x))) # Standardize

    #   # Multivariate
    cph.fits <- coxph(meta$OS ~ ., data = eg, x = TRUE, y = TRUE)
    dats <- get.cis(cph.fits)
    dats <- dats[, -c(2,6)]
    plot.cis(dats)

    # Univariate
    # cph.fits <- lapply(eg, function(x) coxph(meta$OS ~ x, x = TRUE, y = TRUE))
    # dats <- do.call(rbind, lapply(cph.fits, get.cis))
    # dats <- dats[, -c(2,6)]
    # plot.cis(dats)

    for (i in seq_len(ncol(eg))) {
      if (col[i] == the.module) {
        eg.i <- eg[, paste0("ME", col[i])]
        eg.high <- eg.i >= mean(eg.i)
        plot(survfit(meta$OS ~ factor(eg.high)), conf.int = TRUE,
             main = "", col = c(col[i], "black"), lwd = 2,
             axes = FALSE,
             xlab = "years", ylab = "Survival proportion")
        axis(1); axis(2)
        legend("bottom", bty = "n", lwd = 2,
               legend = c("High eigengene", "Low eigengene"),
               col = c(col[i], "black"), horiz = TRUE)
        legend("topright", legend = paste(cleanName(col[i]), "eigengene"),
               bty = "n", text.col = col[i])
      }
    }
    title(switch(j, "GSE10846 CHOP", "GSE10846 R-CHOP"), xpd = TRUE)
  }
}
dev.off()
## ---- end ----



## ---- dlbcl_the_module ----

# Refit the model on the THE module only
the.genes <- names(which(dlbcl.modules == the.module))

if (!exists("the.rcm") || recompute) {

  gep.sub <- lapply(gep, function(x) exprs(x)[the.genes, ])
  dlbcl.ns  <- sapply(gep.sub, ncol)
  dlbcl.S   <- lapply(gep.sub, function(x) correlateR::scatter(t(x)))
  nu <- sum(dlbcl.ns) + ncol(dlbcl.S[[1]]) + 1
  psi <- c(nu - ncol(dlbcl.S[[1]]) - 1)*correlateR:::pool(dlbcl.S, dlbcl.ns)

  the.rcm <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                        Psi.init = psi, nu.init = nu, eps = 0.01,
                        max.ite = 1500)
  dimnames(the.rcm$Psi) <- dimnames(dlbcl.S[[1]])

  resave(the.rcm, file = "saved.RData")
}

# Fit model with random genes of the same size a the THE module
set.seed(10)
if (!exists("rand.rcm") || recompute) {
  rand.rcm <- vector("list", 500)

  for (i in seq_along(rand.rcm)) {
    rand.genes <- sample(names(var.pool), length(the.genes))
    gep.sub <- lapply(gep, function(x) exprs(x)[rand.genes, ])
    dlbcl.ns  <- sapply(gep.sub, ncol)
    dlbcl.S   <- lapply(gep.sub, function(x) correlateR::scatter(t(x)))
    nu <- sum(dlbcl.ns) + ncol(dlbcl.S[[1]]) + 1
    psi <- nu*correlateR:::pool(dlbcl.S, dlbcl.ns)

    rand.rcm.i <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                          Psi.init = psi, nu.init = nu, eps = 0.01,
                          max.ite = 1500)
    dimnames(rand.rcm.i$Psi) <- dimnames(dlbcl.S[[1]])
    rand.rcm[[i]] <- rand.rcm.i
    cat(i, "\n"); flush.console()
  }

  resave(rand.rcm, file = "saved.RData")
}


# Do the test for homogeneitiy
set.seed(100)
if (!exists("homogeneity.rcm") || recompute) {
  homogeneity.rcm <- vector("list", 500)
  dat <- do.call(rbind, lapply(gep, function(x) t(exprs(x)[the.genes, ])))
  dat <- data.frame(dat)
  class.lab <- rep(names(gep), sapply(gep, ncol))

  for (i in seq_along(homogeneity.rcm)) {

    s.class.lab <- sample(class.lab)
    gep.sub <- split(dat, s.class.lab)
    dlbcl.ns  <- sapply(gep.sub, nrow)
    dlbcl.S   <- lapply(gep.sub, scatter)
    nu <- sum(dlbcl.ns) + ncol(dlbcl.S[[1]]) + 1
    psi <- nu*correlateR:::pool(dlbcl.S, dlbcl.ns)
    rand.rcm.i <- fit.rcm(S = dlbcl.S, ns = dlbcl.ns, verbose = TRUE,
                          Psi.init = psi, nu.init = nu, eps = 0.01,
                          max.ite = 2500)
    dimnames(rand.rcm.i$Psi) <- dimnames(dlbcl.S[[1]])
    homogeneity.rcm[[i]] <- rand.rcm.i
    cat(i, "\n"); flush.console()
  }

  resave(homogeneity.rcm, file = "saved.RData")
}


get.ICC <- function(x) {
  with(x, ICC(nu, nrow(Psi)))
}

get.TestPValue <- function(the.list, the.object) {
  n <- sum(sapply(the.list, "[[", "nu") < the.object$nu) + 1
  N <- length(the.list) + 1
  return(n / N)
}


## ---- end ----
