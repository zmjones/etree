set.seed(1987)

library(mvtnorm)
library(foreach)
library(iterators)

args <- commandArgs(TRUE)
n <- as.integer(args[1])
groups <- as.integer(args[2])
k <- as.integer(args[3])
phi <- as.numeric(args[4])
tau <- as.numeric(args[5])
m <- as.integer(args[6])
fname <- paste0(paste("../data/sim", n, groups, k, phi * 100, tau, sep = "_"), ".rda")

gen_data <- function(alpha, beta, X, groups, n, k, error_model, ...) {
    df <- vector("list", groups)
    mu <- lapply(1:groups, function(i) alpha[i] + X[[i]] %*% beta)
    df <- as.data.frame(do.call("rbind", lapply(1:groups, function(i)
        cbind(mu[[i]] + as.numeric(arima.sim(error_model, n, ...)), i, X[[i]]))))
    colnames(df) <- c("y", "group", paste0("X", 1:k))
    df
}

alpha <- rnorm(groups, 0, tau)
beta <- runif(k, -5, 5)
X <- lapply(1:groups, function(i) rmvnorm(n / groups, rep(0, k), diag(k), method = "chol"))
df <- foreach(icount(m)) %do% gen_data(alpha, beta, X, groups, n / groups, k,
                                         list(order = c(1, 0, 0), ar = phi), sd = 1)
save(df, file = fname)
