args <- commandArgs(TRUE)
n <- as.integer(args[1])

set.seed(1988)

library(party)
library(ggplot2)
library(reshape2)
library(foreach)
library(iterators)

source("R/boost.R")

x <- runif(n, -4, 4)
y <- sin(x)
df <- data.frame(x, y)

fit_tree <- ctree(y ~ x, df)
df$tree <- as.numeric(predict(fit_tree))

p <- ggplot(df, aes(x, y))
p <- p + geom_point(alpha = .5)
p <- p + geom_line(aes(x, tree), colour = "blue")
p <- p + theme_bw()
ggsave("figures/cart_approximation.png", p, width = 8, height = 4)

object <- party:::ctreedpp(y ~ x, df)
fit_bag <- cforest(y ~ x, df, control = cforest_unbiased(mtry = 1))
out <- sapply(fit_bag@ensemble, function(x)
    unlist(party:::R_getpredictions(x, party:::R_get_nodeID(x, object@inputs, 0))))
out <- as.data.frame(out[, sample(1:ncol(out), 20)])
out$obs <- 1:nrow(out)
out <- melt(out, id.vars = "obs")
out$x <- x
out$y <- y
out$bag <- as.numeric(predict(fit_bag))
p <- ggplot(out, aes(x, y, group = variable))
p <- p + geom_point(alpha = .5)
p <- p + geom_line(aes(x, bag), colour = "red")
p <- p + geom_line(aes(x, value), alpha = .25, colour = "blue")
p <- p + theme_bw()
ggsave("figures/bag_approximation.png", p, width = 8, height = 4)

boost_iter <- 100
idx <- lapply(1:boost_iter, function(x) sample(1:nrow(df), nrow(df), TRUE))
fit <- fit_quad_boost_ctree(y ~ x, df, idx, boost_iter, .5)

pred <- matrix(NA, nrow(df), boost_iter)
pred[, 1] <- mean(df$y)
for (i in 2:boost_iter) {
    sub_fit <- list("base" = fit$base[1:i], "iter" = fit$iter, "rho" = fit$rho, "form" = fit$form)
    pred[, i] <- predict_quad_boost_ctree(sub_fit, df)
}

df <- cbind(df, pred)
df$tree <- NULL
colnames(df)[ncol(df)] <- "final"
df <- melt(df, id.vars = c("x", "y", "final"))
out <- df[df$variable %in% as.character(c(1, 10, 15, 25, 30, 50, 75)), ]
p <- ggplot(out, aes(x, y, group = variable, colour = variable))
p <- p + geom_point(colour = "black", alpha = .5)
p <- p + geom_line(aes(x, value))
p <- p + scale_colour_brewer("Iteration", palette = "Spectral")
p <- p + theme_bw()
ggsave("figures/boost_approximation.png", width = 8, height = 4)

mse <- function(y, y_hat) mean((y - y_hat)^2)
out <- apply(pred, 2, function(x) mse(y, x))
out <- data.frame("mse" = out, "iter" = 1:length(out))

p <- ggplot(out, aes(iter, mse))
p <- p + geom_line()
p <- p + labs(x = "boosting iterations", y = "MSE")
p <- p + theme_bw()
ggsave("figures/boost_mse.png", width = 8, height = 4)
