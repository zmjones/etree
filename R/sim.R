args <- commandArgs(TRUE)
cores <- as.integer(args[1])
ntree <- as.integer(args[2])
boost_iter <- as.integer(args[3])

library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(party)
library(nlme)
library(mboost)
library(reshape2)
library(ggplot2)

source("R/tree.R")
source("R/resample.R")
source("R/boost.R")

prefix <- "data/"
files <- list.files(prefix, "^sim")

mse <- function(y, y_hat) mean((y - y_hat)^2)
comb <- function(...) apply(do.call(cbind, list(...)), 1, mean)
comb_all <- function(...) as.data.frame(do.call(rbind, list(...)))
plot_comparison <- function(out, title, fname) {
    out$rho <- paste0("rho == ", out$rho)
    out$tau <- paste0("tau == ", out$tau)
    out <- melt(out, id.vars = c("tau", "rho", "groups"))
    p <- ggplot(out, aes(variable, value, group = variable))
    p <- p + facet_grid(rho ~ tau, scales = "free", labeller = label_parsed)
    p <- p + geom_boxplot()
    p <- p + labs(x = NULL, y = "MSE", title = title)
    p <- p + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0("figures/", fname, ".png"), p, width = 8, height = 8)
}

cl <- makeCluster(cores)
registerDoParallel(cl)

pkgs <- c("party", "nlme", "foreach", "reshape2", "ggplot2", "iterators", "stringr")

bag_mse_all <- foreach(f = files, .packages = pkgs, .combine = comb_all) %dopar% {
    load(paste0(prefix, f))
    bag_mse <- foreach(d = df, .combine = "rbind", .packages = pkgs) %dopar% {
        args <- unlist(str_split(f, "_"))
        groups <- unique(d[, "group"])
        train_idx <- unlist(lapply(groups, function(x) which(d$group == x)[1:(nrow(d) / length(groups) / 2)]))
        train <- d[train_idx, ]
        test <- d[-train_idx, ]
        idx_boot <- foreach(icount(ntree)) %do% bootstrap(train)
        idx_block <- foreach(icount(ntree)) %do% moving_block_bootstrap(train, "group", 4)
        form <- formula(paste0("y ~ ", paste0("X", 1:(ncol(train) - 2), collapse = "+")))
        re_form <- ~ 1 | group
        fit_re_boot <- foreach(idx = idx_boot) %do% fit_re_ctree(train[idx, ], form, re_form, corAR1())
        fit_re_block <- foreach(idx = idx_block) %do% fit_re_ctree(train[idx, ], form, re_form, corAR1())
        pred_boot <- foreach(idx = idx_boot, .combine = comb) %do% predict(ctree(form, train[idx, ]), test)
        pred_block <- foreach(idx = idx_block, .combine = comb) %do% predict(ctree(form, train[idx, ]), test)
        pred_re_boot <- foreach(fit = fit_re_boot, .combine = comb) %do% predict_re_ctree(fit, test, "group")
        pred_re_block <- foreach(fit = fit_re_block, .combine = comb) %do% predict_re_ctree(fit, test, "group")
        c("bootstrap" = mse(test$y, pred_boot),
          "moving block bootstrap" = mse(test$y, pred_block),
          "random effects and boostrap" = mse(test$y, pred_re_boot),
          "random effects and moving block bootstrap" = mse(test$y, pred_re_block),
          "tau" = as.integer(str_extract(args[6], "^[1-2]")),
          "rho" = as.numeric(args[5]) / 100,
          "groups" = as.integer(args[3]))
    }
}
write.csv(bag_mse_all, "data/bag_mse.csv", row.names = FALSE)

plot_comparison(bag_mse_all, "Bagged (50) Trees on LME process with 10 groups, 100 obs. per group", "bag")

boost_mse_all <- foreach(f = files, .combine = comb_all, .packages = pkgs) %dopar% {
    load(paste0(prefix, f))
    boost_mse <- foreach(d = df, .combine = "rbind", .packages = pkgs) %dopar% {
        args <- unlist(str_split(f, "_"))
        groups <- unique(d[, "group"])
        train_idx <- unlist(lapply(groups, function(x) which(d$group == x)[1:(nrow(d) / length(groups) / 2)]))
        train <- d[train_idx, ]
        test <- d[-train_idx, ]
        idx_boot <- foreach(icount(boost_iter)) %do% bootstrap(train)
        idx_move_block <- foreach(icount(boost_iter)) %do% moving_block_bootstrap(train, "group", 4)
        form <- formula(paste0("y ~ ", paste0("X", 1:(ncol(train) - 2), collapse = "+")))
        re_form <- ~ 1 | group
        fit_boot <- fit_quad_boost_ctree(form, train, idx_boot)
        fit_block <- fit_quad_boost_ctree(form, train, idx_move_block)
        fit_re_boot <- fit_quad_boost_re_ctree(form, re_form, train, idx_boot)
        fit_re_block <- fit_quad_boost_re_ctree(form, re_form, train, idx_move_block)
        pred_boot <- predict_quad_boost_ctree(fit_boot, test)
        pred_block <- predict_quad_boost_ctree(fit_block, test)
        pred_re_boot <- predict_quad_boost_re_ctree(fit_re_boot, test, "group")
        pred_re_block <- predict_quad_boost_re_ctree(fit_re_block, test, "group")
        c("boosted ctree bootstrap" = mse(test$y, pred_boot),
          "boosted ctree boost block bootstrap" = mse(test$y, pred_block),
          "boosted re ctree bootstrap" = mse(test$y, pred_re_boot),
          "boosted re ctree moving block bootstrap" = mse(test$y, pred_re_block),
          "tau" = as.integer(str_extract(args[6], "^[1-2]")),
          "rho" = as.numeric(args[5]) / 100,
          "groups" = as.integer(args[3]))
    }
}
write.csv(boost_mse_all, "data/boost_mse.csv", row.names = FALSE)

plot_comparison(boost_mse_all, "Boosted (100) Trees on LME process with 10 groups, 100 obs. per group", "boost")
