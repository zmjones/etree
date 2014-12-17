fit_re_ctree <- function(df, fe_form, re_form, corr = NULL, max_iter = 100, tol = 1e-5, ...) {
    y <- as.character(fe_form[[2]])
    df$y_star <- df[, y]
    rhs <- paste(attr(terms(fe_form), "term.labels"), collapse = "+")
    ll_old <- -Inf
    flag <- ""
    where <- rep.int(0, nrow(df))
    weights <- rep(1, nrow(df))
    storage.mode(where) <- "integer"
    storage.mode(weights) <- "double"
    object <- party:::ctreedpp(formula(paste0("y_star ~ ", rhs)), df)
    fitmem <- party::ctree_memory(object, TRUE)
    for (i in 1:max_iter) {
        tree <- party:::R_TreeGrow(object, weights, fitmem, party::ctree_control(...), where)
        party:::R_remove_weights(tree, TRUE)
        ind <- party:::R_get_nodeID(tree, object@inputs, 0)
        if (length(unique(ind)) == 1) {
            flag <- "no partitions"
            tree <- NULL
            re <- NULL
            break
        } else {
            df$ind <- factor(ind)
            re <- nlme::lme(formula(paste0(y, "~ -1 + ind")),
                            df, re_form, corr, method = "REML", keep.data = FALSE,
                            control = list(maxIter = 500, msMaxEval = 500))
            ll_new <- re$logLik
            if (i <= max_iter - 1 & i != 1 & i != max_iter) {
                if (abs(ll_new - ll_old) < tol) break
                object@responses@y_star <- df[, y] +
                    (fitted(re, level = 0) - fitted(re, level = dim(re$fitted)[2] - 1))
            }
            ll_old <- ll_new
        }
    }
    if (i == max_iter) {
        message("convergence failed: maximum iterations reached")
        flag <- "maximum iterations"
    }
    list("tree" = tree, "re" = re, "iterations" = i, "fe_form" = fe_form, "flag" = flag)
}

predict_re_ctree <- function(fit, df, group, ...) {
    ## extract random effects and use if available
    id <- as.matrix(df[, group])
    re <- ranef(fit$re)
    out <- matrix(NA, nrow(df), length(re))
    for (d in 1:length(re)) {
        for (i in 1:length(unique(id[, d]))) {
            idx <- id[, d] == unique(id)[i, d]
            if (length(re) > 1) {
                re_est <- re[[d]][unique(id)[i, d], ]
                out[idx, d] <- re[[d]][unique(id)[i, d], ]
            } else {
                re_est <- re[[d]][unique(id)[i, d]]
                out[idx, d] <- re[[d]][unique(id)[i, d]]
            }
        }
    }
    ## extract fixed effects from tree
    object <- party:::ctreedpp(fit$fe_form, data = df)
    ind <- party:::R_get_nodeID(fit$tree, object@inputs, 0)
    fe <- unlist(party:::R_getpredictions(fit$tree, ind))
    fe + apply(out, 1, sum, na.rm = TRUE)
}
