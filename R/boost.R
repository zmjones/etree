fit_quad_boost_re_ctree <- function(fe_form, re_form, df, idx, re_iter = 100, boost_iter = 100,
                                    rho = .1, corr = NULL, tol = 1e-5, ...) {
    y <- as.character(fe_form[[2]])
    df$y_star <- df[, y]
    rhs <- paste(attr(terms(fe_form), "term.labels"), collapse = "+")
    boost_form <- formula(paste0("y_star ~ ", rhs))
    ll_old <- -Inf
    for (i in 1:re_iter) {
        fe_fit <- fit_quad_boost_ctree(boost_form, df, idx, boost_iter, rho, alt_form = fe_form)
        fe <- predict_quad_boost_ctree(fe_fit, df)
        df$y_tilde <- df[, y] - fe
        re_fit <- nlme::lme(y_tilde ~ -1, df, re_form, corr, method = "REML", keep.data = FALSE,
                            control = list(maxIter = 500, msMaxEval = 500))
        ll_new <- re_fit$logLik
        if (i <= re_iter - 1 & i != 1 & i != re_iter) {
            if (abs(ll_new - ll_old) < tol) break
            df$y_star <- df[, y] +
                (fitted(re_fit, level = 0) - fitted(re_fit, level = dim(re_fit)$fitted)[2] - 1)
        }
        ll_old <- ll_new
    }
    if (i == re_iter)
        message("convergence failed: maximum iterations reached")
    list("boost" = fe_fit, "re" = re_fit, "re_iterations" = i, "fe_form" = fe_form)
}

fit_quad_boost_ctree <- function(form, df, idx, iter = 100, rho = .1, alt_form = NULL, ...) {
    est <- base <- vector("list", iter)
    y <- as.character(form[[2]])
    rhs <- paste(attr(terms(form), "term.labels"), collapse = "+")
    df$u <- df[, y]
    base_form <- formula(paste0("u ~ ", rhs))
    est <- mean(df[, y])
    where <- rep.int(0, nrow(df))
    storage.mode(where) <- "integer"
    object <- party:::ctreedpp(formula(paste0("u ~ ", rhs)), df)
    fitmem <- party::ctree_memory(object, TRUE)
    for (i in 1:(iter - 1)) {
        object@responses@variables$u <- df[, y] - est
        df$u <- df[, y] - est
        object <- party:::ctreedpp(formula(paste0("u ~ ", rhs)), df)
        weights <- ifelse(1:nrow(df) %in% idx[[i]], 1, 0)
        storage.mode(weights) <- "double"
        base[[i]] <- party:::R_TreeGrow(object, weights, fitmem,
                                        party::ctree_control(mincriterion = 0,
                                                             savesplitstats = FALSE, ...), where)
        party:::R_remove_weights(base[[i]], TRUE)
        ind <- party:::R_get_nodeID(base[[i]], object@inputs, 0)
        fe <- unlist(party:::R_getpredictions(base[[i]], ind))
        est <- est + rho * fe
    }
    if (!is.null(alt_form))
        form <- alt_form
    list("base" = base, "iter" = ifelse(i < iter - 1, i - 1, i), "rho" = rho, "form" = form)
}

predict_quad_boost_ctree <- function(fit, df) {
    object <- party:::ctreedpp(fit$form, data = df)
    out <- sapply(1:(length(fit$base) - 1), function(x) {
        ind <- party:::R_get_nodeID(fit$base[[x]], object@inputs, 0)
        unlist(party:::R_getpredictions(fit$base[[x]], ind))
    })
    apply(out, 1, function(x) sum(x) * fit$rho) - mean(df[, as.character(fit$form[[2]])])
}

predict_quad_boost_re_ctree <- function(fit, df, group) {
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
    fe <- predict_quad_boost_ctree(fit$boost, df)
    fe + apply(out, 1, sum, na.rm = TRUE)
}
