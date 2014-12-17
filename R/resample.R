subsample <- function(df) sample(1:nrow(df), nrow(df) * .632, FALSE)

bootstrap <- function(df) sample(1:nrow(df), nrow(df), TRUE)

block_bootstrap <- function(df, group) {
    groups <- unique(df[, group])
    selected <- sample(groups, length(groups), TRUE)
    which(df[, group] %in% selected)
}

moving_block_bootstrap <- function(df, group, l) {
    groups <- unique(df[, group])
    idx <- vector("integer", nrow(df))
    for (g in groups) {
        tmp <- which(df[, group] == g)
        idx_g <- vector("integer", length(tmp))
        while (any(idx_g == 0)) {
            pt <- which(idx_g == 0)
            tmp_l <- min(length(tmp), l)
            if (length(pt) %% tmp_l != 0)
                tmp_l <- length(pt)
            cand <- sample(tmp, 1)
            if (cand <= max(tmp) - tmp_l & cand >= min(tmp) + tmp_l)
                idx_g[pt[1]:(pt[1] + tmp_l - 1)] <- cand:(cand + tmp_l - 1)
        }
        st_pt <- which(idx == 0)[1]
        idx[st_pt:(st_pt + length(idx_g) - 1)] <- idx_g
    }
    idx
}

k_cv <- function(df, folds) sample(rep(1:folds, length.out = nrow(df)))
