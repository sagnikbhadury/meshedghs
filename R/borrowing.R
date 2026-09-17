# Posterior inference on borrowing: an effect size (borrowing share) and a decoupled selection step.
# These replace thresholding of horseshoe shrinkage weights, which saturates when n_j / tau_j^2 is large.

.group_matrix_for <- function(fit, blocks = NULL, grouping = "design") {
  if (!is.null(fit$ghs) && is.null(blocks)) return(list(G = fit$ghs$group, names = fit$ghs$outcome_names))
  q <- dim(fit$lambda_mcmc)[1]
  if (is.null(blocks)) stop("supply 'blocks' for a fit that was not produced by spmeshed_ghs()")
  list(G = design_groups(blocks, grouping), names = paste0("y", seq_len(q)))
}

borrowing_share <- function(fit, blocks = NULL, grouping = "design", eps = 0.05, probs = c(0.025, 0.5, 0.975)) {
  # omega_{j,r} = sum_{i in B_r, i < j} lambda_ji^2 / sum_h lambda_jh^2, per posterior draw.
  gm <- .group_matrix_for(fit, blocks, grouping); G <- gm$G; nm <- gm$names; L <- fit$lambda_mcmc
  draws <- sapply(seq_len(max(G)), function(g) {
    j <- which(G == g, arr.ind = TRUE)[1, 1]; cols <- which(G[j, ] == g)
    apply(L, 3, function(M) sum(M[j, cols]^2) / sum(M[j, ]^2))
  })
  tab <- do.call(rbind, lapply(seq_len(max(G)), function(g) {
    j <- which(G == g, arr.ind = TRUE)[1, 1]; cols <- which(G[j, ] == g)
    data.frame(group = g, outcome = nm[j], borrows_from = paste(nm[cols], collapse = ", "), size = length(cols),
               t(stats::quantile(draws[, g], probs)), pr_above_eps = mean(draws[, g] > eps), check.names = FALSE)
  }))
  list(table = tab, draws = draws, eps = eps)
}

borrowing_map <- function(fit, blocks = NULL, grouping = "design", eps = 0.05, level = 0.9) {
  # Decoupled selection (Hahn and Carvalho 2015): group soft-thresholding of the posterior mean loadings gives the
  # path; the posterior of the variance share given up, delta_j(S), decides how far along the path to go.
  gm <- .group_matrix_for(fit, blocks, grouping); G <- gm$G; nm <- gm$names; L <- fit$lambda_mcmc
  sh <- borrowing_share(fit, blocks, grouping, eps)$draws; Lbar <- apply(L, 1:2, mean)
  out <- lapply(2:nrow(G), function(j) {
    gs <- sort(unique(G[j, G[j, ] > 0])); if (!length(gs)) return(NULL)
    score <- sapply(gs, function(g) sqrt(sum(Lbar[j, G[j, ] == g]^2) / sum(G[j, ] == g)))
    drop <- integer(0)
    for (g in gs[order(score)]) {                       # drop groups in the order of the group-lasso path
      if (mean(rowSums(sh[, c(drop, g), drop = FALSE]) <= eps) >= level) drop <- c(drop, g) else break
    }
    data.frame(outcome = nm[j], group = gs, borrows_from = sapply(gs, function(g) paste(nm[which(G[j, ] == g)], collapse = ", ")),
               retained = !(gs %in% drop), path_score = score)
  })
  do.call(rbind, out)
}
