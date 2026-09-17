# User interface for the design-based grouped horseshoe on the loadings of a meshed LMC.

design_groups <- function(blocks, grouping = c("design", "element", "none")) {
  # blocks: integer vector of length q; outcomes that share an observation pattern share a label.
  # Returns a q x q integer matrix: 0 = not shrunk, g >= 1 = index of the local scale.
  grouping <- match.arg(grouping)
  q <- length(blocks)
  G <- matrix(0L, q, q)
  if (grouping == "none" || q < 2) return(G)
  if (is.unsorted(blocks)) stop("order the outcomes so that blocks are contiguous and nondecreasing in 'blocks'")
  id <- 0L
  for (j in 2:q) {
    if (grouping == "element") {
      for (l in seq_len(j - 1)) { id <- id + 1L; G[j, l] <- id }
    } else {
      for (b in unique(blocks[seq_len(j - 1)])) { id <- id + 1L; G[j, which(blocks[seq_len(j - 1)] == b)] <- id }
    }
  }
  G
}

spmeshed_ghs <- function(y, x, coords, blocks = NULL, grouping = c("design", "element", "none"), prior = list(), ...) {
  # Meshed LMC with k = q factors and a (grouped) horseshoe on the strictly lower loadings.
  # blocks = NULL infers the blocks from the missingness pattern of the columns of y.
  grouping <- match.arg(grouping)
  y <- as.matrix(y); q <- ncol(y)
  if (is.null(blocks)) {
    pat <- apply(!is.na(y), 2, function(z) paste(which(z), collapse = ","))
    blocks <- match(pat, unique(pat))
  }
  if (length(blocks) != q) stop("'blocks' must have one entry per outcome")
  G <- design_groups(blocks, grouping)
  prior$lambda_group <- G
  fit <- spmeshed(y = y, x = x, coords = coords, k = q, prior = prior, ...)
  fit$ghs <- list(group = G, blocks = blocks, grouping = grouping, n_obs = colSums(!is.na(y)),
                  outcome_names = if (is.null(colnames(y))) paste0("y", seq_len(q)) else colnames(y))
  class(fit) <- c("spmeshed_ghs", class(fit))
  fit
}

zero_distance_cor <- function(fit, probs = c(0.025, 0.5, 0.975)) {
  # Posterior of the zero-distance cross-correlations rho_ij = C_ij(0) / sqrt(C_ii(0) C_jj(0)), C(0) = Lambda Lambda'.
  L <- fit$lambda_mcmc
  R <- vapply(seq_len(dim(L)[3]), function(s) stats::cov2cor(tcrossprod(L[, , s])), matrix(0, dim(L)[1], dim(L)[1]))
  list(mean = apply(R, 1:2, mean), quantiles = apply(R, 1:2, stats::quantile, probs = probs), draws = R)
}

shrinkage_weights <- function(fit, probs = c(0.025, 0.5, 0.975)) {
  # Diagnostic only. Do not threshold these to decide on borrowing: with large n_j / tau_j^2 every nonzero loading
  # gives a weight near zero. Use borrowing_share() and borrowing_map().
  # kappa_g = 1 / (1 + n_j zeta^2 gamma_g^2 / tau_j^2): close to 1 when outcome j does not borrow from the block.
  if (is.null(fit$ghs) || max(fit$ghs$group) == 0) stop("fit has no horseshoe groups")
  G <- fit$ghs$group; nm <- fit$ghs$outcome_names
  out <- lapply(seq_len(max(G)), function(g) {
    j <- which(G == g, arr.ind = TRUE)[1, 1]; cols <- which(G[j, ] == g)
    kap <- 1 / (1 + fit$ghs$n_obs[j] * (fit$ghs_tau_mcmc * fit$ghs_lam_mcmc[g, ])^2 / fit$tausq_mcmc[j, ])
    data.frame(group = g, outcome = nm[j], borrows_from = paste(nm[cols], collapse = ", "), size = length(cols),
               pr_borrow = mean(kap < 0.5), mean = mean(kap), t(stats::quantile(kap, probs)), check.names = FALSE)
  })
  do.call(rbind, out)
}

print.spmeshed_ghs <- function(x, ...) {
  cat("Meshed LMC with", x$ghs$grouping, "horseshoe loadings:", ncol(x$ghs$group), "outcomes,",
      max(x$ghs$group), "groups, blocks", paste(x$ghs$blocks, collapse = " "), "\n")
  if (max(x$ghs$group) > 0) print(borrowing_share(x)$table, digits = 3, row.names = FALSE)
  invisible(x)
}
