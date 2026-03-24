#' Estimate K-cell-type fractions with BayesDeBulk
#'
#' This function runs \code{BayesDeBulk} on bulk expression data and returns
#' estimated cell-type fractions together with the fitted cell-type-specific
#' expression profiles.
#'
#' @param exprB Matrix or data frame of bulk expression values, formatted as
#'   expected by \code{BayesDeBulk}.
#' @param markers Named list of marker genes, where each element contains the
#'   marker set for one cell type.
#' @param seed Random seed.
#' @param n.iter Total number of MCMC iterations passed to \code{BayesDeBulk}.
#' @param burn.in Number of burn-in iterations passed to \code{BayesDeBulk}.
#' @param ... Additional arguments passed to \code{BayesDeBulk}.
#'
#' @return A list with components:
#' \describe{
#'   \item{cell_fraction}{Data frame of estimated cell-type fractions with one
#'     row per sample and one column per cell type.}
#'   \item{cell_expression}{Cell-type-specific expression estimates returned by
#'     \code{BayesDeBulk}.}
#' }
#'
#' @importFrom BayesDeBulk BayesDeBulk
#' @export
pi_estimation_K <- function(exprB,
                            markers,
                            seed   = 1,
                            n.iter = 10000,
                            burn.in = 1000,
                            ...) {
  cell.type<-names(markers)
  index.matrix<-NULL
  for (s in 1:length(cell.type)){
    for (k in 1:length(cell.type)){
      if (s!=k){
        mg<-match(markers[[s]],markers[[k]])
        index.matrix<-rbind(index.matrix,cbind(rep(cell.type[s],sum(is.na(mg))),
                                               rep(cell.type[k],sum(is.na(mg))),markers[[s]][is.na(mg)]))
      }
    }
  }
  set.seed(seed)

  fit <- BayesDeBulk(
    n.iter  = n.iter,
    burn.in = burn.in,
    Y       = list(exprB),
    markers = index.matrix,
    ...
  )

  cell_fraction <- as.data.frame(fit$cell.fraction)

  list(
    cell_fraction  = cell_fraction,
    cell_expression = fit$cell.expression
  )
}
