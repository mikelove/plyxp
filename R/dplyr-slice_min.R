#' @name slice_min
#' @title Slice minimum rows or columns of a PlySummarizedExperiment
#' @description
#' `slice_min()` selects rows or columns with the smallest values of a variable.
#' @param .data a PlySummarizedExperiment object
#' @param ... expressions wrapped in `rows()` or `cols()` giving the ordering
#'   variable. Only the `rows()` and `cols()` contexts are available here.
#' @param n Number of rows/columns to select per group. Cannot be combined
#'   with `prop`.
#' @param prop Fraction of rows/columns to select per group. Cannot be
#'   combined with `n`.
#' @param with_ties Should ties be kept together? Default `TRUE`.
#' @param na_rm Should missing values be removed? Default `FALSE`.
#' @param .preserve When `FALSE` (default), grouping structure is recomputed
#'   based on the resulting data.
#' @examples
#' slice_min(se_simple, rows(length), n = 2)
#'
#' gse <- group_by(se_simple, rows(direction))
#' slice_min(gse, rows(length), n = 1)
#' @export
slice_min.PlySummarizedExperiment <- function(
    .data, ..., n = 1, prop = NULL, with_ties = TRUE, na_rm = FALSE,
    .preserve = FALSE) {
  plyxp(.data, slice_min_se_impl, ...,
    n = n, prop = prop, with_ties = with_ties, na_rm = na_rm,
    .preserve = .preserve
  )
}

slice_min_se_impl <- function(
    .data, ..., n = 1, prop = NULL, with_ties = TRUE, na_rm = FALSE,
    .preserve = FALSE) {
  .env <- caller_env()

  if (!is.null(n) && !is.null(prop)) {
    rlang::abort("Cannot use both `n` and `prop` in `slice_min()`.")
  }

  # Build the transform as a quoted lambda so it can be inlined by plyxp_quos.
  # n/prop/with_ties/na_rm are substituted at construction time via bquote(),
  # leaving `plyxp:::ctx:::n` unsubstituted so it resolves to the current
  # group size at evaluation time (matching the same pattern as slice()'s
  # in_bounds transform).
  take_min <- if (!is.null(prop)) {
    bquote(\(.x) {
      na_last <- if (.(na_rm)) TRUE else NA
      k <- ceiling(.(prop) * `plyxp:::ctx:::n`)
      if (.(with_ties)) {
        which(rank(.x, ties.method = "min", na.last = na_last) <= k)
      } else {
        order(.x, na.last = na_last)[seq_len(min(k, `plyxp:::ctx:::n`))]
      }
    })
  } else {
    bquote(\(.x) {
      na_last <- if (.(na_rm)) TRUE else NA
      k <- min(.(n), `plyxp:::ctx:::n`)
      if (.(with_ties)) {
        which(rank(.x, ties.method = "min", na.last = na_last) <= .(n))
      } else {
        order(.x, na.last = na_last)[seq_len(k)]
      }
    })
  }

  quos <- plyxp_quos(
    ...,
    .ctx = c("assays", "rows", "cols"),
    .trans = list(rows = take_min, cols = take_min)
  )

  ctxs <- vapply(quos, attr, FUN.VALUE = "", which = "plyxp:::ctx")
  if (any(err <- ctxs %in% "assays")) {
    plyxp_assays_cannot(do = "slice_min", review = err)
  }
  if (any(quos_is_named(quos))) {
    plyxp_should_not_named("slice_min")
  }

  names(quos) <- sprintf("..slice_min_%i", seq_along(quos))
  nms <- names(quos)
  quos <- quos_enforce_named(quos, nms)

  mask <- new_plyxp_manager.SummarizedExperiment(obj = .data)
  poke_ctx_local("plyxp:::caller_env", .env)
  poke_ctx_local("plyxp:::manager", mask)
  poke_ctx_local("plyxp:::dplyr_verb", "slice_min")
  mask <- plyxp_evaluate(mask, quos, ctxs, nms, .env)

  chops <- mask$apply(
    .f = function(mask, name) {
      added <- mask$added
      names(added) <- added
      lapply(added, mask$get_chop)
    },
    .on_masks = c("rows", "cols")
  )
  chops <- flatten_slice_chops(chops, groups = group_data_se_impl(.data))
  which_ctx <- unique(ctxs)
  switch(length(which_ctx),
    `1` = switch(which_ctx,
      rows = plyxp_slice_se(.data, chops$rows, .preserve = .preserve),
      cols = plyxp_slice_se(.data, , chops$cols, .preserve = .preserve)
    ),
    `2` = plyxp_slice_se(.data, chops$rows, chops$cols, .preserve = .preserve)
  )
}
