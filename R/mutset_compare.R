#' compare two `BiologyAAmutSet` objects by Chi-Squared test
#'
#' @param x BiologyAAmutSet
#' @param y BiologyAAmutSet
#' @param all return all mutations. if `FALSE`, only return mutations with
#' p < 0.05
#' @param use_numbering use numbering or not
#' @param min_count minimal count of a mutation to perform test
#' @param bysite merge the mutations at the same site
#'
#' @return tibble
#' @export
#'
mutset_compare <- function(x, y, min_count = 2,
                           all = FALSE, bysite = FALSE, use_numbering = FALSE) {
  muts_x <- count_muts(x, use_numbering = use_numbering, bysite = bysite) %>%
    tibble::column_to_rownames("mut_aa") %>%
    rowSums()
  muts_x_flt <- muts_x[muts_x > min_count]
  muts_y <- count_muts(y, use_numbering = use_numbering, bysite = bysite) %>%
    tibble::column_to_rownames("mut_aa") %>%
    rowSums()
  muts_y_flt <- muts_y[muts_y > min_count]

  Vmut <- union(names(muts_x_flt), names(muts_y_flt))
  ord <- stringr::str_extract(Vmut, "\\d+") %>%
    as.numeric() %>%
    order()
  Vmut <- Vmut[ord]

  mut_test <- function(m) {
    x_count <- unname(ifelse(is.na(muts_x[m]), 0, muts_x[m]))
    x_none <- length(x) - x_count
    x_ratio <- round(x_count / length(x), 3)

    y_count <- unname(ifelse(is.na(muts_y[m]), 0, muts_y[m]))
    y_none <- length(y) - y_count
    y_ratio <- round(y_count / length(y), 3)

    mtx <- matrix(c(x_count, x_none, y_count, y_none), nrow = 2)

    pvalue <- suppressWarnings(chisq.test(mtx)$p.value)

    res <- c( # nolint
      "p" = pvalue,
      "x_count" = x_count,
      "x_ratio" = x_ratio,
      "y_count" = y_count,
      "y_ratio" = y_ratio
    )
  }

  all_test <- purrr::map_dfr(Vmut, mut_test) %>%
    dplyr::mutate(mut = Vmut, .before = 1)

  if (all == FALSE) {
    res <- all_test %>% dplyr::filter(.data[["p"]] < 0.05)
  } else {
    res <- all_test
  }
  return(res)
}
