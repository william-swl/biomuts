#' amino acid feature index, raw value
#'
#' @format ## `aa_feature_raw`
#' A data frame with 20 rows and 5 columns:
#' \describe{
#'   \item{hydrophilicity}{Hydrophobicity (Prabhakaran, 1990)}
#'   \item{polarity}{Polarity (Grantham, 1974)}
#'   \item{charge}{Net charge (Klein et al., 1984)}
#'   \item{volume}{Side chain volume (Krigbaum-Komoriya, 1979)}
#' }
#' @source adjusted from [aaindex](https://www.genome.jp/aaindex/)
"aa_feature_raw"

#' amino acid feature index, scaled value from -1 to 1
#'
#' @format ## `aa_feature`
#' A data frame with 20 rows and 5 columns:
#' \describe{
#'   \item{hydrophilicity}{Hydrophobicity (Prabhakaran, 1990)}
#'   \item{polarity}{Polarity (Grantham, 1974)}
#'   \item{charge}{Net charge (Klein et al., 1984)}
#'   \item{volume}{Side chain volume (Krigbaum-Komoriya, 1979)}
#' }
#' @source adjusted from [aaindex](https://www.genome.jp/aaindex/)
"aa_feature"
