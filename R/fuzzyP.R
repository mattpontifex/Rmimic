#' fuzzyP
#'
#' @description Function to properly round p values. Returns several values allowing the user to choose what rounded value is appropriate.
#'
#' @param datain Significance value.
#' @param alpha Critical value for alpha.
#' @param html Boolean parameter for output in html format.
#'
#' @return A list with the following elements:
#' \item{interpret}{Rounded p value for use in interpretations.}
#' \item{report}{Rounded p value for use in summary reports.}
#' \item{exact}{Exact p value.}
#' \item{modifier}{Equals sign unless the p value is less than .001.}
#' \item{significance}{Boolean for if the value is less than or equal to the critical alpha.}
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, July 28, 2017
#'
#' @examples
#'     # Round a randomly generated number between 0 and 1.
#'     outPvalue <- fuzzyP(stats::runif(1))
#'
#' @export

fuzzyP <- function(datain, studywiseAlpha = 0.05, html = FALSE) {
  # --- helpers ---------------------------------------------------------------
  is_scalar_numeric <- function(x) is.numeric(x) && length(x) == 1L
  # count decimals in alpha after removing trailing zeros
  decimals_of <- function(x) {
    sx <- sub("0+$", "", as.character(x))
    parts <- strsplit(sx, ".", fixed = TRUE)[[1]]
    if (length(parts) == 2L) nchar(parts[2]) else 0L
  }
  fmt_fixed <- function(x, d) formatC(x, format = "f", digits = d)
  # trim trailing zeros but keep "0." from becoming "0"
  trim_trailing_zeros <- function(s) {
    if (is.na(s)) return(s)
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
    s
  }
  # make "< 0.001" threshold string & numeric for given digits d
  lower_str_num <- function(d) {
    # e.g., d=3 -> "0.001"
    str <- paste0("0.", paste(rep("0", max(0, d - 1L)), collapse = ""), "1")
    num <- as.numeric(str)
    list(str = str, num = num)
  }
  # choose symbols
  equalityparameters <- c('<', '=', '>')
  if (isTRUE(html)) {
    equalityparameters <- c('&lt;', '&#61;', '&gt;')
  }
  # robust compare epsilon (to dampen float wobbles)
  eps <- .Machine$double.eps * 10
  
  # --- validation ------------------------------------------------------------
  if (!is_scalar_numeric(studywiseAlpha) || !is.finite(studywiseAlpha)) {
    stop("studywiseAlpha must be a single finite numeric.")
  }
  if (!(studywiseAlpha > 0 && studywiseAlpha < 1)) {
    stop("studywiseAlpha must be in (0, 1).")
  }
  if (!is.logical(html) || length(html) != 1L) html <- FALSE
  
  # shape outputs for NA / non-numeric early
  if (!is_scalar_numeric(datain) || is.na(datain) || !is.finite(datain)) {
    return(list(
      interpret    = 1,
      report       = NA_character_,
      exact        = datain,
      modifier     = equalityparameters[2],
      significance = FALSE
    ))
  }
  
  # clamp to [0,1]
  p_exact <- datain
  if (p_exact < 0 - eps) {
    p_exact <- 0
  }
  if (p_exact > 1 + eps) {
    p_exact <- 1
  }
  p_exact <- max(0, min(1, p_exact))
  
  # --- rounding policy -----------------------------------
  # Base precision is tied to alpha's *meaningful* decimals, then bounded
  alpha_decimals <- decimals_of(studywiseAlpha)  # e.g., 0.050 -> 2
  d_interpret <- min(max(alpha_decimals, 2L), 6L)
  d_report <- min(d_interpret + 1L, 6L)
  
  # floor interpret if above alpha after rounding
  p_round <- round(p_exact, digits = d_interpret)
  if (p_exact > studywiseAlpha && p_round > studywiseAlpha) {
    # if it would round *above* alpha, force floor
    p_interpret_num <- floor(p_exact * 10^d_interpret) / 10^d_interpret
  } else {
    p_interpret_num <- p_round
  }
  
  report_str <- fmt_fixed(round(p_exact, digits = d_report), d_report)
  if (studywiseAlpha <= 0.05 + eps) {
    if (p_exact >= 0.3) {
      report_str <- fmt_fixed(round(p_exact, 1), 1)
    }
  }
  
  # extreme small p
  lb <- lower_str_num(d_report)
  modifier <- equalityparameters[2]
  if (p_exact + eps < lb$num) {
    report_str <- lb$str
    modifier <- equalityparameters[1]
  }
  
  # huge p near 1
  if (p_exact > 0.9 + eps) {
    report_str <- "0.9"
    modifier <- equalityparameters[3]
  }
  
  # significance test
  significance <- (p_interpret_num <= studywiseAlpha + eps)
  
  # final tidy of the printed report (trim trailing zeros sensibly)
  report_str <- trim_trailing_zeros(report_str)
  
  # guardrails: avoid printing impossible 0 or 1 snapshots
  if (!is.na(report_str)) {
    if (report_str %in% c("0", "0.", "0.0", "0.00", "0.000", "0.0000", "0.00000", "0.000000")) {
      report_str <- lower_str_num(max(3L, d_report))$str
      modifier <- equalityparameters[1]
    }
    if (report_str %in% c("1", "1.", "1.0", "1.00", "1.000", "1.0000", "1.00000", "1.000000")) {
      report_str <- "0.9"
      modifier <- equalityparameters[3]
    }
  }
  
  res <- list(
    interpret    = p_interpret_num,  # numeric
    report       = report_str,       # character
    exact        = datain,           # original input (pre-clamp), for transparency
    modifier     = modifier,         # "<", "=", ">", HTML-safe if requested
    significance = isTRUE(significance)
  )
  return(res)
}