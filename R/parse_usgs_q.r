#' Parse USGS NGDB qvalue + qualifier into censoring fields
#'
#' @param qvalue Numeric or character vector of qvalues.
#' @param qualifier Optional character vector (e.g. "<", "L", ">", "G").
#'
#' @details
#' Rules implemented (USGS NGDB style):
#' * qvalue < 0  → left-censored at LOD = abs(qvalue)
#' * qvalue text ending ".1111" or ".11111" → right-censored at upper limit, numeric part
#' * qualifier "<", "L", or "N" → left-censored, limit from numeric qvalue
#' * qualifier ">" or "G" → right-censored, limit from numeric qvalue
#'
#' @return A tibble with columns:
#' * detected_val
#' * left_cen, left_limit
#' * right_cen, right_limit
#' Invalid values, unknown qualifiers, and conflicting censoring encodings are
#' excluded from numerical fields and identified in `validation_status`.
#' Only blank/missing, `L`, `N`, `<`, `G`, and `>` qualifiers are accepted.
#' @export
parse_usgs_q <- function(qvalue, qualifier = NULL) {
  q_chr <- stringr::str_trim(as.character(qvalue))
  
  # suffix-based upper limits
  is_upper_suffix <- grepl("\\.1{4,5}$", q_chr)
  upper_num <- suppressWarnings(as.numeric(sub("\\.1{4,5}$", "", q_chr)))
  
  # numeric form
  q_num <- suppressWarnings(as.numeric(q_chr))
  
  left_by_num <- is.finite(q_num) & q_num < 0
  lod_num     <- abs(q_num)
  
  if (is.null(qualifier)) {
    qualifier <- rep(NA_character_, length(qvalue))
  }
  if (length(qualifier) != length(qvalue)) {
    stop("parse_usgs_q: qualifier must be NULL or the same length as qvalue")
  }
  qualifier <- stringr::str_to_upper(stringr::str_trim(as.character(qualifier)))
  
  qual_left  <- qualifier %in% c("L", "N", "<")
  qual_right <- qualifier %in% c("G", ">")
  
  left_cen  <- left_by_num | qual_left
  right_cen <- is_upper_suffix | qual_right
  
  left_limit  <- ifelse(left_by_num, lod_num,
                        ifelse(qual_left & is.finite(q_num), abs(q_num), NA_real_))
  right_limit <- ifelse(is_upper_suffix, upper_num,
                        ifelse(qual_right & is.finite(q_num), q_num, NA_real_))
  
  both_cen <- left_cen & right_cen
  known_qualifier <- is.na(qualifier) | qualifier %in% c("", "L", "N", "<", "G", ">")
  status <- dplyr::case_when(
    both_cen ~ "conflicting_qualifiers",
    is.na(q_chr) | q_chr == "" ~ "missing_value",
    !is.finite(q_num) ~ "invalid_value",
    !known_qualifier ~ "unknown_qualifier",
    left_cen & (!is.finite(left_limit) | left_limit <= 0) ~ "invalid_limit",
    right_cen & (!is.finite(right_limit) | right_limit <= 0) ~ "invalid_limit",
    !left_cen & !right_cen & q_num <= 0 ~ "nonpositive_detection",
    TRUE ~ "usable"
  )
  invalid <- status != "usable"
  left_cen[invalid] <- FALSE
  right_cen[invalid] <- FALSE
  left_limit[invalid] <- NA_real_
  right_limit[invalid] <- NA_real_
  detected_val <- ifelse(left_cen | right_cen | invalid, NA_real_, q_num)

  tibble::tibble(
    detected_val = detected_val,
    left_cen     = left_cen,
    left_limit   = left_limit,
    right_cen    = right_cen,
    right_limit  = right_limit,
    validation_status = status
  )
}
