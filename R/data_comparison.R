#' Comparing population covariates
#'
#' Runs distance tests between multivariate distributions across studies.
#'
#' @param network An `nma_data` object.
#' @param covariates Character vector of covariate names. Cannot except categorical variables
#' @param method Character; "alm", or "energy". Determines the distance calculation method.
#' @param binary If AgD is within the network, binary variables must be stated.
#'
#' @return A list with a `summary` dataframe and `distance_matrix`.
#' @export

population_distance <- function(network,
                                covariates = NULL,
                                method = c("alm", "energy"),
                                binary = NULL) {
  # Check method argument
  method <- match.arg(method)

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }
  # Checks for covariates argument
  if (is.null(covariates)) {
    abort("`covariates` argument must be specified and cannot be NULL.")
  }
  if (!is.character(covariates)) {
    abort("`covariates` argument must be a character vector of covariate names.")
  }
  # checking binary variables are within covariates
  if (!is.null(binary) && !all(binary %in% covariates)) {
    abort("All `binary` must also be included in `covariates`.")
  }
  # Check IPD covariates
  if (nrow(network$ipd) > 0) {
    ipd_covariates <- colnames(network$ipd)
    missing_ipd_covariates <- setdiff(covariates, ipd_covariates)
    if (length(missing_ipd_covariates) > 0) {
      abort(("The following covariates are missing from the IPD data: {paste(missing_ipd_covariates, collapse=', ')}"))
    }
  }
  # Check AGD covariates (exact or with "_mean" suffix)
  if (nrow(network$agd_arm) > 0) {
    # Ensure .sample_size exists
    if (!".sample_size" %in% colnames(network$agd_arm)) {
      abort("Aggregate arm data must contain a '.sample_size' column.")
    }
    agd_covariates <- c(colnames(network$agd_arm))
    agd_covariates <- unique(c(agd_covariates, sub("_mean$", "", agd_covariates)))
    missing_agd_covariates <- setdiff(covariates, agd_covariates)
    if (length(missing_agd_covariates) > 0) {
      abort(("The following covariates are missing from the AGD data: {paste(missing_agd_covariates, collapse=', ')}"))
    }
  }

  # Check AGD contrast covariates (exact or with "_mean" suffix)
  if (nrow(network$agd_contrast) > 0) {
    # Ensure .sample_size exists
    if (!".sample_size" %in% colnames(network$agd_contrast)) {
      abort("Aggregate contrast data must contain a '.sample_size' column.")
    }
    agd_covariates <- c(colnames(network$agd_contrast))
    agd_covariates <- unique(c(agd_covariates, sub("_mean$", "", agd_covariates)))
    missing_agd_covariates <- setdiff(covariates, agd_covariates)
    if (length(missing_agd_covariates) > 0) {
      abort(("The following covariates are missing from the AGD contrast data: {paste(missing_agd_covariates, collapse=', ')}"))
    }
  }

  # Process IPD: Convert logical to numeric and average by study
  if (method != "energy") {
    if (nrow(network$ipd) > 0) {
    ipd_covariate_data <- network$ipd
    ipd_covariate_data[covariates] <- lapply(ipd_covariate_data[covariates], function(x) {
      if (is.logical(x)) as.numeric(x) else x
    })
    ipd_summary <- ipd_covariate_data %>%
      dplyr::group_by(.study) %>%
      dplyr::summarise(
        dplyr::across(all_of(covariates), list(mean = ~ mean(.x, na.rm = TRUE),
                                               sd = ~ sd(.x, na.rm = TRUE))),
        .groups = "drop"
      )
    }
  # AgD function to get means
  extract_agd_means <- function(agd_df, binary = NULL) {
    if (nrow(agd_df) == 0) return(NULL)

    df <- agd_df[, c(".study", ".sample_size"), drop = FALSE]
    retained_covariates <- c()

    for (cov in covariates) {
      mean_col <- paste0(cov, "_mean")
      sd_col <- paste0(cov, "_sd")
      add_covariate <- TRUE

      # Add mean if available
      if (mean_col %in% colnames(agd_df)) {
        df[[paste0(cov, "_mean")]] <- agd_df[[mean_col]]
      } else if (cov %in% colnames(agd_df)) {
        df[[paste0(cov, "_mean")]] <- agd_df[[cov]]
      } else {
        warning(glue::glue("Mean for covariate '{cov}' not found in AgD. Covariate dropped."))
        add_covariate <- FALSE
      }

      # Add SDs if available
      if (add_covariate) {
        if (sd_col %in% colnames(agd_df)) {
          df[[paste0(cov, "_sd")]] <- agd_df[[sd_col]]
        } else if (!is.null(binary) && cov %in% binary) {
          p <- df[[paste0(cov, "_mean")]]
          df[[paste0(cov, "_sd")]] <- sqrt(p * (1 - p))
        } else {
          warning(glue::glue("Missing SD for covariate '{cov}' in AgD and not marked as binary. Covariate dropped."))
          df[[paste0(cov, "_mean")]] <- NULL
          add_covariate <- FALSE
        }
      }
      if (add_covariate) {
        retained_covariates <- c(retained_covariates, cov)
      } else {
        # Remove if SD step failed
        df[[paste0(cov, "_mean")]] <- NULL
        df[[paste0(cov, "_sd")]] <- NULL
      }
    }

    # Update the covariates list inside alm(), if needed
    assign("covariates", retained_covariates, envir = parent.env(environment()))

    return(df)
  }
  agd_contrast_means <- extract_agd_means(network$agd_contrast, binary)
  agd_arm_means <- extract_agd_means(network$agd_arm, binary)

  agd_all <- dplyr::bind_rows(agd_contrast_means, agd_arm_means)

  # Weighted summarisation for AGD
  agd_summary <- agd_all %>%
    dplyr::group_by(.study) %>%
    dplyr::summarise(
      total_n = sum(.data$.sample_size, na.rm = TRUE),
      !!!setNames(
        unlist(
          lapply(covariates, function(cov) {
              list(
                rlang::expr(
                  weighted.mean(!!rlang::sym(paste0(cov, "_mean")), w = .data$.sample_size, na.rm = TRUE)
                ),
                rlang::expr(
                  sqrt(weighted.mean((!!rlang::sym(paste0(cov, "_sd")))^2, w = .data$.sample_size, na.rm = TRUE))
                )
              )
              list(
                rlang::expr(
                  weighted.mean(!!rlang::sym(paste0(cov, "_mean")), w = .data$.sample_size, na.rm = TRUE)
                )
              )
          }),
          recursive = FALSE
        ),
        # Clean column names
        unlist(
          lapply(covariates, function(cov) {
            if (scale) {
              c(paste0(cov, "_mean"), paste0(cov, "_sd"))
            } else {
              paste0(cov, "_mean")
            }
          })
        )
      ),
      .groups = "drop"
    )

  stop_point <- "whatever"
  }
}

