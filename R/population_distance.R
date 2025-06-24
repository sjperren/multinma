#' Comparing population covariates
#'
#' Runs distance tests between multivariate distributions across studies on different subnetworks.
#'
#' @param network An `nma_data` object.
#' @param covariates Character vector of covariate names.
#' @param scale Logical; whether to scale differences using pooled SDs.
#' @param method Character; one of "unscaled", "scaled", or "energy". Determines the distance calculation method.
#' @param binary_covariates Optional character vector of binary covariate names.
#'
#' @return A list with a `summary` dataframe and `distance_matrix`.
#' @export

population_distance <- function(network,
                                covariates = NULL,
                                method = c("unscaled", "scaled", "energy"),
                                binary_covariates = NULL) {

  # Check method argument
  method <- match.arg(method)

  # Map to old 'scale' logic for the Euclidean cases
  scale <- if (method == "scaled") TRUE else FALSE

  # Check network
  if (!inherits(network, "nma_data")) {
    abort("Expecting an `nma_data` object, as created by the functions `set_*`, `combine_network`, or `add_integration`.")
  }

  if (all(purrr::map_lgl(network, is.null))) {
    abort("Empty network.")
  }

  if (method == "energy") {
    if(nrow(network$agd_arm) > 0 && !inherits(network, "mlnmr_data")) {
    abort("Expecting a `mlnmr_data` object, as created by the function `add_integration` when using the `energy` method with aggregate data.")
    }
  }
  # Checks for covariates argument
  if (is.null(covariates)) {
    abort("`covariates` argument must be specified and cannot be NULL.")
  }

  if (!is.character(covariates)) {
    abort("`covariates` argument must be a character vector of covariate names.")
  }

  if (length(covariates) == 0) {
    abort("`covariates` argument must not be an empty character vector.")
  }

  if (!is.null(binary_covariates) && !all(binary_covariates %in% covariates)) {
    abort("All `binary_covariates` must also be included in `covariates`.")
  }
  # Check IPD covariates
  ipd_covariates <- colnames(network$ipd)
  missing_ipd_covariates <- setdiff(covariates, ipd_covariates)
  if (length(missing_ipd_covariates) > 0) {
    abort(("The following covariates are missing from the IPD data: {paste(missing_ipd_covariates, collapse=', ')}"))
  }

  # Check AGD covariates (exact or with "_mean" suffix)
  if (method == "energy" && nrow(network$agd_arm) > 0) {
    missing_covs <- covariates[!covariates %in% network$int_names]
    if (length(missing_covs) > 0) {
      abort(glue::glue("The following covariates are not in `network$int_names`: {paste(missing_covs, collapse = ', ')}. ",
                       "When using aggregate data, ensure all covariates have been included via the `add_integration()` function."))
    }
  }
  agd_covariates <- c(colnames(network$agd_contrast), colnames(network$agd_arm))
  agd_covariates_base <- unique(c(
    agd_covariates,
    sub("_mean$", "", agd_covariates)
  ))

  missing_agd_covariates <- setdiff(covariates, agd_covariates_base)
  if (length(missing_agd_covariates) > 0) {
    abort(("The following covariates are missing from the AGD data: {paste(missing_agd_covariates, collapse=', ')}"))
  }

  # Process IPD: Convert logical to numeric and average by study
  if (method != "energy") {
  ipd_covariate_data <- network$ipd
  ipd_covariate_data[covariates] <- lapply(ipd_covariate_data[covariates], function(x) {
    if (is.logical(x)) as.numeric(x) else x
  })
  ipd_summary <- ipd_covariate_data |>
    dplyr::group_by(.study) |>
    dplyr::summarise(
      dplyr::across(all_of(covariates), list(mean = ~ mean(.x, na.rm = TRUE),
                                             sd = ~ sd(.x, na.rm = TRUE))),
      .groups = "drop"
    )
  # AgD function to get means
  extract_agd_means <- function(agd_df, binary_covariates = NULL) {
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

      # If scaled, handle SDs
      if (scale && add_covariate) {
        if (sd_col %in% colnames(agd_df)) {
          df[[paste0(cov, "_sd")]] <- agd_df[[sd_col]]
        } else if (!is.null(binary_covariates) && cov %in% binary_covariates) {
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
  agd_contrast_means <- extract_agd_means(network$agd_contrast, binary_covariates)
  agd_arm_means <- extract_agd_means(network$agd_arm, binary_covariates)

  agd_all <- dplyr::bind_rows(agd_contrast_means, agd_arm_means)

  # Weighted summarisation for AGD
  agd_summary <- agd_all |>
    dplyr::group_by(.study) |>
    dplyr::summarise(
      total_n = sum(.data$.sample_size, na.rm = TRUE),
      !!!setNames(
        unlist(
          lapply(covariates, function(cov) {
            if (scale) {
              list(
                rlang::expr(
                  weighted.mean(!!rlang::sym(paste0(cov, "_mean")), w = .data$.sample_size, na.rm = TRUE)
                ),
                rlang::expr(
                  sqrt(weighted.mean((!!rlang::sym(paste0(cov, "_sd")))^2, w = .data$.sample_size, na.rm = TRUE))
                )
              )
            } else {
              list(
                rlang::expr(
                  weighted.mean(!!rlang::sym(paste0(cov, "_mean")), w = .data$.sample_size, na.rm = TRUE)
                )
              )
            }
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

  ipd_summary$source <- "IPD"
  agd_summary$source <- "AGD"
  all_summary <- dplyr::bind_rows(ipd_summary, agd_summary)
  g <- igraph::as.igraph(network)
  components <- igraph::components(g)

  treatment_components <- data.frame(
    .trt = names(components$membership),
    subnetwork = components$membership
  )

  study_trt_lookup <- list(
    network$ipd,
    network$agd_contrast,
    network$agd_arm
  ) |>
    purrr::compact() |>
    purrr::map_dfr(~ {
      cols <- colnames(.x)
      if (all(c(".study", ".trt") %in% cols)) {
        dplyr::tibble(
          .study = as.character(.x$.study),
          .trt   = as.character(.x$.trt)
        )
      } else {
        NULL
      }
    }) |>
    dplyr::distinct()

  # Join subnetwork info to each study
  study_components <- study_trt_lookup |>
    dplyr::left_join(treatment_components, by = ".trt") |>
    dplyr::select(-.trt) |>                      # Drop the treatment column
    dplyr::distinct(.study, subnetwork)          # Keep one row per study

  all_summary <- dplyr::left_join(all_summary, study_components, by = ".study")

  # Split the data by subnetwork
  sub1 <- dplyr::filter(all_summary, subnetwork == 1)
  sub2 <- dplyr::filter(all_summary, subnetwork == 2)

  # Prepare empty matrix to store distances
  dist_matrix <- matrix(NA,
                        nrow = nrow(sub1),
                        ncol = nrow(sub2),
                        dimnames = list(sub1$.study, sub2$.study))

  dist_matrix_full <- matrix(NA,
                             nrow = nrow(all_summary),
                             ncol = nrow(all_summary),
                             dimnames = list(all_summary$.study, all_summary$.study))

  # Calculate distances (scale or unscale)
  # Sub network 1  VS Sub network 2
  for (i in seq_len(nrow(sub1))) {
    for (j in seq_len(nrow(sub2))) {

      # Extract covariate means
      vec1 <- as.numeric(sub1[i, paste0(covariates, "_mean")])
      vec2 <- as.numeric(sub2[j, paste0(covariates, "_mean")])

      if (scale) {
        # Extract SDs and compute pooled SDs
        sd1 <- as.numeric(sub1[i, paste0(covariates, "_sd")])
        sd2 <- as.numeric(sub2[j, paste0(covariates, "_sd")])
        pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)

        # Avoid division by zero or NA
        valid <- !is.na(vec1) & !is.na(vec2) & !is.na(pooled_sd) & pooled_sd > 0
        diff_scaled <- (vec1[valid] - vec2[valid]) / pooled_sd[valid]
        dist_matrix[i, j] <- sqrt(sum(diff_scaled^2))
      } else {
        # Unscaled Euclidean distance
        valid <- !is.na(vec1) & !is.na(vec2)
        dist_matrix[i, j] <- sqrt(sum((vec1[valid] - vec2[valid])^2))
      }
    }
  }

  # All studies vs all studies
  for (i in seq_len(nrow(all_summary))) {
    for (j in seq_len(nrow(all_summary))) {
      if (i == j) next   # leave diagonal as NA

      # Extract covariate means
      vec1 <- as.numeric(all_summary[i, paste0(covariates, "_mean")])
      vec2 <- as.numeric(all_summary[j, paste0(covariates, "_mean")])

      if (scale) {
        # Extract SDs and compute pooled SDs
        sd1 <- as.numeric(all_summary[i, paste0(covariates, "_sd")])
        sd2 <- as.numeric(all_summary[j, paste0(covariates, "_sd")])
        pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)

        # Avoid division by zero or NA
        valid <- !is.na(vec1) & !is.na(vec2) & !is.na(pooled_sd) & pooled_sd > 0
        diff_scaled <- (vec1[valid] - vec2[valid]) / pooled_sd[valid]
        dist_matrix_full[i, j] <- sqrt(sum(diff_scaled^2))
      } else {
        # Unscaled Euclidean distance
        valid <- !is.na(vec1) & !is.na(vec2)
        dist_matrix_full[i, j] <- sqrt(sum((vec1[valid] - vec2[valid])^2))
      }
    }
  }

  return(list(
    summary = all_summary,
    distance_matrix = dist_matrix,
    distance_matrix_full = dist_matrix_full
  ))
  }

  if (method == "energy") {
    # Energy distance calculation
  library(energy)
    ipd_data <- network$ipd %>%
      dplyr::select(.study, dplyr::all_of(covariates)) %>%
      dplyr::group_split(.study)

    # Assign names
    names(ipd_data) <- network$ipd %>%
      dplyr::distinct(.study) %>%
      dplyr::pull(.study)

    # Drop .study column from each tibble
    ipd_data <- lapply(ipd_data, function(df) df[, setdiff(names(df), ".study"), drop = FALSE])

    # Subset to .study and .int_* columns
    int_cols <- grep("^\\.int_", names(network$agd_arm), value = TRUE)
    agd_arm_subset <- network$agd_arm[, c(".study", int_cols)]

    # Split by study, then remove the `.study` column from each
    agd_by_study <- lapply(split(agd_arm_subset, agd_arm_subset$.study), function(df) {
      # Remove .study column
      df <- df[ , int_cols, drop = FALSE]

      # For each .int_* column (which is a list of vectors), flatten into one vector
      covariate_data <- lapply(df, function(var_list) {
        unlist(var_list, use.names = FALSE)  # combine all arms into a single vector
      })

      covariate_df <- as.data.frame(covariate_data)
      # Rename columns: remove ".int_" prefix
      names(covariate_df) <- gsub("^\\.int_", "", names(covariate_df))

      return(covariate_df)
    })

    # Remove any studies in agd_by_study that are NULL or have 0 rows
    agd_by_study <- Filter(function(df) {
      !is.null(df) && nrow(df) > 0
    }, agd_by_study)

    all_studies <- c(ipd_data, agd_by_study)

    #Generate all unique study‐pairs
    study_names <- names(all_studies)
    pairs <- combn(study_names, 2, simplify = FALSE)

    # Loop through pairs, bind & test
    energy_results <- lapply(pairs, function(p) {
      # Pull and tag each study’s data
      df1 <- all_studies[[p[1]]] %>% mutate(.study = p[1])
      df2 <- all_studies[[p[2]]] %>% mutate(.study = p[2])

      # Record group sizes
      sizes <- c(nrow(df1), nrow(df2))

      # Combine and compute distance on covariates only
      X <- bind_rows(df1, df2) %>% select(-.study) %>% as.matrix()
      d <- dist(X)

      # Run the energy test
      set.seed(1234)
      et <- eqdist.etest(d, sizes = sizes, distance = TRUE, R = 199)

      # Return a one‐row summary
      data.frame(
        study1    = p[1],
        study2    = p[2],
        size1     = sizes[1],
        size2     = sizes[2],
        statistic = et$statistic,
        p.value   = et$p.value,
        row.names = NULL
      )
    })

    # 4. Tidy up into a single data.frame
    summary_results <- bind_rows(energy_results)
    class(summary_results) <- c("energy_distance_results", class(summary_results))

    return(summary_results)
  }
}
#' Plot distance matrix from `population_distance()`
#'
#' Produces a coloured `gt` table of distances from `population_distance()`.
#'
#' @param x Output from `alm()`.
#'
#' @return A `gt` table.
#' @export

# Create a coloured table from the output of alm()
plot_alm_matrix <- function(x) {
  if (inherits(x, "energy_distance_results")) {
    studies <- sort(unique(c(x$study1, x$study2)))
    n       <- length(studies)

    # 2. Init empty matrices
    stat_mat <- matrix(NA,   n, n, dimnames = list(studies, studies))
    pval_mat <- matrix(NA,   n, n, dimnames = list(studies, studies))

    # look up numeric positions of study1 & study2 in your ordering
    idx1 <- match(x$study1, studies)
    idx2 <- match(x$study2, studies)

    for(i in seq_along(idx1)) {
      # force the larger index into the ROW, the smaller into the COL
      r <- max(idx1[i], idx2[i])
      c <- min(idx1[i], idx2[i])
      stat_mat[r, c] <- x$statistic[i]
      pval_mat[r, c] <- x$p.value[i]
    }

    # 4. Combine into one character matrix
    combo <- matrix("", n, n, dimnames = list(studies, studies))
    for(i in seq_len(n)) for(j in seq_len(n)) {
      if (!is.na(stat_mat[i,j])) {
        combo[i,j] <- sprintf("%.2f\n(p=%.3f)",
                              stat_mat[i,j],
                              pval_mat[i,j])
      }
    }

    # — 5. Fill diagonal with each study’s N —
    sizes_by_study <- sapply(studies, function(s) {
      # pick up the unique size1 or size2 for that study
      sz1 <- unique(x$size1[x$study1 == s])
      sz2 <- unique(x$size2[x$study2 == s])
      sz  <- unique(c(sz1, sz2))
      if (length(sz) > 1) warning("Multiple sizes for ", s, ": ", paste(sz, collapse=","))
      sz[1]
    })
    diag(combo) <- as.character(sizes_by_study)

    # 5. Blank out diagonal & upper triangle
    combo[upper.tri(combo)] <- ""

    tbl <- gridExtra::tableGrob(
      combo,
      rows = studies,
      theme = gridExtra::ttheme_minimal(
        core = list(fg_params = list(fontsize=10)),
        colhead = list(fg_params = list(fontsize=11, fontface="bold"))
      )
    )
    grid::grid.newpage()
    grid::grid.draw(tbl)

    invisible(tbl)
  } else {
  mat <- x$distance_matrix
  mat_full <- x$distance_matrix_full
  summary_df <- x$summary

  # Convert matrix to dataframe with .study column
  df <- as.data.frame(mat)
  df <- tibble::rownames_to_column(df, var = ".study")

  df_full <- as.data.frame(mat_full)
  df_full <- tibble::rownames_to_column(df_full, var = ".study")

  # Extract correct data sources for rows and columns using the 'source' column
  row_colors <- summary_df |>
    dplyr::filter(.study %in% rownames(mat)) |>
    dplyr::select(.study, source)

  col_colors <- summary_df |>
    dplyr::filter(.study %in% colnames(mat)) |>
    dplyr::select(.study, source)

  row_colors_full <- summary_df |>
    dplyr::filter(.study %in% rownames(mat_full)) |>
    dplyr::select(.study, source)

  col_colors_full <- summary_df |>
    dplyr::filter(.study %in% colnames(mat_full)) |>
    dplyr::select(.study, source)

  # Create gt table s1 vs s2
  gt_tbl <- gt::gt(df, rowname_col = ".study") |>
    gt::data_color(
      columns = everything(),
      fn = scales::col_numeric(
        palette = c("white", "#FDB0B0", "red"),
        domain = range(mat, na.rm = TRUE)
      )
    ) |>
    # Style row labels (stub) using IPD/AGD color
    gt::tab_style(
      style = list(cell_fill(color = "navy"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_stub(rows = df$.study %in% row_colors$.study[row_colors$source == "AGD"])
    ) |>
    gt::tab_style(
      style = list(cell_fill(color = "darkgreen"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_stub(rows = df$.study %in% row_colors$.study[row_colors$source == "IPD"])
    ) |>
    # Style column headers using IPD/AGD color
    gt::tab_style(
      style = list(cell_fill(color = "navy"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_column_labels(columns = col_colors$.study[col_colors$source == "AGD"])
    ) |>
    gt::tab_style(
      style = list(cell_fill(color = "darkgreen"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_column_labels(columns = col_colors$.study[col_colors$source == "IPD"])
    )

  # Create gt table for FULL
  gt_tbl_full <- gt::gt(df_full, rowname_col = ".study") |>
    gt::data_color(
      columns = everything(),
      fn = scales::col_numeric(
        palette = c("white", "#FDB0B0", "red"),
        domain = range(mat_full, na.rm = TRUE)
      )
    ) |>
    # Style row labels (stub) using IPD/AGD color
    gt::tab_style(
      style = list(cell_fill(color = "navy"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_stub(rows = df_full$.study %in% row_colors_full$.study[row_colors_full$source == "AGD"])
    ) |>
    gt::tab_style(
      style = list(cell_fill(color = "darkgreen"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_stub(rows = df_full$.study %in% row_colors_full$.study[row_colors_full$source == "IPD"])
    ) |>
    # Style column headers using IPD/AGD color
    gt::tab_style(
      style = list(cell_fill(color = "navy"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_column_labels(columns = col_colors_full$.study[col_colors_full$source == "AGD"])
    ) |>
    gt::tab_style(
      style = list(cell_fill(color = "darkgreen"), cell_text(color = "white", weight = "bold")),
      locations = gt::cells_column_labels(columns = col_colors_full$.study[col_colors_full$source == "IPD"])
    )
  return(list(
    subnetwork_table = gt_tbl,
    full_table = gt_tbl_full)
  )
  }
}




#' @export
baseline_synthesis <- function(network,
                               consistency = c("consistency", "ume", "nodesplit"),
                               trt_effects = c("fixed", "random"),
                               regression = NULL,
                               class_interactions = c("common", "exchangeable", "independent"),
                               class_effects = c("independent", "common", "exchangeable"),
                               class_sd = c("independent", "common"),
                               likelihood = NULL, link = NULL,
                               ...,
                               nodesplit = get_nodesplits(network, include_consistency = TRUE),
                               prior_intercept = .default(normal(scale = 100)),
                               prior_trt = .default(normal(scale = 10)),
                               prior_het = .default(half_normal(scale = 5)),
                               prior_het_type = c("sd", "var", "prec"),
                               prior_reg = .default(normal(scale = 10)),
                               prior_aux = .default(),
                               prior_aux_reg = .default(),
                               prior_class_mean = .default(normal(scale = 10)),
                               prior_class_sd = .default(half_normal(scale = 5)),
                               aux_by = NULL,
                               aux_regression = NULL,
                               QR = FALSE,
                               center = TRUE,
                               adapt_delta = NULL,
                               int_thin = 0,
                               int_check = TRUE,
                               mspline_degree = 3,
                               n_knots = 7,
                               knots = NULL,
                               mspline_basis = NULL,
                               prior_intercept_sd = .default(half_normal(scale = 5)),
                               random_baseline = TRUE) {

  if (.is_default(prior_intercept_sd)) {
    warn(glue::glue("Warning: 'prior_intercept_sd' was left at its default value: {get_prior_call(prior_intercept_sd)}")) }

  check_prior(prior_intercept_sd)

  fit <- nma(network = network,
             consistency = consistency,
             trt_effects = trt_effects,
             regression = regression,
             class_interactions = class_interactions,
             class_effects = class_effects,
             class_sd = class_sd,
             likelihood = likelihood,
             link = link,
             ...,
             nodesplit = nodesplit,
             prior_intercept = prior_intercept,
             prior_trt = prior_trt,
             prior_het = prior_het,
             prior_het_type = prior_het_type,
             prior_reg = prior_reg,
             prior_aux = prior_aux,
             prior_aux_reg = prior_aux_reg,
             prior_class_mean = prior_class_mean,
             prior_class_sd = prior_class_sd,
             aux_by = aux_by,
             aux_regression = aux_regression,
             QR = QR,
             center = center,
             adapt_delta = adapt_delta,
             int_thin = int_thin,
             int_check = int_check,
             mspline_degree = mspline_degree,
             n_knots = n_knots,
             knots = knots,
             mspline_basis = mspline_basis,
             random_baseline = random_baseline,
             prior_intercept_sd = prior_intercept_sd)

  ss <- rstan::summary(fit$stanfit,
                       pars  = c("baseline_new","baseline_mean","baseline_sd","mu"),
                       probs = c(0.025, 0.5, 0.975))$summary

  keep <- grepl("^(baseline_new|baseline_mean|baseline_sd|mu\\[)", rownames(ss))
  summary_df <- as.data.frame(ss[keep, , drop = FALSE])
  summary_df$parameter <- rownames(ss)[keep]
  summary_df <- summary_df[, c("parameter", setdiff(names(summary_df), "parameter"))]
  rownames(summary_df) <- NULL

  fit$baseline_summary <- summary_df

  class(fit) <- c("baseline_synthesis", class(fit))
  fit
}

#' @export
print.baseline_synthesis <- function(x, ...) {
  print(x$baseline_summary)
  invisible(x)
}


