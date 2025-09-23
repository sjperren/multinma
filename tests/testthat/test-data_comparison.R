test_that("arguments entered into population_distance() get correct errors and warnings", {
  m <- ("Expecting an `nma_data` object, as created by the functions `set_\\*`, `combine_network`, or `add_integration`\\.")
  n <- ("`covariates` argument must be specified and cannot be NULL.")
  p <- ("All `binary` must also be included in `covariates`.")
  af_net <- set_agd_arm(atrial_fibrillation[atrial_fibrillation$studyc != "WASPO", ],
                        study = studyc,
                        trt = trtc,
                        r = r,
                        n = n,
                        trt_class = trt_class)

  expect_error(population_distance(atrial_fibrillation, method = "alm"), m)
  expect_error(population_distance(1, method = "alm"), m)
  expect_error(population_distance(af_net, covariates = c(), method = "alm"), n)
  expect_error(population_distance(af_net, covariates = , method = "alm"), n)
  expect_error(population_distance(af_net, covariates = c("strok"), binary = c("stroke"), method = "alm"), p)
}
)
test_that("checking correct errors when missing covariates in data", {
  # --- Dummy data -----------------------------------------
  ipd_df <- rbind(
    data.frame(study="S1", trt="A",
               a = rnorm(30, 50, 10),
               b = rbinom(30, 1, 0.4),
               y = rnorm(30, 0.0, 1.0)),
    data.frame(study="S1", trt="B",
               a = rnorm(28, 52, 10),
               b = rbinom(28, 1, 0.6),
               y = rnorm(28, 0.4, 1.0))
  )

  agd_arm_df <- data.frame(
    study=c("S2","S2"), trt=c("B","C"), n=c(60,55),
    a=c(51,49), a_sd=c(1,1), c=c(10,12), c_sd=c(1,1), y=c(0.3,0.1), se=c(0.20,0.22)
  )

  agd_contrast_df <- data.frame(
    study=c("S3","S3"), trt=c("C","A"), n=c(50,48),
    b=c(0.45,0.40), b_sd=c(1,1), c=c(11,11.5), c_sd=c(1,1),
    y=c(NA, 0.25), se=c(0.25, 0.30)
  )

  net <- combine_network(
    set_ipd(ipd_df, study = study, trt = trt, y = y, trt_ref = "A"),
    set_agd_arm(agd_arm_df, study = study, trt = trt, y = y, se = se, sample_size = n),
    set_agd_contrast(agd_contrast_df, study = study, trt = trt, y = y, se = se, sample_size = n),
    trt_ref = "A"
  )

  net_ipd_agd_arm <- combine_network(
    set_ipd(ipd_df, study = study, trt = trt, y = y, trt_ref = "A"),
    set_agd_arm(agd_arm_df, study = study, trt = trt, y = y, se = se),
    trt_ref = "A"
  )

  net_agd_con <- set_agd_contrast(agd_contrast_df, study = study, trt = trt, y = y, se = se)

  expect_error(population_distance(net, method = "alm", covariates = c("b")), "The following covariates are missing from the AGD arm data: b")
  expect_error(population_distance(net, method = "alm", covariates = c("c")), "The following covariates are missing from the IPD data: c")
  expect_error(population_distance(net, method = "alm", covariates = c("a")), "The following covariates are missing from the AGD contrast data: a")

  expect_error(population_distance(net_ipd_agd_arm, method = "alm", covariates = c("a")), "Aggregate arm data must contain a '.sample_size' column.")
  expect_error(population_distance(net_agd_con, method = "alm", covariates = c("b")), "Aggregate contrast data must contain a '.sample_size' column.")

  # Need to write error checks for when data contains NAs or agd doesn't include mean or sd columns of variables.

  }
)


