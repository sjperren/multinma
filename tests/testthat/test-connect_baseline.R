library(multinma)

# small example network
data(smoking)
net <- set_agd_arm(smoking, study = studyn, trt = trtc, r = r, n = n,
                   trt_ref = "no")


test_that("con helper validation", {
  expect_s3_class(con("fixed", studies = "1"), "nma_connect")
  expect_error(con("random", studies = "1"), "baseline_prior")
})


cb_spec <- con("random", studies = levels(net$studies)[1:2],
               baseline_prior = normal(scale = 10))

fixed_spec <- con("fixed", studies = levels(net$studies)[3])


test_that("nma connects baselines", {
  expect_error(
    nma(net, connect_baseline = cb_spec,
        prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10),
        iter = 1, chains = 1),
    NA
  )
  expect_error(
    nma(net, connect_baseline = list(cb_spec, fixed_spec),
        prior_intercept = normal(scale = 10), prior_trt = normal(scale = 10),
        iter = 1, chains = 1),
    NA
  )
})
