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

test_that("fixed baseline connection merges studies", {
  spec <- con("fixed", studies = levels(net$studies)[1:2])
  fit <- nma(net, connect_baseline = spec,
             prior_intercept = normal(scale = 10),
             prior_trt = normal(scale = 10),
             iter = 1, chains = 1)
  expect_true(any(grepl("&", levels(fit$network$studies))))
})

test_that("fixed baseline connection not allowed for contrast data", {
  contr_net <- set_agd_contrast(parkinsons,
                                study = studyn,
                                trt = trtn,
                                y = diff,
                                se = se_diff,
                                sample_size = n)
  spec <- con("fixed", studies = levels(contr_net$studies)[1:2])
  expect_error(
    nma(contr_net, connect_baseline = spec,
        prior_intercept = normal(scale = 10),
        prior_trt = normal(scale = 10),
        iter = 1, chains = 1),
    "AgD contrast data"
  )
})
