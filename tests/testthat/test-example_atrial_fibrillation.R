# Generated by vignette example_atrial_fibrillation.Rmd: do not edit by hand
# Instead edit example_atrial_fibrillation.Rmd and then run precompile.R

skip_on_cran()



params <-
list(run_tests = FALSE)

## ----code=readLines("children/knitr_setup.R"), include=FALSE-----------------------------------------------------------------

## ----include=FALSE-----------------------------------------------------------------------------------------------------------
set.seed(4783982)


## ----eval = FALSE------------------------------------------------------------------------------------------------------------
# library(multinma)
# options(mc.cores = parallel::detectCores())

## ----setup, echo = FALSE-----------------------------------------------------------------------------------------------------
library(multinma)
nc <- switch(tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")), 
             "true" =, "warn" = 2, 
             parallel::detectCores())
options(mc.cores = nc)


## ----------------------------------------------------------------------------------------------------------------------------
head(atrial_fibrillation)


## ----------------------------------------------------------------------------------------------------------------------------
af_net <- set_agd_arm(atrial_fibrillation[atrial_fibrillation$studyc != "WASPO", ], 
                      study = studyc,
                      trt = trtc,
                      r = r, 
                      n = n,
                      trt_class = trt_class)
af_net


## ----af_network_plot, fig.width=8, fig.height=6, out.width="100%"------------------------------------------------------------
plot(af_net, weight_nodes = TRUE, weight_edges = TRUE, show_trt_class = TRUE) + 
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")


## ----------------------------------------------------------------------------------------------------------------------------
summary(normal(scale = 100))
summary(half_normal(scale = 5))


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
# af_fit_1 <- nma(af_net,
#                 trt_effects = "random",
#                 prior_intercept = normal(scale = 100),
#                 prior_trt = normal(scale = 100),
#                 prior_het = half_normal(scale = 5),
#                 adapt_delta = 0.99)

## ----echo=FALSE--------------------------------------------------------------------------------------------------------------
af_fit_1 <- nma(af_net, 
                seed = 103533305,
                trt_effects = "random",
                prior_intercept = normal(scale = 100),
                prior_trt = normal(scale = 100),
                prior_het = half_normal(scale = 5),
                adapt_delta = 0.99,
                iter = 5000)


## ----------------------------------------------------------------------------------------------------------------------------
af_fit_1


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
# # Not run
# print(af_fit_1, pars = c("d", "mu", "delta"))


## ----af_1_pp_plot, fig.width=8, fig.height=6, out.width="100%"---------------------------------------------------------------
plot_prior_posterior(af_fit_1, prior = c("trt", "het"))


## ----------------------------------------------------------------------------------------------------------------------------
(af_1_releff <- relative_effects(af_fit_1, trt_ref = "Placebo/Standard care"))


## ----af_1_releff_plot--------------------------------------------------------------------------------------------------------
plot(af_1_releff, ref_line = 0)


## ----af_1_ranks--------------------------------------------------------------------------------------------------------------
(af_1_ranks <- posterior_ranks(af_fit_1))
plot(af_1_ranks)

## ----af_1_rankprobs----------------------------------------------------------------------------------------------------------
(af_1_rankprobs <- posterior_rank_probs(af_fit_1))
plot(af_1_rankprobs)

## ----af_1_cumrankprobs-------------------------------------------------------------------------------------------------------
(af_1_cumrankprobs <- posterior_rank_probs(af_fit_1, cumulative = TRUE))
plot(af_1_cumrankprobs)


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
# af_fit_4b <- nma(af_net,
#                  trt_effects = "random",
#                  regression = ~ .trt:stroke,
#                  class_interactions = "common",
#                  QR = TRUE,
#                  prior_intercept = normal(scale = 100),
#                  prior_trt = normal(scale = 100),
#                  prior_reg = normal(scale = 100),
#                  prior_het = half_normal(scale = 5),
#                  adapt_delta = 0.99)

## ----echo=FALSE, eval=!params$run_tests--------------------------------------------------------------------------------------
# af_fit_4b <- nma(af_net,
#                  seed = 579212814,
#                  trt_effects = "random",
#                  regression = ~ .trt:stroke,
#                  class_interactions = "common",
#                  QR = TRUE,
#                  prior_intercept = normal(scale = 100),
#                  prior_trt = normal(scale = 100),
#                  prior_reg = normal(scale = 100),
#                  prior_het = half_normal(scale = 5),
#                  adapt_delta = 0.99)

## ----echo=FALSE, eval=params$run_tests---------------------------------------------------------------------------------------
af_fit_4b <- nowarn_on_ci(nma(af_net, 
                 seed = 579212814,
                 trt_effects = "random",
                 regression = ~ .trt:stroke,
                 class_interactions = "common",
                 QR = TRUE,
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100),
                 prior_reg = normal(scale = 100),
                 prior_het = half_normal(scale = 5),
                 adapt_delta = 0.99,
                 iter = 5000))


## ----------------------------------------------------------------------------------------------------------------------------
af_fit_4b


## ----eval=FALSE--------------------------------------------------------------------------------------------------------------
# # Not run
# print(af_fit_4b, pars = c("d", "mu", "delta"))


## ----af_4b_pp_plot-----------------------------------------------------------------------------------------------------------
plot_prior_posterior(af_fit_4b, prior = c("reg", "het"))


## ----af_4b_releff_plot, fig.height = 16, eval=FALSE--------------------------------------------------------------------------
# # Not run
# (af_4b_releff <- relative_effects(af_fit_4b, trt_ref = "Placebo/Standard care"))
# plot(af_4b_releff, ref_line = 0)


## ----af_4b_releff_01_plot----------------------------------------------------------------------------------------------------
(af_4b_releff_01 <- relative_effects(af_fit_4b, 
                                     trt_ref = "Placebo/Standard care",
                                     newdata = data.frame(stroke = c(0, 1), 
                                                          label = c("stroke = 0", "stroke = 1")),
                                     study = label))
plot(af_4b_releff_01, ref_line = 0)


## ----af_4b_betas-------------------------------------------------------------------------------------------------------------
plot(af_fit_4b, pars = "beta", stat = "halfeye", ref_line = 0)


## ----af_4b_betas_transformed-------------------------------------------------------------------------------------------------
af_4b_beta <- as.array(af_fit_4b, pars = "beta")

# Subtract beta[Control:stroke] from the other class interactions
af_4b_beta[ , , 2:3] <- sweep(af_4b_beta[ , , 2:3], 1:2, 
                              af_4b_beta[ , , "beta[.trtclassControl:stroke]"], FUN = "-")

# Set beta[Anti-coagulant:stroke] = -beta[Control:stroke]
af_4b_beta[ , , "beta[.trtclassControl:stroke]"] <- -af_4b_beta[ , , "beta[.trtclassControl:stroke]"]
names(af_4b_beta)[1] <- "beta[.trtclassAnti-coagulant:stroke]"

# Summarise
summary(af_4b_beta)
plot(summary(af_4b_beta), stat = "halfeye", ref_line = 0)


## ----af_4b_ranks-------------------------------------------------------------------------------------------------------------
(af_4b_ranks <- posterior_ranks(af_fit_4b,
                                newdata = data.frame(stroke = c(0, 1), 
                                                     label = c("stroke = 0", "stroke = 1")), 
                                study = label))
plot(af_4b_ranks)

## ----af_4b_rankprobs, fig.height=12------------------------------------------------------------------------------------------
(af_4b_rankprobs <- posterior_rank_probs(af_fit_4b,
                                         newdata = data.frame(stroke = c(0, 1), 
                                                              label = c("stroke = 0", "stroke = 1")), 
                                         study = label))

# Modify the default output with ggplot2 functionality
library(ggplot2)
plot(af_4b_rankprobs) + 
  facet_grid(Treatment~Study, labeller = label_wrap_gen(20)) + 
  theme(strip.text.y = element_text(angle = 0))

## ----af_4b_cumrankprobs, fig.height=12---------------------------------------------------------------------------------------
(af_4b_cumrankprobs <- posterior_rank_probs(af_fit_4b, cumulative = TRUE,
                                            newdata = data.frame(stroke = c(0, 1), 
                                                                 label = c("stroke = 0", "stroke = 1")), 
                                            study = label))

plot(af_4b_cumrankprobs) + 
  facet_grid(Treatment~Study, labeller = label_wrap_gen(20)) + 
  theme(strip.text.y = element_text(angle = 0))


## ----------------------------------------------------------------------------------------------------------------------------
(af_dic_1 <- dic(af_fit_1))

## ----------------------------------------------------------------------------------------------------------------------------
(af_dic_4b <- dic(af_fit_4b))


## ----af_1_resdev_plot--------------------------------------------------------------------------------------------------------
plot(af_dic_1)


## ----af_4b_resdev_plot-------------------------------------------------------------------------------------------------------
plot(af_dic_4b)


## ----atrial_fibrillation_tests, include=FALSE, eval=params$run_tests---------------------------------------------------------
#--- Test against TSD 2 results ---
library(testthat)
library(dplyr)

tol <- 0.05
tol_dic <- 0.1

# No covariates

Cooper_1_releff <- tribble(
~trt                                       , ~est , ~lower, ~upper,
"Low adjusted dose anti-coagulant"         , -1.08,-1.77  , -0.37 ,
"Standard adjusted dose anti-coagulant"    , -0.76,-1.16  , -0.36 ,
"Fixed dose warfarin"                      , 0.18 ,-0.73  , 1.06  ,
"Low dose aspirin"                         , -0.15,-0.56  , 0.27  ,
"Medium dose aspirin"                      , -0.37,-0.83  , 0.07  ,
"High dose aspirin"                        , -0.25,-1.72  , 1.23  ,
"Alternate day aspirin"                    , -1.67,-4.54  , 0.41  ,
"Ximelagatran"                             , -0.84,-1.50  , -0.18 ,
"Triflusal"                                , -0.11,-1.35  , 1.20  ,
"Indobufen"                                , -0.52,-1.47  , 0.47  ,
"Dipyridamole"                             , -0.18,-1.02  , 0.66  ,
"Fixed dose warfarin + low dose aspirin"   , -0.29,-1.09  , 0.51  ,
"Fixed dose warfarin + medium dose aspirin", 0.13 ,-0.60  , 0.83  ,
"Acenocoumarol"                            , -1.56,-3.31  , 0.06  ,
"Low dose aspirin + copidogrel"            , -0.24,-1.06  , 0.57  ,
"Low dose aspirin + dipyridamole"          , -0.49,-1.38  , 0.38  ,
) %>% 
  mutate(trt = ordered(trt, levels = levels(af_net$treatments))) %>%
  arrange(trt)

af_1_releff_df <- as.data.frame(af_1_releff)

test_that("Relative effects (no covariates)", {
  expect_equivalent(af_1_releff_df$mean, Cooper_1_releff$est, tolerance = tol)
  expect_equivalent(af_1_releff_df$`2.5%`, Cooper_1_releff$lower, tolerance = tol)
  expect_equivalent(af_1_releff_df$`97.5%`, Cooper_1_releff$upper, tolerance = tol)
})

af_1_tau <- as.data.frame(summary(af_fit_1, pars = "tau"))

test_that("Heterogeneity SD (no covariates)", {
  expect_equivalent(af_1_tau$`50%`, 0.28, tolerance = tol)
  expect_equivalent(af_1_tau$`2.5%`, 0.02, tolerance = tol)
  expect_equivalent(af_1_tau$`97.5%`, 0.57, tolerance = tol)
})

test_that("DIC (no covariates)", {
  expect_equivalent(af_dic_1$resdev, 60.22, tolerance = tol_dic)
  expect_equivalent(af_dic_1$pd, 48.35, tolerance = tol_dic)
  expect_equivalent(af_dic_1$dic, 108.57, tolerance = tol_dic)
})

test_that("SUCRAs", {
  af_ranks_1 <- posterior_ranks(af_fit_1, sucra = TRUE)
  af_rankprobs_1 <- posterior_rank_probs(af_fit_1, sucra = TRUE)
  af_cumrankprobs_1 <- posterior_rank_probs(af_fit_1, cumulative = TRUE, sucra = TRUE)
  
  expect_equal(af_ranks_1$summary$sucra, af_rankprobs_1$summary$sucra)
  expect_equal(af_ranks_1$summary$sucra, af_cumrankprobs_1$summary$sucra)
})


# Check construction of all contrasts
af_1_releff_all_contr <- relative_effects(af_fit_1, all_contrasts = TRUE)

# Reconstruct from basic contrasts in each study
dk <- function(study, trt, sims) {
  if (trt == "Placebo/Standard care") return(0)
  else if (is.na(study)) {
    return(sims[ , , paste0("d[", trt, "]"), drop = FALSE])
  } else {
    return(sims[ , , paste0("d[", study, ": ", trt, "]"), drop = FALSE])
  }
}

test_af_1_all_contr <- tibble(
  contr = af_1_releff_all_contr$summary$parameter,
  .trtb = factor(stringr::str_extract(contr, "(?<=\\[)(.+)(?= vs\\.)"), levels = levels(af_net$treatments)),
  .trta = factor(stringr::str_extract(contr, "(?<=vs\\. )(.+)(?=\\])"), levels = levels(af_net$treatments))
) %>%
  rowwise() %>%
  mutate(as_tibble(multinma:::summary.mcmc_array(dk(NA, .trtb, af_1_releff$sims) - dk(NA, .trta, af_1_releff$sims)))) %>%
  select(.trtb, .trta, parameter = contr, mean:Rhat)

test_that("Construction of all contrasts is correct (no covariates)", {
  ntrt <- nlevels(af_net$treatments)
  expect_equal(nrow(af_1_releff_all_contr$summary), ntrt * (ntrt - 1) / 2)
  expect_equal(select(af_1_releff_all_contr$summary, -Rhat),
               select(test_af_1_all_contr, -Rhat),
               check.attributes = FALSE)
})

# With stroke covariate, shared interactions

Cooper_4b_releff <- tribble(
~trt                                       , ~est  , ~lower, ~upper, 
"Low adjusted dose anti-coagulant"         , -1.20 ,-1.89  , -0.54 ,
"Standard adjusted dose anti-coagulant"    , -0.77 ,-1.14  , -0.38 ,
"Fixed dose warfarin"                      , -0.11 ,-0.90  , 0.72  ,
"Low dose aspirin"                         , -0.08 ,-0.47  , 0.30  ,
"Medium dose aspirin"                      , -0.45 ,-0.87  , -0.03 ,
"High dose aspirin"                        , -0.39 ,-1.86  , 1.11  ,
"Alternate day aspirin"                    , -1.74 ,-5.16  , 0.48  ,
"Ximelagatran"                             , -0.86 ,-1.42  , -0.27 ,
"Triflusal"                                , 0.13  ,-1.05  , 1.38  ,
"Indobufen"                                , -1.21 ,-2.26  , -0.13 ,
"Dipyridamole"                             , -0.21 ,-1.01  , 0.58  ,
"Fixed dose warfarin + low dose aspirin"   , 0.54  ,-0.80  , 1.85  ,
"Fixed dose warfarin + medium dose aspirin", 0.12  ,-0.53  , 0.80  ,
"Acenocoumarol"                            , -0.534,-2.67  , 1.38  ,
"Low dose aspirin + copidogrel"            , -0.14 ,-0.82  , 0.53  ,
"Low dose aspirin + dipyridamole"          , -0.53 ,-1.38  , 0.30  ,
) %>% 
  mutate(trt = ordered(trt, levels = levels(af_net$treatments))) %>%
  arrange(trt)

af_4b_releff_Cooper <- as.data.frame(relative_effects(af_fit_4b, 
                                                      newdata = tibble(stroke = 0.27),
                                                      trt_ref = "Placebo/Standard care"))

test_that("Relative effects (common interaction)", {
  expect_equivalent(af_4b_releff_Cooper$mean, Cooper_4b_releff$est, tolerance = tol)
  expect_equivalent(af_4b_releff_Cooper$`2.5%`, Cooper_4b_releff$lower, tolerance = tol)
  expect_equivalent(af_4b_releff_Cooper$`97.5%`, Cooper_4b_releff$upper, tolerance = tol)
})

Cooper_4b_beta <- tribble(
~trt_class       , ~est , ~lower, ~upper,
"Anti-coagulant", -0.71,-1.58  , 0.15  ,
"Anti-platelet" , 0.23 ,-0.45  , 0.93  ,
#"Mixed"          , 3.05 ,-1.26  , 7.30  ,
"Mixed"          , 3.21 ,-0.91  , 7.30  ,
)

af_4b_beta_df <- as.data.frame(summary(af_4b_beta))

test_that("Interaction estimates (common interaction)", {
  expect_equivalent(af_4b_beta_df$mean, Cooper_4b_beta$est, tolerance = tol)
  skip_on_ci()
  expect_equivalent(af_4b_beta_df$`2.5%`, Cooper_4b_beta$lower, tolerance = tol)
  expect_equivalent(af_4b_beta_df$`97.5%`, Cooper_4b_beta$upper, tolerance = tol)
})

af_4b_tau <- as.data.frame(summary(af_fit_4b, pars = "tau"))

test_that("Heterogeneity SD (common interaction)", {
  expect_equivalent(af_4b_tau$`50%`, 0.19, tolerance = tol)
  expect_equivalent(af_4b_tau$`2.5%`, 0.01, tolerance = tol)
  expect_equivalent(af_4b_tau$`97.5%`, 0.48, tolerance = tol)
})

test_that("DIC (common interaction)", {
  expect_equivalent(af_dic_4b$resdev, 58.74, tolerance = tol_dic)
  expect_equivalent(af_dic_4b$pd, 48.25, tolerance = tol_dic)
  expect_equivalent(af_dic_4b$dic, 106.99, tolerance = tol_dic)
})

test_that("SUCRAs", {
  stroke_01 <- data.frame(stroke = c(0, 1), label = c("stroke = 0", "stroke = 1"))
  af_ranks_4b <- posterior_ranks(af_fit_4b, newdata = stroke_01, 
                                study = label, sucra = TRUE)
  af_rankprobs_4b <- posterior_rank_probs(af_fit_4b, newdata = stroke_01, 
                                           study = label, sucra = TRUE)
  af_cumrankprobs_4b <- posterior_rank_probs(af_fit_4b, cumulative = TRUE, newdata = stroke_01,
                                              study = label, sucra = TRUE)
  
  expect_equal(af_ranks_4b$summary$sucra, af_rankprobs_4b$summary$sucra)
  expect_equal(af_ranks_4b$summary$sucra, af_cumrankprobs_4b$summary$sucra)
})

# Check construction of all contrasts
af_4b_releff <- relative_effects(af_fit_4b, trt_ref = "Placebo/Standard care")
af_4b_releff_all_contr <- relative_effects(af_fit_4b, all_contrasts = TRUE)

test_af_4b_all_contr <- tibble(
  contr = af_4b_releff_all_contr$summary$parameter,
  .study = factor(stringr::str_extract(contr, "(?<=\\[)(.+)(?=:)")),
  .trtb = factor(stringr::str_extract(contr, "(?<=\\: )(.+)(?= vs\\.)"), levels = levels(af_net$treatments)),
  .trta = factor(stringr::str_extract(contr, "(?<=vs\\. )(.+)(?=\\])"), levels = levels(af_net$treatments))
) %>% 
  rowwise() %>% 
  mutate(as_tibble(multinma:::summary.mcmc_array(dk(.study, .trtb, af_4b_releff$sims) - dk(.study, .trta, af_4b_releff$sims)))) %>% 
  select(.study, .trtb, .trta, parameter = contr, mean:Rhat)

test_that("Construction of all contrasts is correct (common interaction)", {
  ntrt <- nlevels(af_net$treatments)
  nstudy <- nlevels(test_af_4b_all_contr$.study)
  expect_equal(nrow(af_4b_releff_all_contr$summary), nstudy * ntrt * (ntrt - 1) / 2)
  expect_equal(select(af_4b_releff_all_contr$summary, -Rhat), 
               select(test_af_4b_all_contr, -Rhat), 
               check.attributes = FALSE)
})

# Check construction of all contrasts in target population
af_4b_releff_new <- relative_effects(af_fit_4b, newdata = tibble(stroke = 0.27), trt_ref = "Placebo/Standard care")
af_4b_releff_all_contr_new <- relative_effects(af_fit_4b, newdata = tibble(stroke = 0.27), all_contrasts = TRUE)

test_af_4b_all_contr_new <- tibble(
  contr = af_4b_releff_all_contr_new$summary$parameter,
  .study = factor(stringr::str_extract(contr, "(?<=\\[)(.+)(?=:)")),
  .trtb = factor(stringr::str_extract(contr, "(?<=\\: )(.+)(?= vs\\.)"), levels = levels(af_net$treatments)),
  .trta = factor(stringr::str_extract(contr, "(?<=vs\\. )(.+)(?=\\])"), levels = levels(af_net$treatments))
) %>% 
  rowwise() %>% 
  mutate(as_tibble(multinma:::summary.mcmc_array(dk(.study, .trtb, af_4b_releff_new$sims) - dk(.study, .trta, af_4b_releff_new$sims)))) %>% 
  select(.study, .trtb, .trta, parameter = contr, mean:Rhat)

test_that("Construction of all contrasts in target population is correct (common interaction)", {
  ntrt <- nlevels(af_net$treatments)
  nstudy <- nlevels(test_af_4b_all_contr_new$.study)
  expect_equal(nrow(af_4b_releff_all_contr_new$summary), nstudy * ntrt * (ntrt - 1) / 2)
  expect_equal(select(af_4b_releff_all_contr_new$summary, -Rhat), 
               select(test_af_4b_all_contr_new, -Rhat), 
               check.attributes = FALSE)
})

test_that("Robust to custom options(contrasts) settings", {
  af_fit_4b_SAS <- withr::with_options(list(contrasts = c(ordered = "contr.SAS",
                                                       unordered = "contr.SAS")),
             nowarn_on_ci(nma(af_net, 
                 seed = 579212814,
                 trt_effects = "random",
                 regression = ~ .trt:stroke,
                 class_interactions = "common",
                 QR = TRUE,
                 prior_intercept = normal(scale = 100),
                 prior_trt = normal(scale = 100),
                 prior_reg = normal(scale = 100),
                 prior_het = half_normal(scale = 5),
                 adapt_delta = 0.99,
                 iter = 5000)))

  expect_equal(as_tibble(summary(af_fit_4b_SAS))[, c("parameter", "mean", "sd")],
               as_tibble(summary(af_fit_4b))[, c("parameter", "mean", "sd")],
               tolerance = tol)
  expect_equal(as_tibble(relative_effects(af_fit_4b_SAS))[, c("parameter", "mean", "sd")],
               as_tibble(relative_effects(af_fit_4b))[, c("parameter", "mean", "sd")],
               tolerance = tol)
})



# Force clean up
rm(list = ls())
gc()

