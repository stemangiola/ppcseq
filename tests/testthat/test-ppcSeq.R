context('ppcseq')

test_that("VB post approx no correction",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(ppcseq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      .abundance = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = T,
      approximate_posterior_analysis = T,
      do_correct_approx = F,
      how_many_negative_controls = 50,
      cores=1,pass_fit = T
    )

  expect_equal(

    as.integer(unlist(res[,5])),
    c(0,1,0)
  )

})

test_that("VB post approx yes correction",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(ppcseq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      .abundance = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = T,
      approximate_posterior_analysis = T,
      do_correct_approx = T,
      how_many_negative_controls = 50,
      cores=1,pass_fit = T
    )

  expect_equal(

    as.integer(unlist(res[,5])),
    c(0,1,0)
  )

})

test_that("VB post full",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(ppcseq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      .abundance = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = T,
      approximate_posterior_analysis = F,
      how_many_negative_controls = 50,
      cores=1
    )

  expect_equal(

    as.integer(unlist(res[,5])),
    c(0,1,0)
  )

})

test_that("bayes post approx",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(ppcseq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      .abundance = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = F,
      approximate_posterior_analysis = T,
      how_many_negative_controls = 50,
      cores=1
    )

  expect_equal(

    as.integer(unlist(res[,5])),
    c(0,1,0)
  )

})

test_that("bayes full",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(ppcseq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      .abundance = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = F,
      approximate_posterior_analysis = F,
      how_many_negative_controls = 50,
      cores=1
    )

  expect_equal(

    as.integer(unlist(res[,5])),
    c(0,1,0)
  )

})
