context('ppcSeq')

test_that("VB post approx",{

  res =
    ppcSeq::ppc_seq(
      dplyr::mutate(ppcSeq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
      percent_false_positive_genes = "1%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = T,
      approximate_posterior_analysis = T,
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
    ppcSeq::ppc_seq(
      dplyr::mutate(ppcSeq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
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
    ppcSeq::ppc_seq(
      dplyr::mutate(ppcSeq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
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
    ppcSeq::ppc_seq(
      dplyr::mutate(ppcSeq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
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
