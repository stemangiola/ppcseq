context('ppcSeq')

test_that("Quick test",{

  res =
    ppcSeq::ppc_seq(
      dplyr::mutate(ppcSeq::counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), T, F) ),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
      percent_false_positive_genes = "5%",
      tol_rel_obj = 0.01,
      approximate_posterior_inference = T,
      approximate_posterior_analysis = T
    )

  expect_equal(

    as.integer(unlist(res[c(1, 2, 4, 5, 6),5])),
    c(0,1,0,0,0)
  )

})
