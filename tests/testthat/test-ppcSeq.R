context('ppcSeq')

test_that("Quick test",{

  FDR_threshold = 0.01

  res =
    ppc_seq(
      dplyr::mutate(ppcSeq::counts, is_significant = FDR < FDR_threshold),
      formula = ~ Label,
      significance_column = PValue,
      do_check_column  = is_significant,
      value_column = value,
      percent_false_positive_genes = "5%", tol_rel_obj = 0.01
    )

  expect_equal(

    as.integer(unlist(res[c(1, 2, 4, 5, 6),5])),
    c(0,1,0,0,0)
  )

})
