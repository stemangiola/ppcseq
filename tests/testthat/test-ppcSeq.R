context('ppcseq')

data("counts")

if(Sys.info()[['sysname']] != "Windows"){

test_that("VB post approx no correction",{

  res =
    counts %>%
    dplyr::mutate(  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), TRUE, FALSE) ) %>%
    ppcseq::identify_outliers(
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      percent_false_positive_genes = 1,
      tol_rel_obj = 0.01,
      approximate_posterior_inference =TRUE,
      approximate_posterior_analysis =TRUE,
      how_many_negative_controls = 50,
      cores=1,pass_fit = TRUE
    )

  expect_equal(

    as.integer(unlist(res[,4])),
    c(0,1,0)
  )

})

test_that("VB post full",{

  res =
    ppcseq::identify_outliers(
      dplyr::mutate(counts,  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), TRUE, FALSE) ),
      formula = ~ Label,
      sample, symbol, value,
      .significance = PValue,
      .do_check  = is_significant,
      percent_false_positive_genes = 1,
      tol_rel_obj = 0.01,
      approximate_posterior_inference =TRUE,
      approximate_posterior_analysis = FALSE,
      how_many_negative_controls = 50,
      cores=1
    )

  expect_equal(

    as.integer(unlist(res[,4])),
    c(0,1,0)
  )

})

} else{
  print("tests for windows are temporarily disabled")
}
