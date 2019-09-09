context('ppcSeq')

test_that("dummy",expect_equal(1,1))

#
#
# test_that("Quick test",{
#
  # FDR_threshold = 0.01
  #
  # res =
  #   ppcSeq::counts %>%
  #   mutate(is_significant = FDR < FDR_threshold) %>%
  #   ppc_seq(
  #     formula = ~ Label,
  #     significance_column = PValue,
  #     do_check_column  = is_significant,
  #     value_column = value,
  #     percent_false_positive_genes = "5%", tol_rel_obj = 0.01
  #   )
#
#   expect_equal(
#
#     length(
#       attr(
#         create_ttBulk(
#           ttBulk::counts_mini,
#           sample_column = sample,
#           transcript_column = transcript,
#           counts_column = `read count`
#         ) ,
#         "parameters"
#       )
#     ),
#     3
#   )
#
# })
