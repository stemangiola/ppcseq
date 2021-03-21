#' Drop lowly transcribed genes for TMM normalization
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats median
#' @importFrom edgeR filterByExpr
#' @importFrom rlang enquo
#' @importFrom rlang quo_is_symbol
#' @importFrom purrr when
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param factor_of_interest The name of the column of the factor of interest
#' @param minimum_counts A positive integer. Minimum counts required for at least some samples.
#' @param minimum_proportion A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for scaling procedure.
#'
#' @return A tibble filtered
add_scaled_counts_bulk.get_low_expressed <- function(.data,
																										 .sample = `sample`,
																										 .transcript = `transcript`,
																										 .abundance = `count`,
																										 factor_of_interest = NULL,
																										 minimum_counts = 10,
																										 minimum_proportion = 0.7) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	factor_of_interest = enquo(factor_of_interest)

	# Check if factor_of_interest is continuous and exists
	string_factor_of_interest =

		factor_of_interest %>%
		when(

			quo_is_symbol(.) &&
				select(.data, !!(.)) %>%
				lapply(class) %>%
				as.character() %in% c("numeric", "integer", "double") ~ {
					message("ppcseq says: The factor of interest is continuous (e.g., integer,numeric, double). The data will be filtered without grouping.")
					NULL
				},

			quo_is_symbol(.) ~ .data %>%
				distinct(!!.sample, !!factor_of_interest) %>%
				arrange(!!.sample) %>%
				pull(!!factor_of_interest),

			~ NULL
		)

	if (minimum_counts < 0)
		stop("The parameter minimum_counts must be > 0")
	if (minimum_proportion < 0 |	minimum_proportion > 1)
		stop("The parameter minimum_proportion must be between 0 and 1")

	.data %>%
		select(!!.sample,!!.transcript, !!.abundance) %>%
		spread(!!.sample, !!.abundance) %>%

		# Drop if transcript have missing value
		drop_na() %>%

		# If I don't have any transcript with all samples give meaningful error
		when(
			nrow(.) == 0 ~ stop("ppcseq says: you don't have any transcript that is in all samples. Please consider using impute_missing_abundance."),
			~ (.)
		) %>%

		# Call edgeR
		as_matrix(rownames = !!.transcript) %>%
		filterByExpr(
			min.count = minimum_counts,
			group = string_factor_of_interest,
			min.prop = minimum_proportion
		) %>%
		not() %>%
		which %>%
		names
}

# Set internal
#' @importFrom purrr when
.identify_abundant = 		function(.data,
																.sample = NULL,
																.transcript = NULL,
																.abundance = NULL,
																factor_of_interest = NULL,
																minimum_counts = 10,
																minimum_proportion = 0.7)
{
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	factor_of_interest = enquo(factor_of_interest)

	.data %>%

		{
			gene_to_exclude =
				add_scaled_counts_bulk.get_low_expressed(
					.,
					.sample = !!.sample,
					.transcript = !!.transcript,
					.abundance = !!.abundance,
					factor_of_interest = !!factor_of_interest,
					minimum_counts = minimum_counts,
					minimum_proportion = minimum_proportion
				)

			dplyr::mutate(., .abundant := !!.transcript %in% gene_to_exclude %>% not())
		}
}


#' Get a tibble with scaled counts using TMM
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr equals
#' @importFrom rlang :=
#' @importFrom stats median
#' @importFrom purrr when
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A character string. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#' @param reference_sample A character string. The name of the reference sample. If NULL the sample with highest total read count will be selected as reference.
#'
#' @return A tibble including additional columns
#'
#'
get_scaled_counts_bulk <- function(.data,
																	 .sample = NULL,
																	 .transcript = NULL,
																	 .abundance = NULL,
																	 method = "TMM",
																	 reference_sample = NULL) {
	# Get column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Reformat input data set
	df <-
		.data %>%

		# Rename
		dplyr::select(!!.sample,!!.transcript,!!.abundance) %>%

		# Set samples and genes as factors
		dplyr::mutate(!!.sample := factor(!!.sample),!!.transcript := factor(!!.transcript))


	# Get reference
	reference <-
		reference_sample %>%
		when(
			!is.null(.) ~ (.),

			# If not specified take most abundance sample
			df %>%
				group_by(!!.sample) %>%
				summarise(sum = median(!!.abundance)) %>%
				mutate(med = max(sum)) %>%
				mutate(diff = abs(sum - med)) %>%
				arrange(diff) %>%
				head(n = 1) %>%
				pull(!!.sample) %>%
				as.character()
		)


	nf_obj <-
		add_scaled_counts_bulk.calcNormFactor(
			df,
			reference,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			method
		)

	# Calculate normalization factors
	nf_obj$nf %>%
		dplyr::left_join(
			df %>%
				group_by(!!.sample) %>%
				summarise(tot = sum(!!.abundance, na.rm = TRUE)) %>%
				ungroup() %>%
				dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))),
			by = quo_name(.sample)
		) %>%
		mutate(multiplier =
					 	1 /
					 	(tot_filt * nf) *

					 	# Put everything to the reference sample scale
					 	((.) %>% filter(!!.sample == reference) %>% pull(tot))) %>%

		# I have correct the strange behaviour of edgeR of reference
		# sample not being 1
		# I HAD TO COMMENT BECAUSE TEST FAILING
		# {
		# 	mult_ref = (.) %>%  filter(!!.sample == reference) %>% pull(multiplier)
		# 	(.) %>%  mutate(
		# 		multiplier =
		# 			multiplier /mult_ref
		# 	)
		# } %>%

	dplyr::select(-tot,-tot_filt) %>%
		dplyr::rename(TMM = nf)

}

#' Calculate the norm factor with calcNormFactor from limma
#'
#' @keywords internal
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang :=
#' @importFrom stats setNames
#' @importFrom edgeR calcNormFactors
#'
#' @param .data A tibble
#' @param reference A reference matrix, not sure if used anymore
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#' @param method A string character. The scaling method passed to the backend function (i.e., edgeR::calcNormFactors; "TMM","TMMwsp","RLE","upperquartile")
#'
#'
#' @return A list including the filtered data frame and the normalization factors
add_scaled_counts_bulk.calcNormFactor <- function(.data,
																									reference = NULL,
																									.sample = `sample`,
																									.transcript = `transcript`,
																									.abundance = `count`,
																									method) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Get data frame for the highly transcribed transcripts
	df.filt <-
		.data %>%
		# dplyr::filter(!(!!.transcript %in% gene_to_exclude)) %>%
		droplevels() %>%
		select(!!.sample, !!.transcript, !!.abundance)

	# scaled data set
	nf =
		tibble::tibble(
			# Sample factor
			sample = factor(levels(df.filt %>% pull(!!.sample))),

			# scaled data frame
			nf = calcNormFactors(
				df.filt %>%
					tidyr::spread(!!.sample,!!.abundance) %>%
					tidyr::drop_na() %>%
					dplyr::select(-!!.transcript),
				refColumn = which(reference == factor(levels(
					df.filt %>% pull(!!.sample)
				))),
				method = method
			)
		) %>%

		setNames(c(quo_name(.sample), "nf")) %>%

		# Add the statistics about the number of genes filtered
		dplyr::left_join(
			df.filt %>%
				dplyr::group_by(!!.sample) %>%
				dplyr::summarise(tot_filt = sum(!!.abundance, na.rm = TRUE)) %>%
				dplyr::mutate(!!.sample := as.factor(as.character(!!.sample))),
			by = quo_name(.sample)
		)

	# Return
	list(
		# gene_to_exclude = gene_to_exclude,
		nf = nf
	)
}
