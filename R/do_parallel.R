

do_parallel_start = function(df, cores, partition_by){

	# Only if cores > 1
	if(cores > 1)		cl <- multidplyr::new_cluster(cores)

	df %>%
		dplyr::left_join(
			(.) %>%
				dplyr::select(!!partition_by) %>%
				dplyr::distinct() %>%
				dplyr::mutate(
					`.part` = 1:n() %>%
						magrittr::divide_by(length((.))) %>%
						magrittr::multiply_by(!!cores) %>%
						ceiling
				)
		)  %>%
		group_by(.part) %>%

		# Only if cores > 1
		ifelse_pipe(cores > 1,	~ .x %>% multidplyr::partition(cl))

}


do_parallel_end = function(.){
	(.) %>%
		# Only if cores > 1
		ifelse_pipe((.) %>% class %>% magrittr::equals("multidplyr_party_df") %>% any,	~ .x %>% dplyr::collect()) %>%
		dplyr::ungroup() %>%

		# Only if cores > 1
		dplyr::select(-`.part`)
}
