

do_parallel_start = function(df, cores, partition_by){

	cl <- multidplyr::new_cluster(cores)

	df %>%
		dplyr::left_join(
			(.) %>%
				dplyr::select(!!partition_by) %>%
				dplyr::distinct() %>%
				dplyr::mutate(
					part = 1:n() %>%
						magrittr::divide_by(length((.))) %>%
						magrittr::multiply_by(!!cores) %>%
						ceiling
				)
		)  %>%
		group_by(part) %>%
		multidplyr::partition(cl)
}


do_parallel_end = function(.){
	(.) %>%
		dplyr::collect() %>%
		dplyr::ungroup() %>%
		dplyr::select(-part)
}
