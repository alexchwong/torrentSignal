
#' @export
idxstats <- function(bam_file, n_threads) {
    idxstats_pbam(bam_file, n_threads)
}

c_oneSignal <- function(bam_file, binned_region, rg_df) {
    binned_elem <- tstrsplit(binned_region, split = ":", fixed=T)
    CHR <- binned_elem[[1]]
    start_end <- tstrsplit(binned_elem[[2]], split = "-", fixed=T)
    start <- as.numeric(start_end[[1]])
    end <- as.numeric(start_end[[2]])
    
    c_getSignal(
        bam_file, binned_region,
        CHR, start-1, end,
        rg_df$read_group, rg_df$key_sequence, rg_df$flow_order
    )
}
