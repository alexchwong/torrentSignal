#' @useDynLib torrentSignal, .registration = TRUE
#' @import Rcpp
#' @import data.table
#' @importFrom Rsamtools scanBamHeader
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BSgenome getSeq
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
NULL