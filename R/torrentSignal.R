
#' @export
bamGetHeader <- function(bam_file) {
    if(!file.exists(bam_file)) return(NULL)
    header_info <- scanBamHeader(bam_file)
    return(header_info)
}

#' @export
getVariantSignal <- function(
    bam_file,
    variant = NULL,
    # expect in the format: "chr20-31022441-A-AG"
    # where A is the reference sequence and AG is the alternate
    # e.g. ASXL1 c.1934del = "chr20-31022441-A-AG"
    genome = "hg19" # choices are hg19 or hg38
) {
    # sanity check
    if(!file.exists(bam_file)) stop(bam_file, " not found")
    if(!(genome %in% c("hg19", "hg38"))) stop("genome must be either hg19 or hg38")

    var_elems <- tstrsplit(variant, split = "-")
    if(length(var_elems) != 4) stop(
        "variant must be of the format chrX-123456-G-GA, ",
        "where 123456 is the genomic locus, G is the reference sequence ",
        "and GA is the variant sequence"
    )
    
    # check variant is valid
    CHR <- var_elems[[1]]
    POS <- as.numeric(var_elems[[2]])
    REF <- var_elems[[3]]
    ALT <- var_elems[[4]]

    header <- bamGetHeader(bam_file)
    if(!( CHR %in% names(header[[1]]$targets) )) stop(
        "chromosome not found in BAM file - check the variant name is correct"
    )
    idx_chr <- which( names(header[[1]]$targets) == CHR )
    chr_len <- header[[1]]$targets[idx_chr]
    if(!( POS > 0 & POS < chr_len)) stop(
        "variant position is invalid - check variant name is correct"
    )
    
    len_ref <- nchar(REF)
    gr_ref <- GRanges(
        CHR,
        IRanges(POS, POS + len_ref - 1)
    )
    if(genome == "hg19") {
        genome_obj <- BSgenome.Hsapiens.UCSC.hg19
    } else {
        genome_obj <- BSgenome.Hsapiens.UCSC.hg38   
    }
    seq_ref <- as.character(getSeq(genome_obj, gr_ref))
    if(seq_ref != REF) {
        "reference sequence incorrect - check genome and variant is correct"
    }
    
    # Identify coordinates of the bin of interest
    # consider the sequence as bins of unique nucleotides
    # i.e. homopolymers are considered a single bin
    
    if(nchar(REF) == 1 & nchar(ALT) == 1) stop(
        "Only indels are supported at this stage"
    )
    if(nchar(REF) > 1 & nchar(ALT) > 1) stop(
        "Only indels are supported at this stage"
    )
    if(nchar(REF) > 2) stop(
        "Only 1-nt inserts and deletions are supported at this stage"
    )
    if(nchar(ALT) > 2) stop(
        "Only 1-nt inserts and deletions are supported at this stage"
    )
    
    # determine length of homopolymer
    gr_hp <- GRanges(
        CHR,
        IRanges(POS + 1, POS + 100)
    )
    seq_hp <- as.character(getSeq(genome_obj, gr_hp))
    matches <- gregexpr("(.)\\1*", seq_hp)
    blocks <- regmatches(seq_hp, matches)[[1]]
    
    homoLen <- nchar(blocks[1])
    gr_bin <- GRanges(
        CHR,
        IRanges(POS + 1, POS + homoLen)
    )
    seq_bin <- as.character(getSeq(genome_obj, gr_bin))
    
    binned_region <- paste0(CHR, ":", 
        as.character(POS+1), "-", as.character(POS+homoLen) )
    
    header_text_lines <- header[[1]]$text
    header_text_lines <- header_text_lines[names(header_text_lines) == "@RG"]    
    
    df_RG <- NULL
    for(i in seq_len(length(header_text_lines))) {
        line <- header_text_lines[[i]]
        df_line <- data.frame(
            read_group = "",
            key_sequence = "",
            flow_order = ""
        )
        for(j in seq_len(length(line))) {
            if(substr(line[j], 1, 3) == "ID:") {
                df_line$read_group <- substr(line[j], 4, nchar(line[j]))
            }
            if(substr(line[j], 1, 3) == "KS:") {
                df_line$key_sequence <- substr(line[j], 4, nchar(line[j]))
            }
            if(substr(line[j], 1, 3) == "FO:") {
                df_line$flow_order <- substr(line[j], 4, nchar(line[j]))
            }
        }
        if(is.null(df_RG)){
            df_RG <- df_line
        } else {
            df_RG <- rbind(df_RG, df_line)
        }
    }
    
    return(
        c_oneSignal(
            bam_file, binned_region, df_RG
        )
    )
}