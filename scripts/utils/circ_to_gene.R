#' Function to create BED files
#' 
#' @param df data.frame. Expression matrix
#' 
IDstoBED <- function(df) {

    coords <- str_split(circIDs, "\\.", simplify = T) |> as.data.frame() 
    colnames(coords) <- c("chr", "start", "end")
    # catch odd cases
    coords <- coords[-grep("e", coords$start),]
    coords <- coords[as.numeric(coords$end) >= as.numeric(coords$start) - 1, ] |>
        suppressWarnings()
    coords <-  coords[complete.cases(coords),]
    # make GRanges
    bed <- makeGRangesFromDataFrame(coords)
    mcols(bed)$circID <- paste(coords$chr, coords$start, coords$end, sep = ".")
    
    return(bed)
}

#' Function to map circRNA to genes
#' 
#' @param gr GRanges object. Output of IDStoBED()
#' 
circToGene <- function(gr) {

    print("---Starting circToGene mapping")

    circ_genes <- list()

    for (i in 1:length(gr)) {
        genes <- c()
        circRNA = gr[i]

        # get exons in circRNA
        overlaps <- findOverlaps(unlist(exons), circRNA, type = "within")
        exons_within <- unlist(exons)[queryHits(overlaps)]
        n_exons <- length(exons_within)

        # skip if there are no exons
        if (n_exons == 0) {
            circ_genes[[circRNA$circID]] <- character(0)
            next
        }

        # table of number of exons per gene
        gene_counts <- table(names(exons_within))

        # keep genes where number of overlapping exons >= half number of exons in circRNA
        genes <- names(gene_counts[gene_counts >= floor(n_exons / 2)])
        circ_genes[[circRNA$circID]] <- unname(mcols(exons[genes])$gene_name)
    }

    # creating mapping df of circ to gene from circ_genes list
    mapping <- do.call(rbind, lapply(names(circ_genes), function(col) {
        if (length(circ_genes[[col]]) > 0) {
            data.frame(old_name = col, gene = circ_genes[[col]])
        } else {
            NULL
        }
    })) |> as.data.frame(stringsAsFactors = FALSE)

    return(mapping)
}