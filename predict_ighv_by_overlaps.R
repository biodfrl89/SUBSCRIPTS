# Name: predict_ighv_by_overlaps.R
# Description: This script is designed to predict IGHV genes and their respective exons, using the results from exonerate,
# reducing the sequences and finding overlaps between gene and exon anotations.
# Author: David Felipe Rendón Luna 
# Date: September-2023


# CHECK LIBRARIES AND ARGUMENTS -------------------------------------------

# Check for optparse library to load arguments from command line
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# Check for aditional libraries
if(nzchar(system.file(package = "rtracklayer"))) {
  cat("-rtracklayer library found.\n")
} else {
  stop("rtracklayer was not found. Exiting.\n")
}

if(nzchar(system.file(package = "GenomicRanges"))) {
  cat("-GenomicRanges library found.\n")
} else {
  stop("GenomicRanges was not found. Exiting.\n")
}

# Load parser
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-q", "--query"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Exonerate filtered file ", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-s", "--subject"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "", 
                            metavar = "[FILENAME]"))



# PREPARE QUERY AND SUBJECT -----------------------------------------------

# Import BED file of manual anotation
print("Reading gene annotation by exonerate.")
exonerate_gene_anot <- rtracklayer::import("./DATA/Arja_exonerate_gene_protein.gff", format = "GFF")

# Make ranges reduction
reduced_exonerate_gene_anot <- GenomicRanges::reduce(exonerate_gene_anot)

# Edit names and get the element code
#gene_names <- sapply(exonerate_gene_anot@elementMetadata@listData$sequence, function(x) gsub('\\|.*', '', x ), USE.NAMES = FALSE)

# Substitute each element with only the element code
#manual_anotation@elementMetadata@listData$name <- gene_names 

# Import reduced GFF file
print("Reading exon annotation by exonerate.")
exonerate_exon_anot <- rtracklayer::import("./DATA/Arja_exonerate_exon_protein.gff", format = "GFF")

# Make ranges reduction
reduced_exonerate_exon_anot <- GenomicRanges::reduce(exonerate_exon_anot)


# MAKE OVERLAPS AND FIND INDEX OF TWO OR MORE -----------------------------

# Make overlaps
overlaps <- GenomicRanges::findOverlaps(query = reduced_exonerate_gene_anot, subject = reduced_exonerate_exon_anot )

# Produce a table of number of overlaps in each query
overlaps_table <- table(overlaps@from)

# Get the exonerate "genes" that have two or more overlaps with exonerate "exons"
overlaps_two_more <- as.numeric(names(overlaps_table[overlaps_table >= 2]))

# Use names to filter the desired overlaps
overlaps[overlaps_two_more]

# Convert Granges to Dataframe
df_overlaps <- as.data.frame(overlaps)

# Find the records of the overlap dataframe which subject having two or more overlaps with query
df_overlaps_two <- df_overlaps[df_overlaps$queryHits %in% overlaps_two_more,]


# ADD ID IN THE QUERY OBJECT (GENE) ---------------------------------------

g = 1
p = 1
ighv_genes_detected <- c()
for (i in 1:length(reduced_exonerate_gene_anot)) {
  if (i %in% unique(df_overlaps_two$queryHits)) {
    ighv_genes_detected <- c(ighv_genes_detected, paste0("Arja_IGHV_gene_", g))
    g = g + 1
  } else {
    ighv_genes_detected <- c(ighv_genes_detected, paste0("Arja_IGHV_pseudogene_", p))
    p = p + 1
  }
}

# Add the created names to the reduced gene anotation object
reduced_exonerate_gene_anot$ID <- ighv_genes_detected


# ADD ID IN THE SUBJECT OBJECT (EXON) -------------------------------------

# Begin empty ID column of metadata
reduced_exonerate_exon_anot$ID <- NA

# Rename the numbers in increasing way
renamed_exons_numbers <- dplyr::dense_rank(df_overlaps_two$queryHits)

# Subset the subjectHits
temp_df <- df_overlaps_two["subjectHits"]

# Paste
temp_df$renamed <- paste0("Arja_IGHV_exon_", renamed_exons_numbers)

# Add the corresponding name of exon following the index of subjectHits
cont = 1
for (i in temp_df$subjectHits) {
  reduced_exonerate_exon_anot[i]$ID <- temp_df$renamed[cont]
  cont = cont + 1
}

# Any remaining with NA is renamed to Sequence
reduced_exonerate_exon_anot[is.na(reduced_exonerate_exon_anot$ID)]$ID <- "Sequence"
