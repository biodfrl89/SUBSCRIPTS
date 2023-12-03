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
                            help = "Exonerate filtered file with only exons", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-s", "--subject"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Exonerate filtered file with only genes", 
                            metavar = "[FILENAME]"))

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# Check for arguments -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$query)) {
  print_help(opt_parser)
  stop("A query file must be submitted", call.=FALSE)
}

if (is.null(opt$subject)) {
  print_help(opt_parser)
  stop("A subject file must be submitted", call.=FALSE)
}

# REASIGMENT -----------------------------------------------
# Gene
subject <- opt$subject
# Exons
query <- opt$query

# PREPARE QUERY AND SUBJECT -----------------------------------------------

# Import genes
red_gene_anot_subj <- rtracklayer::import(subject, format = "GFF")

# Import exons
red_exon_anot_query <- rtracklayer::import(query, format = "GFF")

# MAKE OVERLAPS AND FIND INDEX OF TWO OR MORE -----------------------------

# Make overlaps
overlaps <- GenomicRanges::findOverlaps(query = red_exon_anot_query, subject = red_gene_anot_subj )

# Produce a table of number of overlaps in each query
overlaps_table <- table(overlaps@to)

# Get the exonerate "genes" that have two or more overlaps with exonerate "exons"
overlaps_two_more <- as.numeric(names(overlaps_table[overlaps_table >= 2]))

# Use names to filter the desired overlaps
overlaps_filtered <- overlaps[which(overlaps@to %in% overlaps_two_more)]

# Convert Granges to Dataframe
df_overlaps <- as.data.frame(overlaps_filtered)

# ADD ID IN THE QUERY OBJECT (GENE) ---------------------------------------

g = 1
p = 1
ighv_genes_detected <- c()
for (i in seq_along(red_gene_anot_subj )) {
  if (i %in% unique(df_overlaps$queryHits)) {
    ighv_genes_detected <- c(ighv_genes_detected, paste0("IGHV_gene_", g))
    g = g + 1
  } else {
    ighv_genes_detected <- c(ighv_genes_detected, paste0("IGHV_pseudogene_", p))
    p = p + 1
  }
}

# Add the created names to the reduced gene anotation object
red_gene_anot_subj$ID <- ighv_genes_detected

# Save as gff3
rtracklayer::export.gff3(red_gene_anot_subj, "gene_prediction.gff")

# ADD ID IN THE SUBJECT OBJECT (EXON) -------------------------------------

# Begin empty ID column of metadata
red_exon_anot_query$ID <- NA

# Rename the numbers in increasing way
renamed_exons_numbers <- dplyr::dense_rank(df_overlaps$queryHits)

# Subset the queryHits
temp_df <- df_overlaps["queryHits"]

# Paste
temp_df$renamed <- paste0("IGHV_exon_", renamed_exons_numbers)

# Add the corresponding name of exon following the index of queryHits
cont = 1
for (i in temp_df$queryHits) {
  red_exon_anot_query[i]$ID <- temp_df$renamed[cont]
  cont = cont + 1
}

# Any remaining with NA is renamed to Sequence
red_exon_anot_query[is.na(red_exon_anot_query$ID)]$ID <- "Sequence"

# Save as gff3
rtracklayer::export.gff3(red_exon_anot_query, con = "exon_prediction.gff")