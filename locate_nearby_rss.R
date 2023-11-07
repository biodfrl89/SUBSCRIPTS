# Name: locate_nearby_rss.R
# Description: This script is designed to detect which rss where detected in the vecinity of a exon or gene sequence.
# This is achieved by increasing the range of exon or gene both at the 5' and 3' ends of V segments, and then trying to make an overlap,
# with the detected RSS sequences from nhmmer fro V segments
# Author: David Felipe Rendón Luna 
# Date: October-2023

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

# CHECK ARGUMENTS -----------------------------------------------------
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

# SCRIPTS -----------------------------------------------

# Load files from overlap and overlaped exonerate results
gff_nhmmer <- rtracklayer::import.gff("./DATA/IGV/Nueva carpeta/nhmmer_RSS_IGHV_Roae.gff")
gff_overlap <- rtracklayer::import.gff("./DATA/IGV/Nueva carpeta/overlap_gene_prediction_IGHV_vs_Roae.gff")
gff_overlap_raw <- rtracklayer::import.gff("./DATA/IGV/Nueva carpeta/overlap_gene_prediction_IGHV_vs_Roae.gff")

# Increase the range of the exons 
gff_overlap@ranges@start <- as.integer(gff_overlap@ranges@start - 50)
gff_overlap@ranges@width <- as.integer(gff_overlap@ranges@width + 100)

# Make overlaps
overlaps <- GenomicRanges::findOverlaps(query = gff_nhmmer, subject = gff_overlap)

# Copy original df of subject to make the selection
gff_overlap_raw_mod <- gff_overlap_raw 

# Make a loop to cycle trough all the located querys
for (i in 1:length(overlaps@from)) {
  
  # Select from and to using position index
  j_query <- overlaps@from[i]
  j_subject <- overlaps@to[i]
  
  # Get all four positions (begin and end of query and subject)
  positions <- c(gff_nhmmer@ranges@start[j_query],
                 gff_nhmmer@ranges@start[j_query] + gff_nhmmer@ranges@width[j_query] -1,
                 gff_overlap_raw@ranges@start[j_subject],
                 gff_overlap_raw@ranges@start[j_subject] + gff_overlap_raw@ranges@width[j_subject] -1)
  
  # Get max and min position
  min_pos <- min(positions)
  max_pos <- max(positions)
  
  # Get the width of the range
  range_pos <- max_pos - min_pos
  
  # Security if
  if (is.na(min_pos) || is.na(range_pos) ) {
    next
  }
  
  # Replace values
  gff_overlap_raw_mod@ranges@start[j_subject] <- as.integer(min_pos)
  gff_overlap_raw_mod@ranges@width[j_subject] <- as.integer(range_pos)
}

# Save results of overlaped genes and RSS
rtracklayer::export.gff3(gff_overlap_raw_mod, "./Roae_overlaping_IGHV_RSSV.gff")

# Save only matching RSS
rtracklayer::export.gff3(gff_nhmmer[overlaps@from], "./Roae_matching_RSSV.gff")
