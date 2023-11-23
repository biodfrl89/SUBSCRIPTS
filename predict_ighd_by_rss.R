# Name: predict_ighd_by_rss.R
# Description: This script is designed to predict IGHD genes location 
# using the location of the 5´ and 3' RSS for D segments
# Author: David Felipe Rendón Luna 
# Date: November-2023

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
opt_list = list(make_option(opt_str = c("-f", "--five"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "hmmer gff file for 5' RSS", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-t", "--three"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "hmmer gff file for 3' RSS", 
                            metavar = "[FILENAME]"))

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# CHECK FOR ARGUMENTS -----------------------------------------------------
if (length(opt) == 1) {
  print_help(opt_parser)
  stop("No arguments submitted")
}

if (is.null(opt$five)) {
  print_help(opt_parser)
  stop("A hmmer gff file of the 5' RSS must be submitted", call.=FALSE)
}

if (is.null(opt$three)) {
  print_help(opt_parser)
  stop("A hmmer gff file of the 3' RSS must be submitted", call.=FALSE)
}

# REASSIGNMENTS -------------------------------------------------------------

five_gff <- opt$five
three_gff <- opt$three

# SCRIPT ------------------------------------------------------------------

# Read 5' RSS
gr_5 <- rtracklayer::import(five_gff, format = "gff")

# Read 3' RSS
gr_3 <- rtracklayer::import(three_gff , format = "gff")

# Merge the results from both gff
gr_merged <- c(gr_5, gr_3)
gr_merged <- setNames(gr_merged, rep("RSS_signal",length(gr_merged)))

# Make reduction ignoring strand
gr_reduced <- GenomicRanges::reduce(gr_merged, ignore.strand = TRUE)

# Sort the genomic ranges object. This is to mix the ranges without considering the strand
gr_reduced_sorted <- gr_reduced[order(gr_reduced@ranges@start)]

# Calculate the distances between the end of one element and the start of the next
distances <- (gr_reduced_sorted[2:length(gr_reduced_sorted)]@ranges@start) - (gr_reduced_sorted@ranges@start + gr_reduced_sorted@ranges@width) 

# Detect the distances between 0 and 50 bp
dist_log <- distances < 50 & distances > 0

# Extract the detected segment lengths
length_segments <- distances[dist_log]

# Extract the records asociated with such distances
gr_segments <- gr_reduced_sorted[dist_log]

# Create a new GR object. The seqname correspond to the detected scaffold. 
# Iranges is constructed considering the end of one RSS and the same but adding the distance of the D segment, minus one position
D_segmentes_gr <- GenomicRanges::GRanges(
  seqnames = rep("MW800879.1", length(gr_segments)),
  ranges = IRanges::IRanges(gr_segments@ranges@start + gr_segments@ranges@width, 
                   end = gr_segments@ranges@start + gr_segments@ranges@width + distances[dist_log] -1, 
                   names = paste0("D_segment_", 1:length(gr_segments))),
  source = "predicted",
  type = "D_segment", 
  score = ".", 
  phase = ".",
  Query = "Predicted_D_segment")

# Merge RSS data with predicted D_segments
final_merged <- c(gr_merged, D_segmentes_gr)

# Save the resulting GR object as gff
rtracklayer::export(final_merged , "RSS_D_IGH_D_segments.gff", format = "gff3") 