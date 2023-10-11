# Name: gff_disambiguation.R
# Description: This script is designed to eliminate the redundancy in the results obtained
# by BLAST and EXONERATE, presesnted in GFF format. The output is also a GFF file, with no
# overlapping sequences.
# Author: David Felipe Rendón Luna 
# Date: July-2023


# Check for optparse library to load arguments from command line ----------
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

# Load parser -------------------------------------------------------------
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-f", "--file"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "GFF files to reduce the ambiguity", 
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

if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("A file must be submitted", call.=FALSE)
}

# SCRIPT
library(utils)

# CREATE OUTPUT
filename <- opt$file
out_filename <- paste0("reduced_", strsplit(filename, "/")[[1]][5])
out_path <- paste0(strsplit(filename, "/")[[1]][1:4], collapse = "/")
out_path_complete <- paste(c(out_path, out_filename), collapse = "/")

# Load GFF file using the name provided in the file option
print("Reading GFF and transforming to Granges.")
gff_file <- rtracklayer::import.gff(opt$file)


# Reduce the GFF
print("Reducing ranges.")
gff_reduced <- GenomicRanges::reduce(gff_file )

# Saving the reduced GFF
print("Saving reduced GFF.")
rtracklayer::export.gff(gff_reduced, out_path_complete)

# Load and resave to eliminate first rows
readr::write_lines(readr::read_lines(out_path_complete, skip = 3),
                   out_path_complete)
