# Name: vulgar_to_table.R
# Description: 
# Author: David Felipe Rendón Luna 
# Date: Septembre-2023


# CHECK LIBRARIES ----------

# Check for optparse library to load arguments from command line
if(suppressMessages(!require("optparse"))) {
  stop("optparse was not found. Exiting.")
}

# LOAD PARSER AND OPTIONS ----------

# Load parser
library("optparse")
# Create list of arguments and help asociadted to receive.
opt_list = list(make_option(opt_str = c("-f", "--file"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Exonerate vulgar anotation file", 
                            metavar = "[FILENAME]")
                )

# Make the parser
opt_parser = OptionParser(option_list = opt_list)

# Load the arguments in the parser into an object
opt = parse_args(opt_parser)

# CHECK ARGUMENTS -----------------------------------------------------
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

# Read vulgar data
df <- read.delim(opt$file, sep = ";", header = FALSE, strip.white = TRUE, 
                 col.names = c("Query", "Query_begin", "Query_end", "Query_strand", "Target", 
                               "Target_begin", "Target_end", "Target_strand", "Raw_score", "Vulgar_notation") )

# Create df for deposit the counts
df_final <- data.frame("Match" = numeric(),
           "Codon" = numeric(),
           "Gap" = numeric(),
           "Non_equivalent" = numeric(),
           "Splice_5" = numeric(),
           "Intron" = numeric(),
           "Splice_3" = numeric(),
           "Split" = numeric(),
           "Frameshift" = numeric())

for (COUNT in 1:nrow(df)) {
  # Split last column into individual elements
  vec_all_sections <- strsplit(df$Vulgar_notation[COUNT], " ")[[1]]

  # Get only the elements multiple of 3
  vec_sections <- vec_all_sections[seq(1, length(vec_all_sections), 3)]

  # Convert to factor with levels predefined
  fac_sections <- factor(vec_sections, levels = (c("M", "C", "G", "N", "5" ,"I", "3" , "S", "F")), ordered = TRUE)

  # Deposit counts of all categories in the final df
  df_final[COUNT,] <- unname(table(fac_sections))
}

# Get path of the file
in_path <- opt$file

# Separate by slash, and get the last element, the filename
out_filename_old <- strsplit(in_path, "/")[[1]][5]

# Separate by underscore, remove first element and replace it with "table"
out_filename_new <- c("table", strsplit(out_filename_old, "_")[[1]][-1])

# Collapse the previous elements with underscore
out_filename_new <- paste(out_filename_new, collapse = "_")

# From the original path, split, remove last element and replace it with the new filename. Collapse with slash
out_path <- paste(c(strsplit(in_path, "/")[[1]][-5], out_filename_new), collapse = "/" )

# Write the final dataframe in the created destination
write.table(cbind(df[,1:9], df_final), 
            file = out_path , quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


