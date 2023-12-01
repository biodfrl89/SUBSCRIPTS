# Name: locate_nearby_rss.R
# Description: This script is designed to detect which rss are detected in the vecinity of a exon or gene sequence.
# This is achieved by increasing the range of exon or gene both at the 5' and 3' ends of V segments, and then trying to make an overlap,
# with the detected RSS sequences from nhmmer fro V segments.
# The coordinates of overlapping V segment and RSS are mutualy corrected.
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
opt_list = list(make_option(opt_str = c("-n", "--nhmmer"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "nhmmer gff file corresponding to the founded RSS of the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-t", "--tbl"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "nhmmer tbl file corresponding to the founded RSS of the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-v", "--overlap"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "Overlap or reduced file of exonerate gene file from the desired region (V-D-J)", 
                            metavar = "[FILENAME]"),
                make_option(opt_str = c("-r", "--range"), 
                            action="store",
                            type = "numeric", 
                            default = NULL, 
                            help = "The range of extension of the V segment in order to detect any RSS signal", 
                            metavar = "[INT]"),
                make_option(opt_str = c("-m", "--mode"), 
                            action="store",
                            type = "character", 
                            default = NULL, 
                            help = "The segments in wich the analysis is going to be performed, i.e., [VDJ]", 
                            metavar = "[STRING]")
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

if (is.null(opt$nhmmer)) {
  print_help(opt_parser)
  stop("A nhmmer file must be submitted", call.=FALSE)
}

if (is.null(opt$tbl)) {
  print_help(opt_parser)
  stop("A tbl file must be submitted", call.=FALSE)
}

if (is.null(opt$overlap)) {
  print_help(opt_parser)
  stop("An overlap file must be submitted", call.=FALSE)
}

if (is.null(opt$range)) {
  print_help(opt_parser)
  stop("A range must be submitted", call.=FALSE)
}

if (is.null(opt$mode)) {
  print_help(opt_parser)
  stop("A mode must be submitted", call.=FALSE)
} 

if (!(opt$mode %in% c("V_segments", "J_segments"))) {
  print_help(opt_parser)
  stop("Incorrect mode submitted", call.=FALSE)
} 

# SCRIPT  -----------------------------------------------------------------

# V MODE  -----------------------------------------------------------------

if (opt$mode == "V_segments" ) {
  # Load files from overlap and overlaped exonerate results
  print("ENTERING V MODE")
  gff_nhmmer <- rtracklayer::import.gff(opt$nhmmer)
  gff_overlaps <- rtracklayer::import.gff(opt$overlap)
  gff_overlaps_raw <- gff_overlaps
  hmm_tbl <- read.table(opt$tbl, 
                          col.names = c("target_name", "accession", "query_name", "accession", "hmm_from", "hmm_to", "ali_from", "ali_to", "enf_from", "env_to", 
                                        "modlen", "strand", "E-value", "score", "bias", "description_target"))
  detection_range <- opt$range

  # Increase the range of the exons 
  gff_overlaps@ranges@start <- as.integer(gff_overlaps@ranges@start - detection_range)
  gff_overlaps@ranges@width <- as.integer(gff_overlaps@ranges@width + detection_range * 2)

  # Separate plus and minus from overlap gene prediction file
  gff_overlaps_plus <- gff_overlaps[which(as.character(gff_overlaps@strand) == "+")]
  gff_overlaps_minus <- gff_overlaps[which(as.character(gff_overlaps@strand) == "-")]

  # Separate plus and minus from nhmmer gff
  gff_nhmmer_plus <- gff_nhmmer[which(as.character(gff_nhmmer@strand) == "+")]
  gff_nhmmer_minus <- gff_nhmmer[which(as.character(gff_nhmmer@strand) == "-")]

  # Separate plus and minus from nhmmer tbl
  hmm_tbl_plus <- hmm_tbl[which(hmm_tbl$strand == "+"),]
  hmm_tbl_minus <- hmm_tbl[which(hmm_tbl$strand == "-"),]

  # Make overlaps in both strands
  overlaps_plus <- GenomicRanges::findOverlaps(query = gff_nhmmer_plus, subject = gff_overlaps_plus)
  overlaps_minus <- GenomicRanges::findOverlaps(query = gff_nhmmer_minus, subject = gff_overlaps_minus)

  # Checkpoint for overlaps and make analysis in plus strand
  if (length(overlaps_plus) < 1) {
    print("No overlaps detected for plus strand")
  } else {
    # Copy original df of subject to make the selection
    gff_overlaps_plus <- gff_overlaps_raw
    gff_overlaps_plus <- gff_overlaps_plus[which(as.character(gff_overlaps_plus@strand) == "+")]
    
    # Make a loop to cycle trough all the located querys
    for (i in seq_along(overlaps_plus)) {
      # Select from and to using position index
      j_query <- overlaps_plus@from[i]
      j_subject <- overlaps_plus@to[i]
      
      # Calculate RSS real begining for the record in nhhmer query
      rss_real_begin <- hmm_tbl_plus[j_query, "ali_from"] + 1 - hmm_tbl_plus[j_query, "hmm_from"] 
      
      # Replace RSS real begining at the start of RSS
      gff_nhmmer_plus@ranges@start[j_query] <- as.integer(rss_real_begin)
      
      # Calculate the real RSS end
      #rss_real_end <- rss_real_begin + hmm_tbl[j_query,"hmm_to"] + (39 - hmm_tbl[j_query,"hmm_to"]) - 1
      rss_real_end <- rss_real_begin + hmm_tbl[j_query, "modlen"] 
      
      # Replace width for gff_nhmmer
      gff_nhmmer_plus@ranges@width[j_query] <- as.integer(rss_real_end - rss_real_begin) 
      
      # Replace RSS real begining at the end of V segment
      v_start <- gff_overlaps_plus[j_subject]@ranges@start
      
      v_new_width <- abs(rss_real_begin - v_start)
      
      gff_overlaps_plus[j_subject]@ranges@width <- as.integer(v_new_width)
    }

    # Initialice column of metadata for RSS
    gff_nhmmer_plus$ID <- NA
    # Fill metadata
    if (length(gff_nhmmer_plus[overlaps_plus@from]) >= 1) gff_nhmmer_plus[overlaps_plus@from]$ID <- paste0("RSS_V-associated_", seq_along(gff_nhmmer_plus[overlaps_plus@from]))
    if (length(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]) >= 1) gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]$ID <- paste("RSS_V-non-associated_", seq_along(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]))
    
    # Select V segments with nearby RSS signals, with corrected coordinates from plus strand
    final_plus <- c(gff_overlaps_plus[overlaps_plus@to], gff_nhmmer_plus[overlaps_plus@from])
  }

  # Checkpoint for overlaps and make analysis in minus strand
  if (length(overlaps_minus) < 1) {
    print("No overlaps detected for minus strand")
  } else {
    # Copy original df of subject to make the selection
    gff_overlaps_minus <- gff_overlaps_raw
    gff_overlaps_minus <- gff_overlaps_minus[which(as.character(gff_overlaps_minus@strand) == "-")]
    
    # Make a loop to cycle trough all the located querys
    for (i in seq_along(overlaps_minus)) {
      # Select from and to using position index
      j_query <- overlaps_minus@from[i]
      j_subject <- overlaps_minus@to[i]
      
      # Calculate RSS real end for the record in nhhmer query
      rss_real_end <- hmm_tbl_minus[j_query, "ali_from"] + 1 - hmm_tbl_minus[j_query, "hmm_from"] 
      
      # Calculate RSS real begin according to length
      rss_real_begin <- rss_real_end - hmm_tbl_minus[j_query, "modlen"] + 1
      
      # Replace RSS real begining in the record
      gff_nhmmer_minus@ranges@start[j_query] <- as.integer(rss_real_begin)
      
      # Replace RSS real width in in the record
      gff_nhmmer_minus@ranges@width[j_query] <- as.integer(hmm_tbl_minus[j_query, "modlen"])
      
      # Replace V segment start with new location
      gff_overlaps_minus[j_subject]@ranges@start <- as.integer(rss_real_end + 1)
    }

    # Initialice column of metadata for RSS
    gff_nhmmer_minus$ID <- NA
    # Fill metadata
    if (length(gff_nhmmer_minus[overlaps_minus@from]) >= 1) gff_nhmmer_minus[overlaps_minus@from]$ID <- paste0("RSS_J-associated_", seq_along(gff_nhmmer_minus[overlaps_minus@from]))
    if (length(gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]) >= 1) gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]$ID <- paste("RSS_J-non-associated_", seq_along(gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]))
    
    # Select V segments with nearby RSS signals, with corrected coordinates from minus strand
    final_minus <- c(gff_overlaps_minus[overlaps_minus@to], gff_nhmmer_minus[overlaps_minus@from])
  }

# Export gff with all the final data
rtracklayer::export.gff3(c(final_minus, final_plus), "v_rss_analysis.gff")
}

# J MODE  -----------------------------------------------------------------

if (opt$mode == "J_segments" ) {
  print("ENTERING J MODE")
  # Load files from overlap and overlaped exonerate results
  gff_nhmmer <- rtracklayer::import.gff(opt$nhmmer)
  gff_overlaps <- rtracklayer::import.gff(opt$overlap)
  gff_overlaps_raw <- gff_overlaps
  hmm_tbl <- read.table(opt$tbl, 
                          col.names = c("target_name", "accession", "query_name", "accession", "hmm_from", "hmm_to", "ali_from", "ali_to", "enf_from", "env_to", 
                                        "modlen", "strand", "E-value", "score", "bias", "description_target"))
  detection_range <- opt$range

  # Increase the range of the exons 
  gff_overlaps@ranges@start <- as.integer(gff_overlaps@ranges@start - detection_range)
  gff_overlaps@ranges@width <- as.integer(gff_overlaps@ranges@width + detection_range * 2)

  # Separate plus and minus from overlap gene prediction file
  gff_overlaps_plus <- gff_overlaps[which(as.character(gff_overlaps@strand) == "+")]
  gff_overlaps_minus <- gff_overlaps[which(as.character(gff_overlaps@strand) == "-")]

  # Separate plus and minus from nhmmer gff
  gff_nhmmer_plus <- gff_nhmmer[which(as.character(gff_nhmmer@strand) == "+")]
  gff_nhmmer_minus <- gff_nhmmer[which(as.character(gff_nhmmer@strand) == "-")]

  # Separate plus and minus from nhmmer tbl
  hmm_tbl_plus <- hmm_tbl[which(hmm_tbl$strand == "+"),]
  hmm_tbl_minus <- hmm_tbl[which(hmm_tbl$strand == "-"),]

  # Make overlaps in both strands
  overlaps_plus <- GenomicRanges::findOverlaps(query = gff_nhmmer_plus, subject = gff_overlaps_plus)
  overlaps_minus <- GenomicRanges::findOverlaps(query = gff_nhmmer_minus, subject = gff_overlaps_minus)

  ### PLUS STRAND ANALYSIS ###
  if (length(overlaps_plus) < 1) {
    print("No overlaps detected for plus strand")
    final_plus <- c()
  } else {
    print("ENTERING PLUS")
    # Copy original df of subject to make the selection
    gff_overlaps_plus <- gff_overlaps_raw
    gff_overlaps_plus <- gff_overlaps_plus[which(as.character(gff_overlaps_plus@strand) == "+")]
    
    # Make a loop to cycle trough all the located querys
    for (i in seq_along(overlaps_plus)) {
      # Select from and to using position index
      j_query <- overlaps_plus@from[i]
      j_subject <- overlaps_plus@to[i]
      
      # Calculate RSS real begining for the record in nhhmer query
      rss_real_begin <- hmm_tbl_plus[j_query, "ali_from"] + 1 - hmm_tbl_plus[j_query, "hmm_from"] 
      
      # Replace RSS real begining at the start of RSS
      gff_nhmmer_plus@ranges@start[j_query] <- as.integer(rss_real_begin)
      
      # Calculate the real RSS end
      #rss_real_end <- rss_real_begin + hmm_tbl[j_query,"hmm_to"] + (39 - hmm_tbl[j_query,"hmm_to"]) - 1
      rss_real_end <- rss_real_begin + hmm_tbl[j_query, "modlen"] 
      
      # Replace width for gff_nhmmer
      gff_nhmmer_plus@ranges@width[j_query] <- as.integer(rss_real_end - rss_real_begin) 
      
      # Calculate the diferrence in width
      v_start <- gff_overlaps_plus[j_subject]@ranges@start
      v_new_width_dif <- abs(rss_real_end - v_start)
      
      # Modify width
      gff_overlaps_plus[j_subject]@ranges@width <- as.integer(gff_overlaps_plus[j_subject]@ranges@width - v_new_width_dif)
      
      # Replace V segment start
      gff_overlaps_plus[j_subject]@ranges@start <- as.integer(rss_real_end)
    }
    
    # Initialice column of metadata for J segments
    gff_overlaps_plus$ID <- NA
    # Fill metadata
    if(length(gff_overlaps_plus[overlaps_plus@to]) >= 1) gff_overlaps_plus[overlaps_plus@to]$ID <- paste0("IGHJ_gene_", seq_along(gff_overlaps_plus[overlaps_plus@to]))
    if(length(gff_overlaps_plus[which(is.na(gff_overlaps_plus$ID))]) >= 1) gff_overlaps_plus[which(is.na(gff_overlaps_plus$ID))]$ID <- paste0("IGHJ_pseudogene_", seq_along(gff_overlaps_plus[which(is.na(gff_overlaps_plus$ID))]))
    
    # Initialice column of metadata for RSS
    gff_nhmmer_plus$ID <- NA
    # Fill metadata
    if (length(gff_nhmmer_plus[overlaps_plus@from]) >= 1) gff_nhmmer_plus[overlaps_plus@from]$ID <- paste0("RSS_J-associated_", seq_along(gff_nhmmer_plus[overlaps_plus@from]))
    if (length(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]) >= 1) gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]$ID <- paste("RSS_J-non-associated_", seq_along(gff_nhmmer_plus[which(is.na(gff_nhmmer_plus$ID))]))
    
    # Select V segments with nearby RSS signals, with corrected coordinates from plus strand
    #final_plus <- c(gff_overlaps_plus[overlaps_plus@to], gff_nhmmer_plus[overlaps_plus@from])
    final_plus <- c(gff_overlaps_plus, gff_nhmmer_plus)
  }
  
  ### MINUS STRAND ANALYSIS ###
  if (length(overlaps_minus) < 1) {
    print("No overlaps detected for minus strand")
    final_minus <- c()
  } else {
    print("ENTERING MINUS")
    # Copy original df of subject to make the selection
    gff_overlaps_minus <- gff_overlaps_raw
    gff_overlaps_minus <- gff_overlaps_minus[which(as.character(gff_overlaps_minus@strand) == "-")]
    
    # Make a loop to cycle trough all the located querys
    for (i in seq_along(overlaps_minus)) {
      # Select from and to using position index
      j_query <- overlaps_minus@from[i]
      j_subject <- overlaps_minus@to[i]
      
      # Calculate RSS real end for the record in nhhmer query
      rss_real_end <- hmm_tbl_minus[j_query, "ali_from"] + 1 - hmm_tbl_minus[j_query, "hmm_from"] 
      
      # Calculate RSS real begin according to length
      rss_real_begin <- rss_real_end - hmm_tbl_minus[j_query, "modlen"] + 1
      
      # Replace RSS real begining in the record
      gff_nhmmer_minus@ranges@start[j_query] <- as.integer(rss_real_begin)
      
      # Replace RSS real width in in the record
      gff_nhmmer_minus@ranges@width[j_query] <- as.integer(hmm_tbl_minus[j_query, "modlen"])
      
      # Replace V segment start with new location
      gff_overlaps_minus[j_subject]@ranges@start <- as.integer(rss_real_begin - gff_overlaps_minus[j_subject]@ranges@width + 1)
      
      # Calculate the distance width between the real start of the RSS and the start of v segment
      v_new_width <- (rss_real_begin - gff_overlaps_minus[j_subject]@ranges@start) #- gff_overlaps_minus[j_subject]@ranges@start
      
      # Replace J segment width
      gff_overlaps_minus[j_subject]@ranges@width <- as.integer(v_new_width)
    }
    
    # Initialice column of metadata for J segments
    gff_overlaps_minus$ID <- NA
    # Fill metadata
    if(length(gff_overlaps_minus[overlaps_minus@to]) >= 1) gff_overlaps_minus[overlaps_minus@to]$ID <- paste0("IGHJ_gene_", seq_along(gff_overlaps_minus[overlaps_minus@to]))
    if(length(gff_overlaps_minus[which(is.na(gff_overlaps_minus$ID))]) >= 1) gff_overlaps_minus[which(is.na(gff_overlaps_minus$ID))]$ID <- paste0("IGHJ_pseudogene_", seq_along(gff_overlaps_minus[which(is.na(gff_overlaps_minus$ID))]))
    
    # Initialice column of metadata for RSS
    gff_nhmmer_minus$ID <- NA
    # Fill metadata
    if (length(gff_nhmmer_minus[overlaps_minus@from]) >= 1) gff_nhmmer_minus[overlaps_minus@from]$ID <- paste0("RSS_J-associated_", seq_along(gff_nhmmer_minus[overlaps_minus@from]))
    if (length(gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]) >= 1) gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]$ID <- paste("RSS_J-non-associated_", seq_along(gff_nhmmer_minus[which(is.na(gff_nhmmer_minus$ID))]))

    # Merge data
    final_minus <- c(gff_overlaps_minus, gff_nhmmer_minus)
  }
rtracklayer::export.gff3(c(final_minus, final_plus), "j_rss_analysis.gff")
}

