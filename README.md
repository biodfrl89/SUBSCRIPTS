# DESCRIPTION
This repository contains all the scripts required for the correct functioning of the master script of inmunoglobulin detection tool.

# STRUCTURE
0. Variables confirmation
1. Make CDHIT reductions.
2. Align with CLUSTALW2 the RSS reductions.
3. Make BLAST analysis.
4. Convert BLAST m6 table to gff.
5. Extract scaffolds of interest.
6. Make EXONERATE analysis.
7. Clean vulgar format.
8. Make HMMER analysis for RSS.
9. Filter EXONERATE files for exons and genes hits.
10. Reduce to eliminate redundancy in filtered EXONERATE files.
11. Make overlap analysis to detect V segments with exon and SP
12. Make MINIPROT analysis.
13. Correct V segments and RSS-J coordinates.
14. Predict D segments based in founded RSS-D.

# MAIN MODULES

<!---TRANSFORM BLAST M6 TO GFF-->

## blast_m6_to_gff.py
This script transforms blast format 6 (tabular with 11 columns) into a gff file, along with some filters.

### Syntax

blastm6_to_gff.py --file [FILE] --source [SOURCE] --bitscore [BITSCORE]

### Options

| Options | Description |
| --- | --- |
| file | The file produced by blast with output format 6 (tabular). |
| source | The mode in which blast was performed, e.g. blastn, blastp, tblastx, ... |
| bitscore | Bitscore obtained by blast. Used to add a level of filter to the result. If 0, all results will be maintained. |

<!---TRANSFORM VULGAR TO TABLE-->

## vulgar_to_table.R
This script takes the vulgar format from the exonerate analysis and transforms the vulgar syntax into a table of condensed results.

### Syntax
vulgar_to_table.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The filtered exonerate result file containing only the records with vulgar formats. |

<!---TRANSFORM HMMR TBL FORMAT TO GFF-->

## hmmer_tbl_to_gff.py 
This script transforms hmmer tbl format into a gff file.

### Syntax
hmmer_tbl_to_gff.py --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The file produced by HMMER analysis, in tbl format. |

<!---REDUCE REDUNDANCY FROM MATCHES IN THE SAME COORDINATES-->

## gff_disambiguation.R
This script takes the filtered exonerate files (genes or exons) and perform a reduction of the overlaping sequences in order to reduce 
ambiguity produced from matches located at the same coordinates. The result is one genomic range per overlapping individual ranges.

### Syntax
gff_disambiguation.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The files resul |

<!---DETECT AND NAME IGHV SEGMENTS AS GENES OR PSEUDOGENES-->

## predict_ighv_by_overlaps.R
This script takes two filtered files, the reduced exons and genes records from EXONERATE results. In order to detect if a V segment is a 
has its structural exon and its signal peptide, at least two exons need to be detected in each gene record. To do this, the number of exons per gene
is counted and every gene that has two or more exons are annotated as genes, whereas genes with one or less associated exons are annotated as
pseudogenes. Every record is numbered in a unique fashion to give each segment a unique name. 

### Syntax
SCRIPTS/SUBSCRIPTS/predict_ighv_by_overlaps.R --query [FILE]  --subject [FILE]

### Options

| Options | Description |
| --- | --- |
| query | The reduced gff file from filtered exon annotated records from EXONERATE |
| subject | The reduced gff file from filtered gene annotated records from EXONERATE |

<!---DETECT RSS NEAR TO V AND J SEGMENTS. 
CORRECT COORDINATES.-->

## locate_nearby_rss.R
This script detect which RSS is locate nearby each V/J segment. This is achieved by increasing the coordinates of RSS segments
in both the start and the end of the match, in order to try to make an artificial overlap with the nearby V/J segment

### Syntax
locate_nearby_rss.R -n [FILE] -t [FILE] -v [FILE] -r [INT] -m [STRING] 

### Options

| Options | Description |
| --- | --- |
| n | The HMMER gff file |
| t | The HMMER tbl file |
| v | The overlap gff file prediction containing the names of genes and pseudogenes |
| r | The number of positions to increase in both sides of the HMMER matches, in which to detect possible overlaps with the V segments |
| m | Mode. Depends on the input files. One of the followings: [V_segments/J_segments] |

<!---PREDICT THE LOCATION OF TRUE D SEGMENTS-->

## predict_ighd_by_rss.R
This script detects the probable location of true D segments. To do this, the reasoning was that true D segments would by flanked by D RSS signal
in both 5' and 3' ends. Every genomic region comprising both ends that have less tha 50 bp lenght is considered as a potential true D segment.

### Syntax
predict_ighd_by_rss.R -f [FILE] -t [FILE]

### Options

| Options | Description |
| --- | --- |
| f | The HMMMER gff file produced with the 5' RSS database |
| t | The HMMMER gff file produced with the 3' RSS database |

# DOCKER
The script uses many dependencies that could be annoying to install, test, and put the required locations. To prevent this, we designed a docker image in order to have 
a clean and ready-to-use environment to perform the genomic analysis using this script.  

# DEPENDENCIES
This are all the dependencies and the versions used in order to make the script to make all the analysis. Here are considered only the programs that could not be
preinstalled in a new linux installation.

| Programming Language | Library, package | Use in | Version | 
| --- | --- | --- | --- |
| Conda | --- | --- | --- |
| python | pandas | --- | --- |
| R | rtracklayer | --- | --- |
| R | GenomicRanges | --- | --- |


