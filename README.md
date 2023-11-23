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
This script takes the filtered exonerate files (genes or exons) and perform a reduction of the overlaping sequences in order to reduce ambiguity produced from matches located at the same coordinates. The result is one genomic range per overlapping individual ranges.

### Syntax
gff_disambiguation.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The files resul |

<!---DETECT AND NAME IGHV SEGMENTS AS GENES OR PSEUDOGENES-->

## predict_ighv_by_overlaps.R
DESCRIPTION

### Syntax
gff_disambiguation.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The input file produced by blast with output format 6 (tabular). |

<!---DETECT RSS NEAR TO V AND J SEGMENTS. 
CORRECT COORDINATES.-->

## locate_nearby_rss.R
DESCRIPTION

### Syntax
gff_disambiguation.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The input file produced by blast with output format 6 (tabular). |

<!---PREDICT THE LOCATION OF TRUE D SEGMENTS-->

## predict_ighd_by_rss.R
DESCRIPTION

### Syntax
gff_disambiguation.R --file [FILE]

### Options

| Options | Description |
| --- | --- |
| file | The input file produced by blast with output format 6 (tabular). |