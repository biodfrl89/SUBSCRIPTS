# SUBSCRIPTS
This repository contains all the scripts required for the correct functioning of the master script of inmunoglobulin detection tool.

## blast_m6_to_gff.py
This script is designed to transform blast format 6 (tabular with 11 columns) into a gff file, along with some filters.

### Syntax

blastm6_to_gff.py --file [FILE] --source [SOURCE] --bitscore [BITSCORE]

### Options

| Options | Description |
| --- | --- |
| file | The input file produced by blast with output format 6 (tabular). |
| source | The mode in which blast was performed, e.g. blastn, blastp, tblastx, ... |
| bitscore | Bitscore obtained by blast. Used to add a level of filter to the result. If 0, all results will be maintained. |

