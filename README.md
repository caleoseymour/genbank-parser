# various useful tools

Various scripts and functions for bioinformatic analyses, written in R. Optimized for some degree of speed (often at the expense of resource utilization or flexibility.) Usage is typically just via the source command.



## io.align.R
Various functions for fast and simple reading of .fna and .faa files. Makes heavy use of fread and fwrite for the sake of speed.

## align.utils.R
Functions for dealing with alignment data (as read by io.align.R). Typically operates on alignments that are stored as a single-level named list or a named vector.

## parse-genbank.R
A fast genbank flat file (.gb) to tsv parser. Omits the ORIGIN and FEATURES sections.
(Bother me about it and I might add that functionality)

usage:
parse-genbank.R [input file (.gb)] > [output file (.txt)]

