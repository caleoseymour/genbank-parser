# various useful tools

Various scripts and functions for bioinformatic analyses, written in R. Optimized for some degree of speed (often at the expense of resource utilization or flexibility.) Usage is typically just via the source command.

A fast genbank flat file (.gb) to tsv parser. Omits the ORIGIN and FEATURES sections.
(Bother me about it and I might add that functionality)

usage:
parse-genbank.R [input file (.gb)] > [output file (.txt)]

