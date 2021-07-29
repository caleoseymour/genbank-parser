#!/usr/bin/env Rscript

## Cale Seymour
## University of Nevada Las Vegas
## Updated August 2020
## Description: Parse a genbank formatted file into a data table, then print to
##              stdout. Does not parse ORIGIN or FEATURES sections.
## Usage: Rscript parse_genbank.R genbank_flat_file.gb > out.tsv

logger = function(msg, collapse = ' ')
{
    msg = paste0(msg, collapse = collapse)
    string = paste0('[@', Sys.time(), ']: ',msg)
    message(string)
}

## Require: stringr, stringi, data.table
args = commandArgs(trailingOnly = TRUE)

#infile = 'amoA_sequences.gb'
infile = args[1]

sentence_case = function(string)
## Take a string and make it "Sentence case."
{
    paste0(toupper(stringi::stri_sub(string, to = 1)),
    tolower(stringi::stri_sub(string, from = 2)))
}

trim_outside_whitespace = function(string)
## Remove leading and trailing whitespace on a character vector.
{
    stringi::stri_replace_all_regex(
        stringi::stri_replace_all_regex(string,
            pattern = ' *$', replacement = ''),
        pattern = '^ *', replacement = '')
}

count_leading_spaces = function(string)
{   
    nchar(string) - nchar(stringi::stri_replace_all_regex(string,
        pattern = '^ *', replacement = ''))
}

## Read in raw text
logger(c('Reading records from file:', infile))
gbtxt = data.table::fread(infile, sep='', strip.white = FALSE, header = FALSE)

## Assign a record number based on the occurrence of the record delimeter '//'.
logger('Parsing records...')
gbtxt[, record := cumsum(stringi::stri_sub(V1, to=2)  == '//') + 1]
gbtxt[stringi::stri_sub(V1, to = 9) == 'ACCESSION',
    margin :=  count_leading_spaces(stringi::stri_sub(V1, from = 10)) + 9]
gbtxt[,margin := na.omit(margin), by = record]
gbtxt[, toss := cumsum(stringi::stri_sub(V1, to=8)  == 'FEATURES') + 1]

logger('Identifying field name for each datapoint...')
datadt = gbtxt[toss == record, .(field = trim_outside_whitespace(stringi::stri_sub(V1, to = margin)),
        data = stringi::stri_sub(V1, from = margin + 1)), by = record]
headers = datadt[field != '', field]
    #headers[which(headers == 'FEATURES'):length(headers)] = NA
    datadt[, header := headers[cumsum(field != '')]]
    
logger('Parsing non-comment fields...')
GENBANK = na.omit(datadt)[header != 'COMMENT',]
GB = GENBANK[, paste0(data, collapse = ' '), by = .(sentence_case(header), record)]

logger('Parsing comment fields...')
## Comment section parsing.
COMMENT = na.omit(datadt)[header == 'COMMENT',]
CM = COMMENT[stringi::stri_count_fixed(data, '##') == 0, .(data), by=record]
groups = trim_outside_whitespace(CM
    [stringi::stri_detect_fixed(data, '::'),
        stringr::str_remove(data, '::.*$')])
CM[, datum := stringr::str_remove(data, '^.*::')]
CM[, groupid := cumsum(stringi::stri_detect_fixed(data, '::'))]
CM[groupid > 0, group := groups[groupid]]
CM[groupid == 0, group := 'Comment']
CM[, groupid := NULL]
CMDataCols = CM[,
    paste0(stringi::stri_remove_empty(trim_outside_whitespace(datum)),
    collapse = '//'), by = .(record, group)]
CMtab = data.table::dcast(CMDataCols, record ~ group, value.var='V1')

logger('Combining into a table...')
GBtab = data.table::dcast(GB, record ~ sentence_case, value.var='V1')
GBtab[,'//' := NULL]
header_order = c('record',
    sentence_case(unique(headers[headers %in% toupper(colnames(GBtab))])))
GBtabsort = GBtab[,..header_order]

## Index according to record.
data.table::setkey(GBtabsort, record)
data.table::setkey(CMtab, record)

## We want the columns from the genbank table to come first...
columns = c(colnames(GBtabsort)[-1], colnames(CMtab)[-1])

## ...But we pull indices from the comment metadata table.
out = CMtab[GBtabsort,][,-1][,..columns]

## Fix the extended whitespace in the Locus field
out[, Locus := stringi::stri_replace_all_regex(Locus, ' +',' ')]

logger('Writing to stdout...')
## Write to file.
data.table::fwrite(out, sep = '\t')
#data.table::fwrite(out, 'amoA.tsv', sep = '\t')
