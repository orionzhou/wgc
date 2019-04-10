#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(tidyverse))

parser <- ArgumentParser(description = 'rename blast (tabular) output seqids')
parser$add_argument("fi", nargs=1, help="blast tabular otuput")
parser$add_argument("fm", nargs=1, help="seqid mapping file")
parser$add_argument("fo", nargs=1, help="output tabular file")
args <- parser$parse_args()

fi = args$fi
fm = args$fm
fo = args$fo
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))
if( file.access(fm) == -1 )
    stop(sprintf("Mapping file ( %s ) cannot be accessed", fm))

cols12 = c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
ti = read_tsv(fi, col_names = cols12)

tm = read_tsv(fm, col_names = 'col1') %>%
    separate(col1, c('oid','nid'), sep = ': ')

to = ti %>% inner_join(tm, by = c('qseqid'='oid')) %>%
    mutate(qseqid = nid) %>% select(-nid) %>%
    inner_join(tm, by = c('sseqid' = 'oid')) %>%
    mutate(sseqid = nid) %>% select(-nid)
write_tsv(to, fo, col_names = F)

