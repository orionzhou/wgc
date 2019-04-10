#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(BBmisc))

parser <- ArgumentParser(description = 'Split BED file into N pieces with equal total size then extract fasta')
parser$add_argument("bed", nargs=1, help="Input (BED) file")
parser$add_argument("db", nargs=1, help="fasta file to extract sequence from")
parser$add_argument("dbsize", nargs=1, help="fasta size file")
parser$add_argument("outdir", nargs=1, help="Output directory")
parser$add_argument("-n", "--N", type="integer", default=10,
                    help="num. piecies to split [default: %(default)s]")
parser$add_argument("--chain", default='out.chain',
                    help="output chain [default: %(default)s]")
args <- parser$parse_args()

fi = args$bed
db = args$db
dbsize = args$dbsize
outdir = args$outdir
n = args$N
chain = args$chain
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))
if( !dir.exists(outdir) )
    dir.create(outdir, showWarnings = T, recursive = T)

source("~/projects/wgc/src/functions.R")
ti = read_tsv(fi, col_names = F, col_types = 'cii') %>%
    transmute(chrom = X1, start = X2, end = X3)

tp = ti %>%
    mutate(size = end - start,
           tchrom = sprintf("%s-%d-%d", chrom, start+1, end),
           tstart = 0, tend = size, srd = '+',
           cid = 1:length(size),
           bin = binPack(size, sum(as.numeric(size))/n + 1))

bedsize = file.path(outdir, 'bed.sizes')
to = tp %>% select(tchrom, size)
write_tsv(to, bedsize, col_names = F)

chainbed = file.path(outdir, 'chain.bed')
to = tp %>% select(tchrom, tstart, tend, srd, chrom, start, end, cid)
#to = tp %>% select(chrom, start, end, srd, tchrom, tstart, tend, cid)
write_tsv(to, chainbed, col_names = F)

cmd = sprintf("chain.py fromBed %s %s %s > %s", chainbed, bedsize, dbsize, chain)
system(cmd)
system(sprintf("rm %s %s", chainbed, bedsize))

bed = file.path(outdir, 'tmp.bed')
for (i in 1:n) {
    to = tp %>% filter(bin == i) %>%
        select(chrom, start, end)
    write_tsv(to, bed, col_names = F)
    fo = sprintf("%s/part.%d.fna", outdir, i)
    cmd = sprintf("fasta.py extract %s %s > %s", db, bed, fo)
    system(cmd)
}
system(sprintf("rm %s", bed))
