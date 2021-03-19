#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(tidyverse))

ps <- ArgumentParser(description = 'wrapper to liftover btw. genome assemblies')
ps$add_argument("fi", help="BED file in query genome coordinate")
#ps$add_argument("fc", help="input chain file")
ps$add_argument("fo", help="output summary tsv")
ps$add_argument("--qry", default='Zmays_Mo17', help="query genome [default: %(default)s]")
ps$add_argument("--tgt", default='Zmays_B73',  help="target genome [default: %(default)s]")
ps$add_argument("--diri", default='/home/springer/zhoux379/projects/wgc/data/raw', help="directory containing genome comparison files (chain/vcf) [default: %(default)s]")
args <- ps$parse_args()

fi = args$fi
fo = args$fo
diri = args$diri; qry = args$qry; tgt = args$tgt
t2q = F
if(qry %in% c("Zmays_B73")) {
    qry = args$tgt
    tgt = args$qry
    t2q = T
}
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))
dirw = sprintf("%s/%s-%s", diri, qry, tgt)
fc = file.path(dirw, ifelse(t2q, '10.t2q.chain', '10.q2t.chain'))
if( file.access(fc) == -1 )
    stop(sprintf("Chain file ( %s ) cannot be accessed", fc))
fv = file.path(dirw, ifelse(t2q, "11.t.vnt.bed", '11.q.vnt.bed'))
if( file.access(fv) == -1 )
    stop(sprintf("Variant file ( %s ) cannot be accessed", fv))

sum_liftover <- function(f_map, f_vnt, max_indel=Inf) {
    #{{{
    #{{{ read in
    ti = read_tsv(f_map, col_names=c('chrom1','beg1','end1','id','opt','chrom2','beg2','end2','id2')) %>%
        select(-id2)
    ti0 = ti %>% distinct(id, chrom1,beg1,end1) %>% mutate(size1=end1-beg1)
    ti = ti %>% select(-chrom1,-beg1,-end1)
    tv = read_tsv(f_vnt, col_names=c('chrom','beg','end','id',
                                     'chrom1','beg1','end1','srd',
                                     'chrom2','beg2','end2','cid',
                                     'type','ref','alt','bp')) %>%
        select(-chrom,-beg,-end,-bp) %>%
        group_by(id) %>% nest() %>% rename(vnt=data)
    #}}}
    #
    tm1 = ti %>% filter(opt %in% c("Unmap","->")) %>%
        mutate(opt = ifelse(opt == '->', 1, 0)) %>%
        rename(n_block=opt) %>% mutate(aln=end2-beg2,size2=end2-beg2)
    ti2 = ti %>% filter(!opt %in% c("Unmap","->"))
    #
    tm2 = ti2 %>% mutate(opt = str_replace_all(opt, "[()]", '')) %>%
        separate(opt, c('opt','chrom1','beg1','end1','srd'), sep=":", fill='right') %>%
        mutate(beg1 = as.double(beg1), end1=as.double(end1)) %>%
        select(id,chrom1,beg1,end1,chrom2,beg2,end2)
    tm2s = tm2 %>% group_by(id) %>%
        summarise(n_block=n(), aln1=sum(end1-beg1), aln2=sum(end2-beg2),
                  chrom2=chrom2[1], beg2=min(beg2), end2=max(end2)) %>%
        ungroup() %>% mutate(size2=end2-beg2)
    stopifnot(sum(tm2s$aln1-tm2s$aln2)==0)
    tm2s = tm2s %>% rename(aln=aln1) %>% select(-aln2)
    #
    tm3a = tm2 %>% group_by(id) %>% slice(-n()) %>% ungroup()
    tm3b = tm2 %>% group_by(id) %>% slice(-1) %>% ungroup() %>%
        select(id, chrom2b=chrom2, beg1b=beg1, end1b=end1, beg2b=beg2, end2b=end2)
    identical(tm3a$id, tm3b$id)
    tm3 = tm3a %>% bind_cols(tm3b[,-1]) %>%
        mutate(srd = ifelse(beg2b < beg2, '-', '+'))
    ids_rm = tm3 %>%
        mutate(inchain = (chrom2b==chrom2 & abs(beg2b-beg2) < max_indel)) %>%
        filter(!inchain) %>% distinct(id) %>% pull(id)
    tm3 = tm3 %>% filter(!id %in% ids_rm)
    #
    tm4a = tm3 %>% filter(srd == "+") %>%
        mutate(is_del = beg1b!=end1, is_ins = beg2b!=end2,
               del=beg1b-end1, ins=beg2b-end2)
    tm4b = tm3 %>% filter(srd == "-") %>%
        mutate(is_del = beg1b!=end1, is_ins = beg2!=end2b,
               del=beg1b-end1, ins=beg2-end2b)
    tm4 = tm4a %>% bind_rows(tm4b) %>%
        group_by(id) %>%
        summarise(n_srd=length(unique(srd)), n_del=sum(is_del), n_ins=sum(is_ins),
            del = sum(del), ins = sum(ins)) %>% ungroup()
    stopifnot(sum(tm4$n_srd != 1)==0)
    tm5 = tm4 %>% select(-n_srd)
    #
    tm = tm2s %>% left_join(tm5, by='id') %>%
        mutate(n_block = ifelse(id %in% ids_rm, 0, n_block)) %>%
        bind_rows(tm1) %>% arrange(id) %>%
        replace_na(list(n_del=0,n_ins=0,del=0,ins=0))
    to = tm %>% left_join(tv, by='id') %>%
        inner_join(ti0, by='id') %>%
        group_by(id,size1,size2,n_block,aln,n_del,n_ins,del,ins,vnt) %>%
        nest() %>% rename(coord = data) %>% ungroup() %>%
        select(id,size1,size2,aln,n_block,n_del,n_ins,del,ins,coord,vnt)
    to
    #}}}
}

cmd = sprintf("CrossMap.py bed %s %s > tmp1.bed", fc, fi)
system(cmd)

cmd = sprintf("intersectBed -wo -sorted -a %s -b %s > tmp2.bed", fi, fv)
system(cmd)

x = sum_liftover('tmp1.bed', 'tmp2.bed')
saveRDS(x, fo)
