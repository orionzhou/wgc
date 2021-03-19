#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(argparse))

ps <- ArgumentParser(description = 'find syntenic regions and variants using crossmap')
ps$add_argument("fi", nargs=1, help="input regions (4-col BED)")
ps$add_argument("fo", nargs=1, help="output file (*.rds)")
ps$add_argument("--aln", default="/home/springer/zhoux379/projects/wgc/data/raw/Zmays_Mo17-Zmays_B73/aln.bed",
                help="synteny alignment BED [default: %(default)s]")
ps$add_argument("--vnt", default='/home/springer/zhoux379/projects/wgc/data/raw/Zmays_Mo17-Zmays_B73/vnt.bed',
                help="variant BED [default: %(default)s]")
ps$add_argument("--minsize", type="integer", default=50000,
                 help="output chain [default: %(default)s]")
args <- ps$parse_args()

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(glue))

fi = args$fi
fo = args$fo
aln = args$aln
vnt = args$vnt
if( file.access(fi) == -1 )
    stop(sprintf("Input file ( %s ) cannot be accessed", fi))

seed = sample(1:1e5,1)
fm = glue('map.{seed}.bed')
fv = glue('vnt.{seed}.bed')

system(glue("intersectBed -sorted -wao -a {fi} -b {aln} > {fm}"))
system(glue("intersectBed -sorted -wo -a {fi} -b {vnt} > {fv}"))

sum_crossmap <- function(fm, fv) {
    #{{{
    #{{{ map
    ti = read_tsv(fm,
        col_names=c('c1o','b1o','e1o','id','c1','b1','e1','srd','c2','b2','e2','cid','bp'),
        col_types = 'ciicciicciiii')
    #
    ti1 = ti %>% filter(c1 == '.') %>%
        select(id, c1=c1o, b1=b1o, e1=e1o) %>% mutate(size1=e1-b1)
    ti2 = ti %>% filter(c1 != '.') %>%
        mutate(trimU = pmax(0, b1o-b1)) %>%
        mutate(trimD = pmax(0, e1-e1o)) %>%
        mutate(b1 = b1+trimU, e1 = e1-trimD) %>%
        mutate(b2 = ifelse(srd == '-', b2+trimD, b2+trimU)) %>%
        mutate(e2 = ifelse(srd == '-', e2-trimU, e2-trimD)) %>%
        mutate(size=bp)
    ti2 %>% filter(size != e1-b1)
    ti2 %>% filter(size != e2-b2)
    ti2 %>% filter(c1o != c1)
    ti3 = ti2 %>% group_by(c1,b1o,e1o,id,c2,cid,srd) %>%
        summarise(srd = srd[1], b2=min(b2), e2=max(e2), size=sum(size)) %>%
        ungroup() %>%
        arrange(id, cid, desc(aln)) %>%
        group_by(id) %>% slice(1) %>% ungroup() %>%
        select(c1,b1=b1o,e1=e1o,id,cid,srd,size,c2,b2,e2)
    ti4 = ti2 %>% select(id,cid, rb1=b1,re1=e1,rb2=b2,re2=e2) %>%
        inner_join(ti3, by=c('id','cid')) %>%
        mutate(rb1=rb1-b1, re1=re1-b1) %>%
        mutate(rb2=rb2-b2, re2=re2-b2) %>%
        group_by(id,cid,srd,size,c1,b1,e1,c2,b2,e2) %>%
        nest() %>% rename(aln=data) %>% ungroup() %>%
        mutate(size1=e1-b1,size2=e2-b2) %>%
        select(id,c1,c2,size1,size2,size,cid,srd,b1,e1,b2,e2,aln)
    tm = ti4 %>% bind_rows(ti1) %>%
        replace_na(list(size2=0, size=0))
    tms = tm %>% select(id,cid,mb1=b1,me1=e1,mb2=b2,me2=e2)
    #}}}
    #{{{ tv
    tv = read_tsv(fv, col_names=c('c0','b0','e0','id',
                                  'c1','b1','e1','srd',
                                  'c2','b2','e2','cid','vtype',
                                  'ref','alt','bp')) %>%
        select(-c0,-b0,-e0,-bp)
    isum1 <- function(x) ifelse(is.null(x), 0, sum(x$ins))
    isum2 <- function(x) ifelse(is.null(x), 0, sum(x$del))
    inrow <- function(x) ifelse(is.null(x), 0, nrow(x))
    tvs = tv %>% filter(vtype=='snp') %>%
        select(id,cid,pos1=e1,pos2=e2,ref,alt) %>%
        inner_join(tms, by=c('id','cid')) %>%
        mutate(pos1 = pos1 - mb1, pos2 = pos2 - mb2) %>%
        select(id,cid,pos1,pos2,ref,alt) %>%
        group_by(id,cid) %>% nest() %>% rename(snp=data) %>%
        mutate(n_snp = map_dbl(snp, inrow))
    tvi = tv %>% filter(vtype=='indel') %>% select(-vtype) %>%
        inner_join(tms, by=c('id','cid')) %>%
        mutate(b1 = pmax(b1, mb1)) %>%
        mutate(e1 = pmin(e1, me1)) %>%
        mutate(b2 = pmax(b2, mb2)) %>%
        mutate(e2 = pmin(e2, me2)) %>%
        mutate(rb1 = b1 - mb1, re1 = e1 - mb1) %>%
        mutate(rb2 = b2 - mb2, re2 = e2 - mb2) %>%
        select(id, cid, rb1, re1, rb2, re2) %>%
        mutate(del=re1-rb1, ins=re2-rb2) %>%
        filter(ins+del > 0) %>%
        group_by(id,cid) %>% nest() %>% rename(idl = data) %>%
        mutate(ins = map_dbl(idl, isum1)) %>%
        mutate(del = map_dbl(idl, isum2)) %>%
        mutate(n_idl = map_dbl(idl, inrow))
    #}}}
    to = tm %>%
        left_join(tvs, by=c('id','cid')) %>%
        left_join(tvi, by=c('id','cid')) %>%
        replace_na(list(n_snp=0,n_idl=0,ins=0,del=0)) %>%
        select(id,c1,c2,size1,size2,size,ins,del,n_snp,n_idl,
               cid,srd,b1,e1,b2,e2,aln,snp,idl) %>%
        arrange(c1,b1,e1)
    to %>% filter(!is.na(c2), size1 != size + del) %>%
        print(width=Inf) %>% pluck('aln',1)
    to %>% filter(size2 != size + ins)
    to %>% filter(size1 != e1 - b1)
    to %>% filter(c1!='B99',c2!='M99') %>%
        mutate(off1 = e1-b1 - size1) %>%
        filter(off1 != 0) %>%
        print(n=10,width=Inf)
    to
    #}}}
}
to = sum_crossmap(fm, fv)
saveRDS(to, fo)
to %>% print(n=20, width=Inf)
system(glue("rm {fm} {fv}"))
