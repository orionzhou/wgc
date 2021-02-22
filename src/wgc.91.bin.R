source('functions.R')
dirw = glue('{dird}/91_genome_bins')
#{{{ functions
sum_liftover <- function(f_map, f_snp, max_indel=1e5) {
    #{{{
    #{{{ read in
    ti = read_tsv(f_map, col_names=c('chrom1','beg1','end1','bin','opt','chrom2','beg2','end2','bin2')) %>%
        select(-chrom1,-beg1,-end1,-bin2)
    tv = read_tsv(f_snp, col_names=c('chrom','beg','end','bin','n_snp')) %>%
        select(bin,n_snp)
    #}}}
    #
    tm1 = ti %>% filter(opt %in% c("Unmap","->")) %>%
        mutate(opt = ifelse(opt == '->', 1, 0)) %>%
        rename(n_block=opt)
    ti2 = ti %>% filter(!opt %in% c("Unmap","->"))
    #
    tm2 = ti2 %>% mutate(opt = str_replace_all(opt, "[()]", '')) %>%
        separate(opt, c('opt','chrom1','beg1','end1','srd'), sep=":", fill='right') %>%
        mutate(beg1 = as.double(beg1), end1=as.double(end1)) %>%
        select(bin,chrom1,beg1,end1,chrom2,beg2,end2)
    tm2s = tm2 %>% group_by(bin) %>%
        summarise(n_block=n(), aln1=sum(end1-beg1), aln2=sum(end2-beg2),
                  chrom2=chrom2[1], beg2=min(beg2), end2=max(end2)) %>%
        ungroup() %>% mutate(size2=end2-beg2)
    stopifnot(sum(tm2s$aln1-tm2s$aln2)==0)
    tm2s = tm2s %>% rename(aln=aln1) %>% select(-aln2)
    #
    tm3a = tm2 %>% group_by(bin) %>% slice(-n()) %>% ungroup()
    tm3b = tm2 %>% group_by(bin) %>% slice(-1) %>% ungroup() %>%
        select(bin, chrom2b=chrom2, beg1b=beg1, end1b=end1, beg2b=beg2, end2b=end2)
    identical(tm3a$bin, tm3b$bin)
    tm3 = tm3a %>% bind_cols(tm3b[,-1]) %>%
        mutate(srd = ifelse(beg2b < beg2, '-', '+'))
    bins_rm = tm3 %>%
        mutate(inchain = (chrom2b==chrom2 & abs(beg2b-beg2) < max_indel)) %>%
        filter(!inchain) %>% distinct(bin) %>% pull(bin)
    tm3 = tm3 %>% filter(!bin %in% bins_rm)
    #
    tm4a = tm3 %>% filter(srd == "+") %>%
        mutate(is_del = beg1b!=end1, is_ins = beg2b!=end2,
               del=beg1b-end1, ins=beg2b-end2)
    tm4b = tm3 %>% filter(srd == "-") %>%
        mutate(is_del = beg1b!=end1, is_ins = beg2!=end2b,
               del=beg1b-end1, ins=beg2-end2b)
    tm4 = tm4a %>% bind_rows(tm4b) %>%
        group_by(bin) %>%
        summarise(n_srd=length(unique(srd)), n_del=sum(is_del), n_ins=sum(is_ins),
            del = sum(del), ins = sum(ins)) %>% ungroup()
    stopifnot(sum(tm4$n_srd != 1)==0)
    tm5 = tm4 %>% select(-n_srd)
    #
    tm = tm2s %>% left_join(tm5, by='bin') %>%
        mutate(n_block = ifelse(bin %in% bins_rm, 0, n_block)) %>%
        bind_rows(tm1) %>% arrange(bin)
    stopifnot(identical(tm$bin, tv$bin))
    to = tm %>% bind_cols(tv[,-1]) %>%
        select(bin,chrom2,beg2,end2,size2,n_block,aln,n_snp,n_del,n_ins,del,ins)
    to
    #}}}
}
sum_gene_itv <- function(fb, fg, min_fraction=.3) {
    #{{{
    fi = sprintf("tmp.bed")
    cmd = sprintf("intersectBed -wo -a %s -b %s > %s", fb, fg, fi)
    system(cmd)
    ti = read_tsv(fi, col_names=c('chrom','start','end','bin','b.chrom','b.start','b.end','ttype','etype','bp')) %>%
        mutate(size=end-start) %>%
        group_by(bin,size,ttype,etype) %>%
        summarise(bp=sum(bp)) %>% ungroup() %>%
        arrange(bin, desc(bp), etype) %>%
        group_by(bin) %>% slice(1) %>% ungroup() %>%
        filter(bp >= size * min_fraction) %>%
        mutate(type = str_c(ttype,etype, sep='_')) %>%
        select(bin, type)
    #
    tb = read_tsv(fb, col_names=c("chrom",'start','end','bin')) %>%
        select(bin) %>%
        left_join(ti, by='bin') %>%
        replace_na(list(type='intergenic'))
    tb
    #}}}
}
gts = tibble(gt0 = c('B73','Mo17','W22','PH207','Oh43'), gt=c('B','M','W','P','O')) %>%
    mutate(gt0 = str_c("Zmays", gt0, sep='_'))
comps = tibble(tgt = "B", qry = c('M','W','P','O'))
#}}}

#{{{ create 100-bp tiles
for (i in 5) {
gt0 = gts$gt0[i]; gt = gts$gt[i]
cfg = read_genome_conf(gt0)
tw = cfg$chrom %>% mutate(start=1) %>% select(chrom, start, end=size) %>%
    mutate(res = map2(start, end, make_tile, winsize=100, winstep=100)) %>%
    select(chrom, res) %>% unnest(res) %>%
    mutate(size = end-beg+1) %>%
    filter(size == 100) %>%
    mutate(bin = sprintf("%s_%07d", chrom, ceiling(beg/100))) %>%
    select(bin, chrom, beg, end, size)
#
fo = sprintf("%s/%s/01.rds", dirw, gt)
saveRDS(tw, fo)
to = tw %>% mutate(beg=beg-1) %>% select(chrom,beg,end,bin)
fo = sprintf("%s/%s/01.bed.gz", dirw, gt)
write_tsv(to, fo, col_names=F)
}
#}}}

#{{{ read lifted bins & SNPs
for (i in c(3)) {
#tgt='B'; qry='W'
tgt = comps$tgt[i]; qry = comps$qry[i]
f_map = glue("{dirw}/{tgt}to{qry}/11.map.bed")
f_snp = glue("{dirw}/{tgt}to{qry}/12.snp.bed")
tm = sum_liftover(f_map, f_snp)
fo = glue("{dirw}/{tgt}to{qry}/20.rds")
saveRDS(tm, fo)
fo = glue("{dirw}/{tgt}to{qry}/20.tsv.gz")
write_tsv(tm, fo, na='')
#
x = tm %>% mutate(n_block=ifelse(!is.na(ins) & !is.na(del) & (n_ins+n_del>4 | ins+del>100), 0, n_block))
x %>% count(n_block)
to = x %>% filter(n_block>0) %>% select(chrom2,beg2,end2,bin) %>%
    arrange(chrom2, beg2)
fo = glue("{dirw}/{tgt}to{qry}/21.bed.gz")
write_tsv(to, fo, col_names=F, na='')
}
#}}}

#{{{ assign gene features to bins
for (i in 1:nrow(gts)) {
    gt = gts$gt[i]
    fb = sprintf("%s/%s/01.bed.gz", dirw, gt)
    fg = sprintf("%s/%s/05.gene.bed", dirw, gt)
    fo = sprintf("%s/%s/31.gene.tsv.gz", dirw, gt)
    to = sum_gene_itv(fb, fg)
    write_tsv(to, fo)
}

for (i in 1:nrow(comps)) {
    tgt = comps$tgt[i]; qry = comps$qry[i]
    fb = sprintf("%s/%sto%s/20.bed.gz", dirw, tgt, qry)
    fg = sprintf("%s/%s/05.gene.bed", dirw, qry)
    fo = sprintf("%s/%sto%s/31.gene.tsv.gz", dirw, tgt, qry)
    to = sum_gene_itv(fb, fg)
    write_tsv(to, fo)
}
#}}}

