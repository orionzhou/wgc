source('functions.R')
qry = 'Bdistachyon'; tgt = 'Osativa_MSU'
qry = 'Alyrata'; tgt = 'Athaliana_Col0'
subdir = 'common'
diri = glue('~/projects/s3/zhoup-wgc/raw/{subdir}/{qry}-{tgt}')
fi = glue("{diri}/xref.pairs")
dirw = glue('~/projects/s3/zhoup-wgc/21_{qry}-{tgt}')
#{{{ prepare data
gconf1 = read_genome_conf(qry, dirg=dirg)
gconf2 = read_genome_conf(tgt, dirg=dirg)
#
size1 = gconf1$chrom %>% select(chrom,size=end)
size2 = gconf2$chrom %>% select(chrom,size=end)
tz1 = flattern_gcoord_prepare(size1, gap=0)
tz2 = flattern_gcoord_prepare(size2, gap=0)
#
tl1 = gconf1$gene.loc %>% filter(ttype=='mRNA',etype=='exon') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup() %>%
    arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
tl2 = gconf2$gene.loc %>% filter(ttype=='mRNA',etype=='exon') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup() %>%
    arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
#
tz1 = tl1 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
tz2 = tl2 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
#}}}

fi = glue("{diri}/xref.t.tsv")
ti = read_tsv(fi, col_names=c('rice','brachy','type'))
ti %>% count(type)


