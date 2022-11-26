source('functions.R')
qry='Alyrata'; tgt='Athaliana_Col0'
subdir = 'common'
diri = glue('~/projects/s3/zhoup-wgc/raw/{subdir}/{qry}-{tgt}')
fi = glue("{diri}/xref.pairs")
fp = glue("{diri}/05.dotplot.pdf")

#{{{
gconf1 = read_genome_conf(qry)
gconf2 = read_genome_conf(tgt)
#
size1 = gconf1$chrom %>% select(chrom,size=end)
size2 = gconf2$chrom %>% select(chrom,size=end)
tz1 = flattern_gcoord_prepare(size1, gap=0)
tz2 = flattern_gcoord_prepare(size2, gap=0)
#
#fl1 = sprintf("%s/50_annotation/15.bed", genome_dir(qry))
#fl2 = sprintf("%s/50_annotation/15.bed", genome_dir(tgt))
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

#{{{ dotplot
ti = read_tsv(fi) %>%
    inner_join(tl1, by=c('gid1'='tid')) %>% select(-pos) %>%
    rename(tid1=gid1, gid1=gid, chrom1=chrom, idx1=idx) %>%
    inner_join(tl2, by=c('gid2'='tid')) %>% select(-pos) %>%
    rename(tid2=gid2, gid2=gid, chrom2=chrom, idx2=idx) %>%
    arrange(bid, idx1)

bids18 = ti %>% count(bid) %>% arrange(desc(n)) %>%
    filter(row_number()<=18) %>% pull(bid)
tp = ti %>% mutate(bid = ifelse(bid %in% bids18, bid, 'other')) %>%
    mutate(bid = factor(bid, levels=c(bids18,'other')))
cols19 = c(pal_igv()(18), 'black')
wd=6; ht=6
p1 = ggplot(tp) +
  #geom_segment(aes(x=idx1s,xend=idx1e, y=idx2s,yend=idx2e,color=ftype), size=.6) +
  geom_point(aes(x=idx2,y=idx1, color=bid), size=.1) +
  geom_vline(xintercept = tz2$end, alpha=.1) +
  geom_hline(yintercept = tz1$end, alpha=.1) +
  scale_x_continuous(name=tgt, breaks=tz2$pos, labels=tz2$chrom, expand=c(0,0)) +
  scale_y_continuous(name=qry, breaks=tz1$pos, labels=tz1$chrom, expand=c(0,0)) +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(name = 'synteny block', values=cols19) +
  otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
         legend.pos='none', legend.dir='v', legend.title=T)
ggsave(p1, filename=fp, width=wd, height=ht)
#}}}

