#{{{
source('functions.R')
tgt = 'B73'
qry = 'W22'
gconf1 = read_genome_conf(qry)
gconf2 = read_genome_conf(tgt)
#
size1 = gconf1$chrom %>% select(chrom,size=end)
size2 = gconf2$chrom %>% select(chrom,size=end)
tz1 = flattern_gcoord_prepare(size1, gap=0)
tz2 = flattern_gcoord_prepare(size2, gap=0)
#
fl1 = sprintf("%s/50_annotation/15.bed", genome_dir(qry))
fl2 = sprintf("%s/50_annotation/15.bed", genome_dir(tgt))
tl1 = gconf1$loc.gene %>% filter(ttype=='mRNA',etype=='exon') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup() %>%
    arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
tl2 = gconf2$loc.gene %>% filter(ttype=='mRNA',etype=='exon') %>%
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
diri = sprintf('%s/raw_output/%s_%s/20_synteny', dird, qry, tgt)
dirw = sprintf('%s/21_%s_%s', dird, qry, tgt)
if (!file.exists(dirw)) dir.create(dirw)

#{{{ dotplot
fi = file.path(diri, '05.pairs')
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
fp = sprintf("%s/01.pdf", dirw)
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
         legend.pos='bottom.right', legend.dir='v', legend.title=T)
ggsave(p1, filename=fp, width=wd, height=ht)
#}}}

fs1 = file.path(diri, 'q.t.last.filtered')
ts1 = read_tsv(fs1, col_names = blast_cols12)
fs2 = file.path(diri, 'q.t.rbh')
ts2 = read_tsv(fs2, col_names=c('qseqid','sseqid','pident'))

fp = file.path(diri, '05.pairs')
tp = read_tsv(fp)
tp %>% count(gid1) %>% count(n)
tp %>% count(gid2) %>% count(n)
gid1r = ti$tid1[!tj$tid1 %in% tp$gid1]
gid2r = ti$tid2[!ti$tid2 %in% tp$gid2]

t_rbh = ts1 %>% filter(qseqid %in% gid1r, sseqid %in% gid2r)
length(unique(t_rbh$qseqid))
length(unique(t_rbh$sseqid))


blast_cols12=c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
fs = file.path(diri, 'q.t.last')
t_sim = read_tsv(fs, col_names=blast_cols12) %>%
    select(qseqid, sseqid, pident, evalue, score=bitscore)

fi = file.path(diri, '07.t.ortholog')
ti = read_tsv(fi, col_names=c('tid2','tid1')) %>%
    mutate(rbh=str_detect(tid1, "\\'$")) %>%
    mutate(tid1 = str_remove(tid1, "\\'$")) %>%
    mutate(type2 = ifelse(tid1 == '.', 'lost1', ifelse(rbh, 'rbh', 'synteny')))
ti %>% count(rbh,type2)

fj = file.path(diri, '07.q.ortholog')
tj = read_tsv(fj, col_names=c('tid1','tid2')) %>%
    mutate(rbh=str_detect(tid2, "\\'$")) %>%
    mutate(tid2 = str_remove(tid2, "\\'$")) %>%
    mutate(type1 = ifelse(tid2 == '.', 'lost2', ifelse(rbh, 'rbh', 'synteny')))
tj %>% count(rbh,type1)

tj2 = tj %>% select(tid1,tid2b=tid2,type1)
to = ti %>% full_join(tj2, by='tid1') %>%
    replace_na(list(tid2='.', type1='lost1', type2='lost2')) #%>%
    #left_join(t_sim, by=c('tid1'='qseqid','tid2'='sseqid'))

to %>% count(type1,type2)
#to %>% count(type1,type2,is.na(score))

ti2 = ti %>% select(tid2,tid1b=tid1,type2)
to = tj %>% full_join(ti2, by='tid2') %>%
    replace_na(list(tid1='.', type1='lost1', type2='lost2'))
to %>% count(type1,type2)
