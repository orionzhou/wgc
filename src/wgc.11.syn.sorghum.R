#{{{
source('functions.R')
qry = 'Mo17'
qry = 'Sbicolor'
tgt = 'B73'
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
tl1 = gconf1$gene %>% filter(ttype=='mRNA') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    mutate(pos=(start+end)/2) %>%
    arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
tl2 = gconf2$gene %>% filter(ttype=='mRNA') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    mutate(pos=(start+end)/2) %>%
    arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
    group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
#
tz1 = tl1 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
tz2 = tl2 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
#}}}
diri = sprintf('%s/raw/%s_%s/20_synteny', dird, qry, tgt)
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
fi = file.path(diri, '06.q.blocks')
#tm = tj %>% inner_join(tl1, by=c('gid1'='gid')) %>% mutate(idx=cidx)
tm = read_tsv(fi,col_names=c('gid1','gid2a','gid2b')) %>%
    inner_join(tl1, by=c('gid1'='tid')) %>% mutate(idx=cidx)
fi = file.path(diri, '06.q.blocks.tracks')
tb = read_tsv(fi, col_names='bids') %>%
    mutate(ftype=sprintf('maize%d', 1:n())) %>%
    mutate(bid = str_split(bids, ',')) %>%
    select(-bids) %>% unnest() %>%
    mutate(bid = sprintf("b%02d", as.numeric(bid)+1))

tp = ti %>%
    #group_by(bid) %>%
    #summarise(idx1s = idx1[1], idx1e = idx1[n()],
    #          idx2s = idx2[1], idx2e = idx2[n()]) %>% ungroup() %>%
    left_join(tb, by='bid') %>% replace_na(list(ftype='unused'))
fp = sprintf("%s/01.pdf", dirw)
wd=8; ht=8
p1 = ggplot(tp) +
  #geom_segment(aes(x=idx1s,xend=idx1e, y=idx2s,yend=idx2e,color=ftype), size=.6) +
  geom_point(aes(x=idx1,y=idx2,color=ftype), size=.1) +
  geom_vline(xintercept = tz1$end, alpha=.1) +
  geom_hline(yintercept = tz2$end, alpha=.1) +
  scale_x_continuous(name=qry, breaks=tz1$pos, labels=tz1$chrom, expand=c(0,0)) +
  scale_y_continuous(name=tgt, breaks=tz2$pos, labels=tz2$chrom, expand=c(0,0)) +
  scale_color_brewer(palette = "Set1") +
  otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
         legend.pos='top.center.out', legend.dir='h')
ggsave(p1, filename=fp, width=wd, height=ht)
#}}}

#{{{ sorghum-maize fractionation plot
diff_frac <- function(chr, idxb, idxe, tm) {
    #{{{
    tm %>% filter(chrom == chr, idx >= idxb, idx <= idxe) %>%
        group_by(1) %>%
        summarise(n_sb=n(),
                  n_maize1=sum(gid2a!='.'), n_maize2=sum(gid2b!='.'),
                  n_both = sum(gid2a!='.' & gid2b!='.')) %>%
        ungroup() %>% select(-`1`)
    #}}}
}

#{{{ double check with James's pan-grass table
tl2s = tl2 %>% select(gid, tid)
tt = tm %>% select(chrom,gid,gid2a,gid2b) %>% rename(gid1=gid) %>%
    left_join(tl2s, by=c('gid2a'='tid')) %>%
    select(chrom,gid1,gid,gid2b) %>% rename(gid2a=gid) %>%
    left_join(tl2s, by=c('gid2b'='tid')) %>%
    select(chrom,gid1,gid2a,gid) %>% rename(gid2b=gid) %>%
    replace_na(list(gid2a='No Gene', gid2b='No Gene'))

fj = '~/projects/genome/data/B73/gene_mapping/sorghum3_intell_plusteff.csv'
tj = read_csv(fj) %>%
    select(gid1=sorghum3, gid2a=maize1_v4, gid2b=maize2_v4) %>%
    mutate(gid1 = str_to_upper(str_replace(gid1, 'Sobic\\.', 'SORBI_3'))) %>%
    mutate(gid2a = str_replace(gid2a, 'No Gene', '\\.')) %>%
    mutate(gid2b = str_replace(gid2b, 'No Gene', '\\.'))
tj2 = tj %>% rename(gid2a.j=gid2a, gid2b.j = gid2b)

tp = tt %>% inner_join(tj2, by='gid1') %>%
    mutate(acon=(gid2a!=gid2a.j & gid2a!='No Gene' & gid2a.j != 'No Gene'),
           bcon=(gid2b!=gid2b.j & gid2b!='No Gene' & gid2b.j != 'No Gene'))

tp %>% count(acon, bcon)
tp %>% filter(acon) %>% count(chrom)
tp %>% filter(bcon) %>% count(chrom)
tp %>% filter(acon & bcon) %>% count(chrom)


tp = tj %>%
    rename(subgenome1=gid2a, subgenome2=gid2b) %>%
    gather(ftype, gid2, -gid1) %>%
    filter(gid2 != '.') %>%
    inner_join(tl1, by=c('gid1'='gid')) %>% select(-pos,-tid) %>%
    rename(chrom1=chrom, idx1=idx) %>%
    inner_join(tl2, by=c('gid2'='gid')) %>% select(-pos,-tid) %>%
    rename(chrom2=chrom, idx2=idx)

fp = sprintf("%s/01.james.pdf", dirw)
wd=8; ht=8
p1 = ggplot(tp) +
  geom_point(aes(x=idx1,y=idx2,color=ftype), size=.2) +
  geom_vline(xintercept = tz1$end, alpha=.1) +
  geom_hline(yintercept = tz2$end, alpha=.1) +
  scale_x_continuous(name=qry, breaks=tz1$pos, labels=tz1$chrom, expand=c(0,0)) +
  scale_y_continuous(name=tgt, breaks=tz2$pos, labels=tz2$chrom, expand=c(0,0)) +
  scale_color_brewer(palette = "Set1") +
  otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
         legend.pos='top.center.out', legend.dir='h')
ggsave(p1, filename=fp, width=wd, height=ht)

#}}}

winsize=100; winstep=20
tw = tl1 %>% arrange(idx) %>% filter(chrom != 'Sb99') %>%
    group_by(chrom) %>% mutate(idx=1:n()) %>%
    summarise(ib=min(idx), ie=max(idx)) %>%
    ungroup() %>% mutate(data=map2(ib, ie, make_tile, winsize, winstep)) %>%
    select(chrom,data) %>% unnest() %>%
    mutate(data = pmap(list(chr=chrom, idxb=beg, idxe=end), diff_frac, tm)) %>%
    unnest()

ftypes = c('maize1','maize2','both')
tp = tw %>% mutate(pos = (beg+end)/2) %>%
    mutate(maize1 = n_maize1/n_sb, maize2 = n_maize2/n_sb,
           both = n_both/n_sb) %>%
    select(chrom,pos,maize1,maize2,both) %>%
    gather(ftype, prop, -chrom, -pos) %>%
    mutate(ftype = fct_relevel(ftype, ftypes))
#
p1 = ggplot(tp) +
  geom_line(aes(x=pos, y=prop, color=ftype, group=ftype), size=.5) +
  scale_x_continuous(name='Gene index', expand=c(0,0)) +
  scale_y_continuous(name='Proportion retained in maize', expand=expand_scale(mult=c(0,.05))) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~chrom, ncol=2, dir='v') +
  otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T, ygrid=T,
         legend.pos='top.center.out', legend.dir='h')
fp = file.path(dirw, '04.retained.maize.pdf')
ggsave(p1, filename=fp, width=8, height=8)
#}}}

#{{{ create maize subgenome gene ID list
tid2gid = tl2$gid; names(tid2gid) = tl2$tid
to = tm %>% select(sorghum=gid1, maize1=gid2a, maize2=gid2b) %>%
    mutate(maize1= ifelse(maize1 == '.', '.', tid2gid[maize1])) %>%
    mutate(maize2= ifelse(maize2 == '.', '.', tid2gid[maize2])) %>%
    mutate(ftype = ifelse(maize1=='.' & maize2=='.', 'non-syntenic', 'syntenic')) %>%
    mutate(ftype = ifelse(maize1!='.', 'maize1-retained', ftype)) %>%
    mutate(ftype = ifelse(maize2!='.', 'maize2-retained', ftype)) %>%
    mutate(ftype = ifelse(maize1!='.' & maize2!='.', 'both-retained', ftype))
to %>% count(ftype)

fo = file.path(dirw, '10.maize.pairs.tsv')
write_tsv(to, fo)

to = tm %>% select(gid1, subgenome1=gid2a, subgenome2=gid2b) %>%
    filter(subgenome1 != '.' | subgenome2 != '.') %>%
    mutate(ftype=ifelse(subgenome1=='.'|subgenome2=='.', 'fractionated', 'retained')) %>%
    gather(subgenome, gid, -gid1, -ftype) %>% filter(gid != '.') %>%
    distinct(gid, ftype, subgenome) %>%
    mutate(ftype=sprintf("%s_%s", subgenome, ftype)) %>%
    rename(tid=gid) %>% inner_join(tl2, by='tid') %>%
    select(gid, ftype)

fo = file.path(dirw, '10.maize.subgenome.tsv')
write_tsv(to, fo)
#}}}

ti2 = ti %>% group_by(bid) %>%
    summarise(chrom1 = chrom1[1], chrom2 = chrom2[1],
              start1 = min(idx1), end1 = max(idx1),
              start2 = min(idx2), end2 = max(idx2),
              srd = ifelse(idx2[n()]>idx2[1], '+', '-'),
              ng1 = length(unique(gid1)), ng2 = length(unique(gid2)),
              span1 = idx1[n()]-idx1[1]+1,
              span2 = ifelse(srd=='-', idx2[1]-idx2[n()]+1, idx2[n()]-idx2[1]+1)) %>%
    ungroup()

