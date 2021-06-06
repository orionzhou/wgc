require(tidyverse)
require(glue)
require(ggsci)
source("/datalus/weiyu/git/rmaize/R/genome.R")
source("/datalus/weiyu/git/rmaize/R/plot.R")

tgt = 'Taestivum_D'
qry = 'Atauschii_AS60'
qry = 'Atauschii_AY61'
gconf1 = read_genome_conf(qry)
gconf2 = read_genome_conf(tgt)

#{{{ make a pairwise dotplot
fl1 = sprintf("%s/50_annotation/15.bed", genome_dir(qry))
fl2 = sprintf("%s/50_annotation/15.bed", genome_dir(tgt))
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
#}}}
#
dirw = '~/projects/wgc/data/31_wheat_synteny'
#{{{ dotplot
tz1 = tl1 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
tz2 = tl2 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
#
diri = sprintf('/datalus/weiyu/projects/wgc/data/raw/%s-%s', qry, tgt)
fi = sprintf("%s/xref.pairs", diri)
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
fp = glue("{dirw}/05.dotplot.pdf")
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

#{{{ make a one-to-one synteny plot
size1 = gconf1$chrom %>% select(chrom, size)
size2 = gconf2$chrom %>% select(chrom, size)
tz1 = flattern_gcoord_prepare(size1, gap=5e7)
tz2 = flattern_gcoord_prepare(size2, gap=5e7)
#
tl1 = gconf1$gene.loc %>% filter(ttype=='mRNA',etype=='exon') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup() %>%
    inner_join(tz1 %>% select(chrom,offset), by='chrom') %>%
    mutate(cpos = pos + offset) %>% select(-offset)
tl2 = gconf2$gene.loc %>% filter(ttype=='mRNA',etype=='exon') %>%
    arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
    summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup() %>%
    inner_join(tz2 %>% select(chrom,offset), by='chrom') %>%
    mutate(cpos = pos + offset) %>% select(-offset)
#}}}
#
#{{{ synteny block plot
ty = tibble(y = 1:2, gt = c(qry, tgt))
diri = sprintf('/datalus/weiyu/projects/wgc/data/raw/%s-%s', qry, tgt)
fi = sprintf("%s/xref.pairs", diri)
ti = read_tsv(fi) %>%
    inner_join(tl1, by=c('gid1'='tid')) %>% select(-pos,-gid1) %>%
    rename(gid1=gid, chrom1=chrom, pos1=cpos) %>%
    inner_join(tl2, by=c('gid2'='tid')) %>% select(-pos,-gid2) %>%
    rename(gid2=gid, chrom2=chrom, pos2=cpos) %>%
    mutate(qry=!!qry, tgt=!!tgt) %>%
    arrange(qry, tgt, bid, pos2)
#
ti2 = ti %>% group_by(qry,tgt,bid, chrom1, chrom2) %>%
    summarise(qBeg=pos1[1], qEnd=pos1[n()], tBeg=pos2[1], tEnd=pos2[n()]) %>%
    ungroup() %>% mutate(i = 1:n(), srd = ifelse(qBeg < qEnd, '+', '-')) %>%
    inner_join(ty,by=c('qry'='gt')) %>% rename(y1=y) %>% 
    inner_join(ty,by=c('tgt'='gt')) %>% rename(y2=y) 
#
tpy = tz1 %>% mutate(gt=qry) %>% bind_rows(tz2 %>% mutate(gt=tgt)) %>%
    select(gt, chrom, start, end, pos) %>%
    inner_join(ty, by='gt') %>% mutate(lab=str_replace(chrom,'chr',''))
tp1 = ti2 %>% select(i, srd, y=y2, tBeg, tEnd) %>%
    gather(type, pos, -i, -y, -srd)
tp2 = ti2 %>% select(i, srd, y=y1, qBeg, qEnd) %>%
    gather(type, pos, -i, -y, -srd)
coordmap = c("tBeg"=1,'tEnd'=2,'qEnd'=3,'qBeg'=4)
tp = tp1 %>% rbind(tp2) %>% mutate(i2 = coordmap[type]) %>% arrange(i, i2) %>%
    mutate(y = ifelse(i2<=2, y-.1, y+.1))
#
fp = glue("{dirw}/06.synteny.pdf")
wd=8; ht=3
p_syn = ggplot(tp) +
    geom_polygon(aes(x=pos,y=y,group=i, fill=srd), alpha=.4,
                 size=0,color=NA) +
    geom_rect(data=tpy, aes(xmin=start,xmax=end,ymin=y-.1,ymax=y+.1),
              fill=NA, color='black', size=.3) +
    geom_text(data=tpy, aes(x=pos,y=y,label=lab), size=2.5) +
    scale_x_continuous(expand=expansion(mult=c(.001, .001))) +
    scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.01))) +
    scale_fill_manual(name = 'synteny block', values=c("royalblue",'red')) +
    otheme(legend.pos='none', ytext=T,ytick=F, panel.border=F,
        margin=c(.2,.2,.2,.2))
ggsave(p_syn, filename=fp, width=wd, height=ht)
#}}}
