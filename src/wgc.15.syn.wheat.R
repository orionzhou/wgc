source("functions.R")
dirw = glue('{dird}/31_wheat_synteny')
fg = glue("{dirw}/wheat.xlsx")
get_comps  <- function(fg) {
    #{{{
    ti = read_xlsx(fg, col_names='gt')
    gts = ti$gt; ngt = length(gts)
    tibble(qry=gts[1:(ngt-1)], tgt=gts[2:ngt])
    #}}}
}
tc = get_comps(fg)
prepare_tz <- function(gconf) {
    #{{{
    size1 = gconf$chrom %>% select(chrom, size)
    tz1 = flattern_gcoord_prepare(size1, gap=5e7)
    tz1
    #}}}
}
prepare_tl <- function(gconf, opt='coord', gap=5e7) {
    #{{{
    tl = gconf$gene.loc %>%
        filter(ttype=='mRNA',etype=='exon') %>%
        arrange(gid, chrom, start, end) %>% group_by(gid,tid) %>%
        summarise(chrom=chrom[1], pos=(start[1]+end[n()])/2) %>% ungroup()
    if (opt == 'rank') {
        tl %>%
            arrange(chrom,pos) %>% mutate(idx=1:n()) %>%
            group_by(chrom) %>% mutate(cidx=1:n()) %>% ungroup()
    } else {
        size1 = gconf$chrom %>% select(chrom, size)
        tz = flattern_gcoord_prepare(size1, gap=gap)
        tl %>% inner_join(tz %>% select(chrom,offset), by='chrom') %>%
        mutate(cpos = pos + offset) %>% select(-offset)
    }
    #}}}
}
tg = read_xlsx(fg, col_names='gt') %>%
    mutate(gconf = map(gt, read_genome_conf)) %>%
    mutate(tz = map(gconf, prepare_tz)) %>%
    mutate(tl = map(gconf, prepare_tl, opt='coord', gap=5e7))

#{{{ write comparison pairs
fo = glue("{dirw}/comps.csv")
write_csv(tc, fo, col_names=F)
#}}}

#{{{ multi-pair syn plot
read_syn <- function(qry, tgt, tl1, tl2) {
    #{{{
    fi = glue('{dird}/raw/{qry}-{tgt}/xref.pairs')
    ti = read_tsv(fi) %>%
        inner_join(tl1, by=c('gid1'='tid')) %>% select(-pos,-gid1) %>%
        rename(gid1=gid, chrom1=chrom, pos1=cpos) %>%
        inner_join(tl2, by=c('gid2'='tid')) %>% select(-pos,-gid2) %>%
        rename(gid2=gid, chrom2=chrom, pos2=cpos) %>%
        #mutate(qry=!!qry, tgt=!!tgt) %>%
        arrange(bid, pos2)
    #
    ti2 = ti %>% group_by(bid, chrom1, chrom2) %>%
        summarise(qBeg=pos1[1], qEnd=pos1[n()], tBeg=pos2[1], tEnd=pos2[n()]) %>%
        ungroup() %>% mutate(i = 1:n(), srd = ifelse(qBeg < qEnd, '+', '-'))
    #}}}
}

tc1 = tc %>% inner_join(tg, by=c('qry'='gt')) %>%
    select(-gconf) %>% rename(tl1=tl, tz1=tz) %>%
    inner_join(tg, by=c('tgt'='gt')) %>% select(-gconf) %>% rename(tl2=tl, tz2=tz) %>%
    mutate(syn = pmap(list(qry,tgt,tl1,tl2),read_syn))

ty = tg %>% select(gt) %>% mutate(y = n():1)
tp0 = tc1 %>% select(qry,tgt,syn) %>% unnest(syn) %>%
    mutate(i = glue("{qry}-{tgt}-{i}")) %>%
    filter(!str_detect(chrom1, 'x'), !str_detect(chrom2, 'x')) %>%
    inner_join(ty,by=c('qry'='gt')) %>% rename(y1=y) %>% 
    inner_join(ty,by=c('tgt'='gt')) %>% rename(y2=y) 
tp1 = tp0 %>% select(i, srd, y=y2, tBeg, tEnd) %>%
    gather(type, pos, -i, -y, -srd)
tp2 = tp0 %>% select(i, srd, y=y1, qBeg, qEnd) %>%
    gather(type, pos, -i, -y, -srd)
coordmap = c("tBeg"=1,'tEnd'=2,'qEnd'=3,'qBeg'=4)
ht_chr=.05
tp = tp1 %>% bind_rows(tp2) %>% mutate(i2 = coordmap[type]) %>% arrange(i, i2) %>%
    mutate(y = ifelse(i2<=2, y+ht_chr, y-ht_chr))
#
tpx = tg %>% select(gt, tz) %>% unnest(tz) %>%
    select(gt, chrom, start, end, pos) %>%
    filter(!str_detect(chrom, 'x')) %>%
    inner_join(ty, by='gt') %>% mutate(lab=str_replace(chrom,'chr',''))

# y-axis color
cols5 = c('black', pal_aaas()(5)[c(1,4,3,2)])
ty2 = ty %>% rename(label=gt) %>% mutate(gt = 'bg') %>%
    mutate(gt = ifelse(str_detect(label, '_A'), 'A', gt)) %>%
    mutate(gt = ifelse(str_detect(label, '_B'), 'B', gt)) %>%
    mutate(gt = ifelse(str_detect(label, '_D'), 'D', gt)) %>%
    mutate(gt = ifelse(str_detect(label, 'Atauschii_'), 'D2', gt)) %>%
    mutate(gt = factor(gt, levels=c('bg','A','B','D','D2'))) %>%
    mutate(col = cols5[as.numeric(gt)])

fp = glue("{dirw}/06.synteny.pdf")
wd=10; ht=8
p = ggplot(tp) +
    geom_polygon(aes(x=pos,y=y,group=i, fill=srd), alpha=.4,
                 size=0,color=NA) +
    geom_rect(data=tpx, aes(xmin=start,xmax=end,ymin=y-ht_chr,ymax=y+ht_chr),
              fill=NA, color='grey', size=.3) +
    geom_text(data=tpx, aes(x=pos,y=y,label=lab), size=2, color='black') +
    scale_x_continuous(expand=expansion(mult=c(.001, .001))) +
    scale_y_continuous(breaks=ty$y, labels=ty$gt, expand=expansion(mult=c(.01,.01))) +
    scale_fill_manual(name = 'synteny block', values=c("royalblue",'red')) +
    otheme(legend.pos='none', ytext=T,ytick=F, panel.border=F,
        margin=c(.2,.2,.2,.2)) +
    theme(axis.text.y = element_text(color=ty2$col, size=7))
ggsave(p, filename=fp, width=wd, height=ht)
#}}}


#{{{ # single pair dot plot [obsolete]
i=6
tgt = tc$tgt[i]
qry = tc$qry[i]
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
#{{{ dotplot
tz1 = tl1 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
tz2 = tl2 %>% group_by(chrom) %>%
    summarise(start=min(idx), end=max(idx), pos=(start+end)/2) %>% ungroup()
#
diri = glue('{dird}/raw/{qry}-{tgt}')
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

#}}}


