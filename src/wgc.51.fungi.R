source("functions.R")
diri = glue("{dird}/Laburnicola")
dirw = glue('{dird}/25_Laburnicola')

gts = c('R11','R22','R44','R19','R8')
gs = c('R11','R22','R44','R19','R11xR22','R8','R11xR44')
tg = tibble(g=gs) %>% mutate(g = factor(g, levels=gs)) %>%
    mutate(conf = map(g, read_genome_conf, subdir="Laburnicola"))
t_gs = tg %>% mutate(t_gs = map(conf, "gene")) %>% unnest(t_gs) %>% select(-conf)
cmps = tibble(
  qry=c('R11','R22','R22','R11xR22','R11xR44'),
  tgt=c('R44','R44','R11','R19',    'R8')
) %>% mutate(cmp = glue("{qry} - {tgt}")) %>% mutate(cmp = as_factor(cmp))

#{{{ plot ancestral karyotype / subgenome 
tg = tibble(g=gts) %>% mutate(fa = glue("{diri}/{g}.anc")) %>%
    mutate(ta = map(fa, read_tsv, col_names=c("chrom","start","end","anc",'sub'))) %>%
    select(-fa) %>%
    unnest(ta) %>% mutate(span = end-start+1) %>%
    mutate(hap = ifelse(g %in% c("R8","R19"), glue("{g}_{LETTERS[sub]}"), g))
thap = tg %>% distinct(g, hap) %>% mutate(g = factor(g, levels=gs)) %>%
    arrange(g, hap) %>% mutate(hap = as_factor(hap))
haps = thap %>% pull(hap)
tg = tg %>% mutate(hap = factor(hap, levels=haps))
tgs = tg %>% group_by(g,chrom,sub) %>% summarize(span=sum(span)) %>% ungroup() %>%
    arrange(g,chrom,desc(span)) %>%
    group_by(g,chrom) %>% slice(1) %>% ungroup() %>% select(g,chrom,sub)
tgl = tg %>% group_by(g,chrom) %>% summarise(size=max(end)) %>% ungroup()
tga = tg %>% filter(g == 'R11') %>% arrange(chrom) %>% mutate(achrom=1:n()) %>%
    select(anc, achrom)
tg2 = tg %>% inner_join(tgs, by=c("g","chrom","sub")) %>% inner_join(tga, by='anc')

tg3 = tg2 %>% group_by(hap,chrom,achrom) %>%
    summarise(span = sum(span)) %>% ungroup()
tg3a = tg3 %>% arrange(hap,achrom, desc(span)) %>%
    group_by(hap,achrom) %>% slice(1) %>% ungroup()
tg3b = tg3 %>% arrange(hap,chrom, desc(span)) %>%
    group_by(hap,chrom) %>% slice(1) %>% ungroup()
tg3a1 = tg3a %>% count(hap, chrom) %>% filter(n==1) %>% select(-n) %>%
    inner_join(tg3a, by=c('hap','chrom'))
tg3a2 = tg3a %>% count(hap, chrom) %>% filter(n>1) %>% select(-n) %>%
    inner_join(tg3b, by=c('hap','chrom'))
tg3a3 = tg3b %>% 
    left_join(tg3a %>% select(hap,chrom,achrom2=achrom), by=c('hap','chrom')) %>%
    filter(is.na(achrom2)) %>% select(-achrom2)
tg3 = rbind(tg3a1, tg3a2, tg3a3)
tg3 %>% select(hap, chrom, achrom) %>%
    group_by(hap,achrom) %>% summarise(chrom=str_c(chrom, collapse=',')) %>%
    spread(hap,chrom) %>% print(n=25)

ty = tg3 %>% inner_join(thap, by='hap') %>%
    arrange(hap,achrom,desc(span)) %>% group_by(hap,achrom,g) %>%
    mutate(i=1:n()) %>% ungroup()
ty2 = ty %>% distinct(achrom,i) %>% arrange(achrom, i) %>% mutate(y = 1:n())
ty = ty %>% inner_join(ty2, by=c('achrom','i')) %>% select(hap,achrom,i,y,chrom,g)

cols21 = c(pal_igv()(18), pal_aaas()(10)[c(10,3,5)])
cols21 = c(pal_ucscgb()(18), pal_igv()(18)[c(8,15,17)])
tys = ty %>% filter(i==1)
tp = tg %>% select(-hap) %>% inner_join(tga, by='anc') %>%
    inner_join(ty %>% select(-achrom), by=c('g','chrom')) %>%
    mutate(achrom = as.factor(achrom))
tgl2 = tgl %>% inner_join(ty, by=c("g","chrom")) %>%
    mutate(chromLabel = as.integer(str_replace(chrom, "^chr", '')))
p = ggplot() +
    geom_rect(data=tgl2, aes(xmin=1,xmax=size,ymin=y-.05,ymax=y+.05),fill='black',size=0) +
    geom_rect(data=tp,aes(xmin=start,xmax=end,ymin=y-.2,ymax=y+.2,fill=achrom), size=0) +
    geom_text(data=tgl2,aes(x=10,y=y-.25,label=chromLabel), size=2.5,hjust=0,vjust=0) +
    facet_wrap(hap~., nrow=1) +
    scale_y_reverse(breaks=tys$y, labels=tys$achrom, expand=c(.05,.01)) +
    scale_fill_manual(values=cols21) +
    otheme(xtitle=F, ytitle=F, xtext=F, ytext=T, xtick=F, ytick=F,
         legend.pos='none', strip.compact=F)
fo = glue("{dirw}/08.karyotype.pdf")
ggsave(p, filename=fo, width=12, height=6)
#}}}

#{{{ collect block/gene pair and write
rcfg2tg <- function(rcfg) {
    #{{{
    rcfg$gene %>% arrange(chrom, start, end) %>%
        mutate(idx = 1:n()) %>%
        mutate(pos = (start+end)/2) %>%
        group_by(chrom) %>% mutate(cidx = 1:n()) %>% ungroup() %>%
        select(gid, tid, chrom, pos, srd, idx, cidx)
    #}}}
}
split_str <- function(block1, block2, ks, sep="_") {
    #{{{
    ks = str_replace(ks, "^_", "") 
    tibble(
           order1 = as.integer(str_split(block1, sep)[[1]]),
           order2 = as.integer(str_split(block2, sep)[[1]]),
           ks = as.numeric(str_split(ks, sep)[[1]])
           )
    #}}}
}
ti = cmps %>%
    mutate(qcfg = glue("{diri}/{qry}-{tgt}/wgdi/00.q.rds")) %>%
    mutate(tcfg = glue("{diri}/{qry}-{tgt}/wgdi/00.t.rds")) %>%
    mutate(fb= glue("{diri}/{qry}-{tgt}/wgdi/14.block.csv")) %>%
    mutate(fk= glue("{diri}/{qry}-{tgt}/wgdi/13.ks.tsv")) %>%
    mutate(qcfg = map(qcfg, readRDS)) %>%
    mutate(tcfg = map(tcfg, readRDS)) %>%
    mutate(rb = map(fb, read_csv, col_names=T)) %>%
    mutate(rk = map(fk, read_tsv, col_names=T)) %>%
    select(-fb, -fk)
wgdi_post <- function(tb0, tk0, qcfg, tcfg, min_len=15, min_homo1=.75) {
    #{{{
    tg1 = rcfg2tg(qcfg)
    tg2 = rcfg2tg(tcfg)
    tb = tb0 %>% rename(bid = id, chrom1=chr1, chrom2=chr2)
    tk = tk0 %>% rename(tid1=id1, tid2=id2)
    block = tb %>% select(-block1, -block2, -ks) %>%
        filter(length >= min_len) %>% mutate(primary=homo1>=min_homo1)
    tb2 = tb %>% select(bid, block1, block2, ks) %>%
        mutate(data = pmap(list(block1, block2, ks), split_str, sep='_')) %>%
        select(bid, data) %>% unnest(data)
    tb3 = tb2 %>% rename(cidx1=order1,cidx2=order2) %>%
        inner_join(block %>% select(bid, chrom1, chrom2, primary), by='bid') %>%
        inner_join(tg1 %>% select(gid1=gid,tid1=tid,chrom1=chrom,idx1=idx,cidx1=cidx), by=c('chrom1','cidx1')) %>%
        inner_join(tg2 %>% select(gid2=gid,tid2=tid,chrom2=chrom,idx2=idx,cidx2=cidx), by=c('chrom2','cidx2')) %>%
        arrange(bid, idx1) %>%
        left_join(tk, by=c("tid1",'tid2')) %>%
        select(bid,chrom1,chrom2,cidx1,cidx2,gid1,gid2,tid1,tid2,primary,
            ka.NG86=ka_NG86, ks.NG86=ks_NG86,
            ka.YN00=ka_YN00, ks.YN00=ks_YN00)
    list(block=block, pair=tb3)
    #}}}
}
to = ti %>%
    mutate(r = pmap(list(rb, rk, qcfg, tcfg), wgdi_post)) %>%
    mutate(block = map(r, 'block'), pair=map(r, 'pair')) %>%
    select(-r)
fo = glue("{dirw}/01.rds")
saveRDS(to, fo)

ti1 = to %>% select(cmp, cmp, block) %>% unnest(block)
ti2 = to %>% select(cmp, pair) %>% unnest(pair)
ti1 %>% count(cmp)
ti2 %>% count(cmp)

to1 = ti1 %>% rename(comparison=cmp,block=bid)
to2 = ti2 %>% rename(comparison=cmp,block=bid,
    gene1=gid1,gene2=gid2, rna1=tid1,rna2=tid2, order1=cidx1, order2=cidx2)

fo1 = glue("{dirw}/03.block.tsv")
fo2 = glue("{dirw}/03.gene_pair.tsv")
write_tsv(to1, fo1)
write_tsv(to2, fo2)
#}}}

#{{{ make xref table and plot composition
fi = glue("{dirw}/01.rds")
ti = readRDS(fi)
ti1 = ti %>% select(qry, tgt, block) %>% unnest(block) %>%
    filter(length >= 10, homo1 >= .5)
ti2 = ti %>% select(qry, tgt, pair) %>% unnest(pair) %>%
    select(qry,tgt,bid,gid1,gid2,ks,order1=cidx1,order2=cidx2) %>%
    inner_join(ti1 %>% select(qry,tgt,bid), by=c('qry','tgt','bid'))
ti1 %>% count(qry, tgt)
ti2 %>% count(qry, tgt)

tg1 = t_gs %>% filter(g=='R11') %>% select(gid1=gid, size1=size.exon)
tg2 = t_gs %>% filter(g=='R44') %>% select(gid2=gid, size2=size.exon)
tp = ti2 %>% filter(qry=='R11', tgt=='R44') %>% distinct(gid1, gid2)
make_xref <- function(tp, gs1, gs2) {
    #{{{
    tp = tp %>% distinct(gid1, gid2)
    tg1 = gs1 %>% select(gid1=gid, size1=size.exon)
    tg2 = gs2 %>% select(gid2=gid, size2=size.exon)
    gids1 = tp %>% count(gid1) %>% filter(n>1) %>% pull(gid1)
    gids2 = tp %>% count(gid2) %>% filter(n>1) %>% pull(gid2)
    tp1 = tp %>% filter(gid1 %in% gids1) %>% mutate(type1='1-m')
    tp2 = tp %>% filter(gid2 %in% gids2) %>% mutate(type2='1-m')
    tp3 = tp %>%
        left_join(tp1, by=c('gid1','gid2')) %>%
        left_join(tp2, by=c('gid1','gid2')) %>%
        replace_na(list(type1='1-1',type2='1-1')) %>%
        mutate(type=ifelse(type1=='1-m', ifelse(type2=='1-m', 'm-m', '1-m'),
                           ifelse(type2=='1-m', 'm-1', '1-1'))) %>%
        inner_join(tg1, by='gid1') %>% 
        inner_join(tg2, by='gid2') %>%
        select(gid1, gid2, type, size1, size2)
    tp3 %>% count(type)
    tp3a = tg1 %>% filter(!gid1 %in% tp3$gid1) %>% mutate(type='1-0')
    tp3b = tg2 %>% filter(!gid2 %in% tp3$gid2) %>% mutate(type='0-1')
    tp4 = tp3 %>% bind_rows(tp3a) %>% bind_rows(tp3b)
    tp4
    #}}}
}

tg0 = t_gs %>% group_by(g) %>% nest() %>% ungroup() %>% rename(gs=data)
tc = ti2 %>% group_by(qry,tgt) %>% nest() %>% ungroup() %>%
    inner_join(tg0, by=c('qry'='g')) %>% rename(gs1=gs) %>%
    inner_join(tg0, by=c('tgt'='g')) %>% rename(gs2=gs) %>%
    mutate(xref = pmap(list(data, gs1, gs2), make_xref)) %>%
    select(qry, tgt, xref)

fo = glue("{dirw}/10.xref.rds")
saveRDS(tc, fo)

types = c('1-1','1-m m-1 m-m','1-0','0-1')
tc2 = tc %>% unnest(xref) %>%
    mutate(size.diff=abs(size1-size2)) %>%
    mutate(type=ifelse(type %in% c('1-m','m-1','m-m'), '1-m m-1 m-m', type)) %>%
    mutate(type = factor(type, levels=types)) %>%
    mutate(size.diff.bin=cut(size.diff, c(0,10,50,100,Inf),
        labels=c("<10",'10,50','50,100','>100'), include.lowest=T))
tc2s = tc2 %>%
    distinct(type, size.diff.bin) %>%
    mutate(lab = ifelse(type=='1-1', glue("{type} ({size.diff.bin})"), as.character(type))) %>%
    arrange(type, size.diff.bin) %>%
    mutate(lab = as_factor(lab))
tp = tc2 %>% inner_join(cmps, by=c('qry','tgt')) %>%
    inner_join(tc2s, by=c('type','size.diff.bin')) %>%
    count(cmp,lab) %>%
    rename(tag1=cmp, tag2=lab)
cols7 = c(brewer.pal(4,'Pastel2'), pal_aaas()(3))
cols7 = c(brewer.pal(4,'Reds')[4:1], pal_simpsons()(1), brewer.pal(3,'Paired')[2:1])

p = cmp_cnt1(tp, ytext=T, ypos='right', legend.title='', fills=cols7,
             legend.pos='top.left', legend.dir='v') +
    o_margin(.1,.3,.1,.3)
fo = glue("{dirw}/15.gene.xref.pdf")
ggsave(p, file=fo, width=6, height=6)
p = cmp_prop1(tp, ytext=T, ypos='right', legend.title='', fills=cols7,
             legend.pos='none', legend.dir='v') +
    o_margin(.1,.3,.1,.3)
fo = glue("{dirw}/15.gene.xref.p.pdf")
ggsave(p, file=fo, width=6, height=6)
#}}}



