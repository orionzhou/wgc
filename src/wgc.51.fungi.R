source("functions.R")
dirw = glue('{dird}/51_fungi')
diri = glue("{dirr}/Laburnicola")

gs=c('R11','R22','R44','R19','R8')
cmps = tibble(
  qry=c('R11','R22','R22','R11xR22','R11xR44'),
  tgt=c('R44','R44','R11','R19',    'R8')
)

#{{{ plot ancestral karyotype / subgenome 
tg = tibble(g=gs) %>% mutate(fa = glue("{diri}/{g}.anc")) %>%
    mutate(ta = map(fa, read_tsv, col_names=c("chrom","start","end","anc",'sub'))) %>%
    select(-fa) %>%
    unnest(ta) %>% mutate(span = end-start+1) %>%
    mutate(hap = ifelse(g %in% c("R8","R19"), glue("{g}_{LETTERS[sub]}"), g))
thap = tg %>% distinct(g, hap) %>% mutate(g = factor(g, levels=gs)) %>%
    arrange(g, hap)
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

ti = cmps %>%
    mutate(fi= glue("{diri}/{qry}-{tgt}/wgdi/20.rds")) %>%
    mutate(r = map(fi, readRDS))

ti1 = ti %>% mutate(block=map(r, 'block')) %>%
    select(qry, tgt, block) %>% unnest(block) %>%
    filter(length >= 10, homo1 >= .5)
ti2 = ti %>% mutate(pair=map(r, 'pair')) %>%
    select(qry, tgt, pair) %>% unnest(pair) %>%
    select(qry,tgt,bid,gid1,gid2,ks,order1=cidx1,order2=cidx2) %>%
    inner_join(ti1 %>% select(qry,tgt,bid), by=c('qry','tgt','bid'))
ti1 %>% count(qry, tgt)
ti2 %>% count(qry, tgt)


to1 = ti1 %>%
    rename(query=qry,target=tgt,block=bid)
to2 = ti2 %>%
    rename(query=qry,target=tgt,block=bid,gene1=gid1,gene2=gid2)

fo1 = glue("{dirw}/03.block.tsv")
fo2 = glue("{dirw}/03.gene_pair.tsv")
write_tsv(to1, fo1)
write_tsv(to2, fo2)
