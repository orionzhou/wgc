source('functions.R')
dirw = file.path(dird, '92_bin_summary')
diri = file.path(dird, '91_genome_bins')
gts = tibble(gt0 = c('B73','Mo17','W22','PH207'), gt=c('B','M','W','P')) %>%
    mutate(gt0 = str_c("Zmays", gt0, sep='_'))
comps = tibble(tgt = c(rep("B", 3)), qry = c('M','W','P'))

#{{{ demo & mapping stats
tgt='B'; qry='P'
fi = sprintf("%s/%sto%s/20.rds", diri, tgt, qry)
tm = readRDS(fi)
tm1 = tm %>% filter(n_block>0) %>% select(-size2) %>%
    replace_na(list(aln=100,n_del=0,n_ins=0,del=0,ins=0)) %>%
    select(-ins) %>% print(n=40,width=Inf)
#
types = c('unmap','mapped w. no indel no snp','mapped w. no indel but >=1 snp',
           'mapped w. >=1 indel')
x = tm %>% mutate(type = ifelse(n_block==0, types[1],
    ifelse(n_block==1, ifelse(n_snp==0,types[2],types[3]), types[4])))
x %>% count(type) %>% mutate(type=factor(type, levels=types)) %>%
    arrange(type) %>% mutate(pct = percent(n/nrow(x), accuracy=1))
#}}}

#{{{ plot gene feature composition
#{{{ read in
ti1 = gts %>% mutate(fi = sprintf("%s/%s/31.gene.tsv.gz", diri, gt)) %>%
    select(ctag = gt, fi) %>%
    mutate(res = map(fi, read_tsv)) %>%
    select(ctag, res)
ti2 = comps %>% mutate(fi = sprintf("%s/%sto%s/31.gene.tsv.gz", diri, tgt, qry)) %>%
    mutate(ctag = sprintf("%s_to_%s", tgt, qry, qry)) %>%
    select(ctag, qry, tgt, fi) %>%
    mutate(res = map(fi, read_tsv)) %>%
    select(ctag, qry, tgt, res)
ti2a = ti2 %>% mutate(ctag = str_c(ctag, qry, sep=' | ')) %>% select(ctag, res)
subset_bins <- function(t1, t2) t1 %>% filter(bin %in% t2$bin)
join_bins <- function(t1, t2) t1 %>% rename(type1=type) %>%
    inner_join(t2, by='bin') %>% rename(type2=type)
ti2b = ti2 %>% select(ctag, tgt, res2=res) %>%
    inner_join(ti1, by=c("tgt"='ctag')) %>% rename(res1=res) %>%
    mutate(res = map2(res1, res2, subset_bins)) %>%
    mutate(ctag = str_c(ctag, tgt, sep=" | ")) %>% select(ctag, res)
tc = ti2 %>% select(ctag, qry, tgt, res2=res) %>%
    inner_join(ti1, by=c("tgt"='ctag')) %>% rename(res1=res) %>%
    mutate(res = map2(res1, res2, join_bins)) %>%
    select(ctag, qry, tgt, res)
#}}}

#{{{ gene feature composition
tp0 = ti1 %>% bind_rows(ti2a) %>% bind_rows(ti2b) %>%
    mutate(x = map(res, get_cnt <- function(r) r %>% count(type))) %>%
    select(ctag, x) %>% unnest(x)

tp = tp0 %>% mutate(type = ifelse(str_detect(type, "^(mRNA)|(intergenic)"), type, 'other_mRNA')) %>%
    group_by(ctag, type) %>% summarise(n=sum(n)) %>% ungroup()
ctags = tp %>% distinct(ctag) %>% arrange(nchar(ctag), ctag) %>% pull(ctag)
types = tp %>% distinct(type) %>% arrange(type) %>% pull(type)
tp = tp %>% mutate(ctag=factor(ctag, levels=ctags)) %>%
    mutate(type=factor(type, levels=types)) %>%
    select(tag1=ctag, tag2=type, n)
#
p = cmp_proportion1(tp, xtitle='', ytitle='', xangle=30, margin.l=.1,
    legend.pos='right', legend.dir='v', legend.title='', barwidth=.8,
    oneline=T, expand.x=c(.01,.01), expand.y=c(.01,.04), pal='simpsons')
fo = file.path(dirw, '05.gene.feature.pdf')
ggsave(p, file=fo,width=10, height=8)
#}}}

#{{{ comparisons
for (i in 1:nrow(tc)) {
ctag = tc$ctag[i]; qry=tc$qry[i]; tgt=tc$tgt[i]
tp = tc$res[[i]] %>% count(type1,type2) %>%
    mutate(type1 = ifelse(str_detect(type1, "^(mRNA)|(intergenic)"), type1, 'other_mRNA')) %>%
    mutate(type2 = ifelse(str_detect(type2, "^(mRNA)|(intergenic)"), type2, 'other_mRNA')) %>%
    group_by(type1,type2) %>% summarise(n=sum(n)) %>% ungroup()
types1 = tp %>% distinct(type1) %>% arrange(type1) %>% pull(type1)
types2 = tp %>% distinct(type2) %>% arrange(type2) %>% pull(type2)
tp = tp %>%
    mutate(type1=factor(type1, levels=types1)) %>%
    mutate(type2=factor(type2, levels=types2)) %>%
    select(tag1=type1,tag2=type2,n)
#
xtit = sprintf("%s features", tgt)
ytit = sprintf("%s features", qry)
p = cmp_proportion1(tp, xtitle=xtit, ytitle='', xangle=0, margin.l=.1,
    legend.pos='right', legend.dir='v', legend.title=ytit, barwidth=.8,
    oneline=T, expand.x=c(.05,.05), expand.y=c(.01,.04), pal='simpsons')
fo = sprintf("%s/06.gene.%s.pdf", dirw, ctag)
ggsave(p, file=fo,width=6, height=6)
}
#}}}
#}}}
