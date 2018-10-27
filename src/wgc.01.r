source('wgc.fun.r')
source('Location.R')
require(GenomicRanges)

qry = 'W22'
tgt = 'B73'
diri = sprintf("%s/raw_output/%s_%s", dird, qry, tgt)
#{{{ read genome config
if(qry %in% search()) detach(qry, character.only = T)
if(tgt %in% search()) detach(tgt, character.only = T)
qcfg = attach(file.path(dirg, qry, '55.rda'), name = qry)
tcfg = attach(file.path(dirg, tgt, '55.rda'), name = tgt)
qcfg$cmap = qcfg$bed.chrom %>% 
    mutate(cs = c(0,cumsum(end-start+1)[-length(chrom)])) %>%
    mutate(gstart = start + cs, gend = end + cs) %>%
    select(chrom, gstart, gend)
tcfg$cmap = tcfg$bed.chrom %>% 
    mutate(cs = c(0,cumsum(end-start+1)[-length(chrom)])) %>%
    mutate(gstart = start + cs, gend = end + cs) %>%
    select(chrom, gstart, gend)
qcfg$gr_chrom = with(qcfg$bed.chrom, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
tcfg$gr_chrom = with(tcfg$bed.chrom, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
qcfg$gr_gap = with(qcfg$bed.gap, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
tcfg$gr_gap = with(tcfg$bed.gap, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
qcfg$gene = qcfg$loc.gene %>% group_by(gid) %>%
    summarise(chrom = chrom[1], start = min(start), end = max(end),
              ttype = ttype[1], size = end - start + 1)
tcfg$gene = tcfg$loc.gene %>% group_by(gid) %>%
    summarise(chrom = chrom[1], start = min(start), end = max(end),
              ttype = ttype[1], size = end - start + 1)
qcfg$gr_gene = with(qcfg$gene, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
tcfg$gr_gene = with(tcfg$gene, GRanges(seqnames = chrom, ranges = IRanges(start, end)))
#}}}

#{{{ #raw chain QC stats 
fi = file.path(diri, '05.bed')
ti = read_tsv(fi, col_names = c('tchrom','tstart','tend','srd','qchrom','qstart','qend','cid','mm'))
ti = ti %>% mutate(alnlen = tend - tstart)
tc = ti %>% group_by(cid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend), 
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend), 
              srd=srd[1], 
              alnlen = sum(alnlen), block = n(), mm = sum(mm)) %>%
    ungroup()

tc %>% mutate(pm = mm/size*100) %>%
    mutate(size_bin = cut(alnlen, 
                          breaks = c(0, 1e4, 5e4, 1e5, 1e6, 1e7, Inf),
                          labels = c('<10Kb', '10-100Kb', '100-500Kb',
                                     '0.5-1Mb', '1-10Mb', '10Mb+'))) %>%
    group_by(size_bin) %>%
    summarise(nchain = n(), alnlen = sum(alnlen),
              block = sum(block), mm = sum(mm),
              pm.mean = mm/size*100,
              pm.q25 = quantile(pm, .25),
              pm.q50 = quantile(pm, .5),
              pm.q75 = quantile(pm, .75)) %>% ungroup()
#}}}

fi = file.path(diri, '10.bed')
ti = read_tsv(fi, col_names = c('tchrom','tstart','tend','srd','qchrom','qstart','qend','cid','mm'))
tc = ti %>% mutate(alnlen = tend - tstart) %>%
    group_by(cid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend), 
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend), 
              srd=srd[1], alnlen = sum(alnlen), block = n(), mm = sum(mm))

#{{{ aln stats & genic coverge stats
t.gr_aln = with(ti, GRanges(seqnames=tchrom, ranges=IRanges(tstart+1, tend))) 
q.gr_aln = with(ti, GRanges(seqnames=qchrom, ranges=IRanges(qstart+1, qend))) 
t.gr_chn = with(tc, GRanges(seqnames=tchrom, ranges=IRanges(tstart+1, tend))) 
q.gr_chn = with(tc, GRanges(seqnames=qchrom, ranges=IRanges(qstart+1, qend)))
t.gr_gene_aln = intersect(tcfg$gr_gene, t.gr_aln)
q.gr_gene_aln = intersect(qcfg$gr_gene, q.gr_aln)
#
x = intersect_basepair(tcfg$gr_gene, t.gr_aln)
gt = tcfg$gene %>% mutate(osize = x) %>%
    mutate(po = osize/size)
gt %>%
    mutate(po.bin = cut(po, breaks = c(0,.5,.75,.95,1,Inf), right = F)) %>%
    count(po.bin) %>%
    mutate(pct = n/sum(n))
#
x = intersect_basepair(qcfg$gr_gene, q.gr_aln)
gq = qcfg$gene %>% mutate(osize = x) %>%
    mutate(po = osize/size) 
gq %>%
    mutate(po.bin = cut(po, breaks = c(0,.5,.75,.95,1,Inf), right = F)) %>%
    count(po.bin) %>%
    mutate(pct = n/sum(n))
#
t.genome = sum(width(tcfg$gr_chrom))
t.gap = sum(width(tcfg$gr_gap))
t.nogap = t.genome - t.gap
t.gene = sum(width(tcfg$gr_gene))
t.aln = sum(width(t.gr_aln))
t.gene.aln = sum(width(t.gr_gene_aln)) 
q.genome = sum(width(qcfg$gr_chrom))
q.gap = sum(width(qcfg$gr_gap))
q.nogap = q.genome - q.gap
q.gene = sum(width(qcfg$gr_gene))
q.aln = sum(width(q.gr_aln))
q.gene.aln = sum(width(q.gr_gene_aln))
t.ngene = nrow(gt)
q.ngene = nrow(gq)
t.ngene.syn = sum(gt$po >= .75)
q.ngene.syn = sum(gq$po >= .75)
#
tt = tibble(name=character(), qnum=numeric(), qpct=numeric(), tnum=numeric(), tpct=numeric()) %>%
    add_row(name = 'genome space', qnum = q.genome, qpct = 1,
            tnum = t.genome, tpct = 1) %>%
    add_row(name = 'gap-free space', qnum = q.nogap, qpct = q.nogap/q.genome,
            tnum = t.nogap, tpct = t.nogap/t.genome) %>%
    add_row(name = 'genic space', qnum = q.gene, qpct = q.gene/q.genome,
            tnum = t.gene, tpct = t.gene/t.genome) %>%
    add_row(name = 'genic space in synteny', 
            qnum = q.gene.aln, qpct = q.gene.aln/q.gene,
            tnum = t.gene.aln, tpct = t.gene.aln/t.gene) %>%
    add_row(name = 'total genes',
            qnum = q.ngene, qpct = 1,
            tnum = t.ngene, tpct = 1) %>%
    add_row(name = 'genes  w. >=75% covered in synteny', 
            qnum = q.ngene.syn, qpct = q.ngene.syn/q.ngene,
            tnum = t.ngene.syn, tpct = t.ngene.syn/t.ngene)
tt %>% mutate(qpct = sprintf("%4.1f%%", qpct*100),
              tpct = sprintf("%4.1f%%", tpct*100))
#}}}

#{{{ variant stats & gene effect
ft = sprintf("%s/15.%s.tsv", diri, tgt)
fq = sprintf("%s/15.%s.tsv", diri, qry)
et = read_tsv(ft, col_names = 
              c('chrom','pos','refl','altl','type','eff','impact','gid','tid'))
eq = read_tsv(fq, col_names = 
              c('chrom','pos','refl','altl','type','eff','impact','gid','tid'))

#{{{ variant num and size distri.
brks = c(0,10,100,1e3,5e3,1e4,Inf)
labs = c('0-10bp','10-100bp','100bp-1kb','1-5kb','5-10kb','10kb+')
nsnp = sum(et$type == 'snp')
ei = et %>% filter(refl==1, altl>1) %>%
    mutate(len = altl - refl) %>%
    mutate(len = cut(len, breaks = brks, labels = labs, right = T)) %>%
    group_by(len) %>%
    summarise(num = n(), base = sum(as.numeric(len))) %>% ungroup() %>%
    mutate(name = as.character(len),
           pct = num/sum(num), bpct = base/sum(base)) %>% 
    select(name, num, pct, base, bpct)
eis = ei %>% group_by(1) %>% 
    summarise(name = 'Insertions', 
              num = sum(num), pct = NA, base = sum(base), bpct = NA) %>%
    ungroup() %>% select(-`1`)
ed = et %>% filter(refl>1, altl==1) %>%
    mutate(len = refl - altl) %>%
    mutate(len = cut(len, breaks = brks, labels = labs, right = T)) %>%
    group_by(len) %>%
    summarise(num = n(), base = sum(as.numeric(len))) %>% ungroup() %>%
    mutate(name = as.character(len),
           pct = num/sum(num), bpct = base/sum(base)) %>% 
    select(name, num, pct, base, bpct)
eds = ed %>% group_by(1) %>% 
    summarise(name = 'Deletions', 
              num = sum(num), pct = NA, base = sum(base), bpct = NA) %>%
    ungroup() %>% select(-`1`)
em = et %>% filter(refl > 1, altl > 1) %>%
    group_by(1) %>%
    summarise(name = 'Mixed', num = n(), pct = NA, base = sum(refl), bpct = NA) %>%
    ungroup() %>% select(-`1`)
#
tv = tibble(name=character(), num=numeric(), pct=numeric(), base=numeric(), bpct=numeric()) %>%
    add_row(name = 'SNP', num = nsnp, pct = NA, base = nsnp, bpct = NA) %>%
    bind_rows(eis) %>%
    bind_rows(ei) %>%
    bind_rows(eds) %>%
    bind_rows(ed) %>%
    bind_rows(em)
#}}}

sum_vnteff <- function(ef, gsyn, syns = c('syntenic', 'non-syntenic'),
              impacts = c('high','moderate','low','modifier','no change')) {
    #{{{
    ef2 = ef %>% filter(eff != 'intergenic_region') %>%
        mutate(gid = str_split_fixed(gid, '-', n=2)[,1],
               eff = str_split_fixed(eff, '&', n=2)[,1])
    ef3 = ef2 %>%
        select(gid, tid, impact, eff) %>%
        mutate(impact = factor(tolower(impact), levels = impacts)) %>%
        group_by(gid) %>%
        arrange(impact) %>%
        filter(row_number() == 1) %>% ungroup()
    ef4 = gsyn %>% select(gid,ttype,po) %>%
        mutate(syn = ifelse(po>=.75, 'syntenic', 'non-syntenic')) %>%
        left_join(ef3, by = 'gid') %>%
        count(syn, impact, eff) %>%
        mutate(impact = as.character(impact)) %>%
        mutate(impact = ifelse(is.na(impact), 'no change', impact)) %>%
        mutate(syn = factor(syn, levels = syns)) %>%
        mutate(impact = factor(impact, levels = impacts)) %>%
        arrange(syn, impact, eff)
    ef4
    #}}}
}
#
et0 = sum_vnteff(et, gt) 
et1 = et0 %>% group_by(syn) %>%
    summarise(tnum = sum(n)) %>% ungroup() %>% mutate(tpct = tnum/sum(tnum))
et2 = et0 %>% group_by(syn, impact) %>%
    summarise(tnum = sum(n)) %>% ungroup() %>% mutate(tpct = tnum/sum(tnum))
et3 = et0 %>% group_by(syn, impact, eff) %>%
    summarise(tnum = sum(n)) %>% ungroup() %>% mutate(tpct = tnum/sum(tnum))
eq0 = sum_vnteff(eq, gq) 
eq1 = eq0 %>% group_by(syn) %>%
    summarise(qnum = sum(n)) %>% ungroup() %>% mutate(qpct = qnum/sum(qnum))
eq2 = eq0 %>% group_by(syn, impact) %>%
    summarise(qnum = sum(n)) %>% ungroup() %>% mutate(qpct = qnum/sum(qnum))
eq3 = eq0 %>% group_by(syn, impact, eff) %>%
    summarise(qnum = sum(n)) %>% ungroup() %>% mutate(qpct = qnum/sum(qnum))
#
ef1 = et1 %>% left_join(eq1, by = 'syn') 
ef2 = et2 %>% left_join(eq2, by = c('syn','impact'))
ef3 = et3 %>% left_join(eq3, by = c('syn','impact','eff'))
#}}}

fo = sprintf("%s/05_stats/%s_%s.rda", dird, qry, tgt)
save(tt, tv, ef1, ef2, ef3, file = fo) 

#{{{ whole genome pairwise alignment dotplot
max_dist = 5000
ta = ti %>% arrange(cid, tchrom, tstart) 
ta1 = ta[-nrow(ta),] %>% 
    transmute(cid1 = cid, srd = srd, ts1 = tstart, te1 = tend, qs1 = qstart, qe1 = qend)
ta2 = ta[-1,] %>% 
    transmute(cid2 = cid, ts2 = tstart, te2 = tend, qs2 = qstart, qe2 = qend)
ta3 = ta1 %>% bind_cols(ta2) %>%
    mutate(tdist = ts2 - te1,
           qdist = ifelse(srd == '-', qs1 - qe2, qs2 - qe1)) %>%
    mutate(flag = ifelse(cid1 == cid2,
                  ifelse(tdist+qdist >= max_dist, 1, 0), 1)) %>%
    mutate(bcid = cumsum(flag))
ta = ta %>% mutate(bcid = c(1,1+ta3$bcid))
tc = ta %>% mutate(alnlen = tend - tstart) %>% 
    group_by(bcid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend), srd=srd[1],
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend), 
              alnlen = sum(alnlen))
tc
describe(tc$alnlen)
#
tp = tc %>% group_by(bcid) %>%
    inner_join(tcfg$cmap, by = c('tchrom'='chrom')) %>%
    mutate(tstart = tstart+gstart-1, tend = tend+gstart-1) %>%
    select(-gstart, -gend) %>%
    inner_join(qcfg$cmap, by = c('qchrom'='chrom')) %>%
    mutate(qstart = qstart+gstart-1, qend = qend+gstart-1) %>%
    select(-gstart, -gend)
tpx = tcfg$cmap %>% mutate(gpos = (gstart+gend)/2)
tpy = qcfg$cmap %>% mutate(gpos = (gstart+gend)/2)
#
p = ggplot(tp) +
    geom_segment(aes(x = tstart, xend = tend, y = qstart, yend = qend), size = .3) +
    geom_hline(data = qcfg$cmap[-1], aes(yintercept = gstart), color='gray', alpha = .7, size = .2) +
    geom_vline(data = tcfg$cmap[-1], aes(xintercept = gstart), color='gray', alpha = .7, size = .2) +
    scale_x_continuous(breaks = tpx$gpos, labels = tpx$chrom, expand = c(0,0)) +
    scale_y_continuous(breaks = tpy$gpos, labels = tpy$chrom, expand = c(0,0)) +
    #scale_color_d3() +
    otheme(xtext = T, ytext = T, xticks = T, yticks = T)
fp = sprintf("%s/07_dotplot/%s_%s.pdf", dird, qry, tgt)
ggsave(p, filename = fp, width = 8, height = 8)
#}}}

#{{{ translocation
tc %>% filter(alnlen < 10000) %>% select(-srd)

cid = 58
idx = which(tc$cid == cid)
tchrom = tc$tchrom[idx]; tbeg = tc$tstart[idx]; tend = tc$tend[idx]
qchrom = tc$qchrom[idx]; qbeg = tc$qstart[idx]; qend = tc$qend[idx]
#tbeg = tbeg - 5000; tend = tend + 5000
#qbeg = qbeg - 5000; qend = qend + 5000
sprintf("%s:%d-%d  %s:%d-%d\n", tchrom, tbeg, tend, qchrom, qbeg, qend)

ta1 = ti %>% filter(tchrom == !!tchrom, tend <!!tbeg) %>%
    arrange(desc(tend)) %>% print(n=3, width=Inf)
ta2 = ti %>% filter(tchrom == !!tchrom, tstart >!!tend) %>%
    arrange(tstart) %>% print(n=3, width=Inf)

#}}}

#{{{ compile wgc stats table
qrys = c('Mo17', "W22", "PH207", "PHB47")
tgt = 'B73'
tz = tibble(qry = qrys, tgt = tgt) %>%
    mutate(s = map2(qry, tgt, .f <- function(x,y) attach(sprintf("%s/05_stats/%s_%s.rda", dird, x, y))))

tz2 = tz %>%
    mutate(tt = map(s, 'tt'), tv = map(s, 'tv'),
           ef1 = map(s, 'ef1'), ef2 = map(s, 'ef2'), ef3 = map(s, 'ef3'))

tt = tz2 %>% select(qry, tgt, tt) %>% unnest() %>%
    mutate(comp = sprintf("%s vs %s", qry, tgt)) %>% 
    select(-qry, -tgt) %>% select(comp, everything())
tv = tz2 %>% select(qry, tgt, tv) %>% unnest() %>%
    mutate(comp = sprintf("%s vs %s", qry, tgt)) %>% 
    select(-qry, -tgt) %>% select(comp, everything())
ef1 = tz2 %>% select(qry, tgt, ef1) %>% unnest() %>%
    mutate(comp = sprintf("%s vs %s", qry, tgt)) %>% 
    select(-qry, -tgt) %>% select(comp, everything())
ef2 = tz2 %>% select(qry, tgt, ef2) %>% unnest() %>%
    mutate(comp = sprintf("%s vs %s", qry, tgt)) %>% 
    select(-qry, -tgt) %>% select(comp, everything())
ef3 = tz2 %>% select(qry, tgt, ef3) %>% unnest() %>%
    mutate(comp = sprintf("%s vs %s", qry, tgt)) %>% 
    select(-qry, -tgt) %>% select(comp, everything())

fo = sprintf("%s/05_stats/all.rda", dird)
save(tt, tv, ef1, ef2, ef3, file = fo)
#}}}


