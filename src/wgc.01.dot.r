source('wgc.fun.r')
require(GenomicRanges)

qry = 'W22'
tgt = 'B73'
#{{{ read genome config
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
qcfg$gr_exon = reduce(with(qcfg$loc.exon, GRanges(seqnames = chrom, ranges = IRanges(start, end))))
tcfg$gr_exon = reduce(with(tcfg$loc.exon, GRanges(seqnames = chrom, ranges = IRanges(start, end))))
#}}}

#{{{ QC filter and save to 10.bed
diri = sprintf("%s/raw_output/%s_%s", dird, qry, tgt)
fi = file.path(diri, '01.bed')
ti = read_tsv(fi, col_names = F)
colnames(ti) = c('tchrom','tstart','tend','srd','qchrom','qstart','qend','cid','mm')
ti = ti %>% mutate(size = tend - tstart)
tis = ti %>% group_by(cid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend), 
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend), 
              srd=srd[1], 
              size = sum(size), block = n(), mm = sum(mm)) %>%
    filter(size >= 50000)
ti = ti %>% filter(cid %in% tis$cid)

tis %>% mutate(pm = mm/size*100) %>%
    mutate(size_bin = cut(size, breaks = c(0, 1e4, 1e5, 1e6, 1e7, Inf))) %>%
    group_by(size_bin) %>%
    summarise(nchain = n(), size = sum(size),
              block = sum(block), mm = sum(mm),
              pm.q25 = quantile(pm, .25),
              pm.q50 = quantile(pm, .5),
              pm.q75 = quantile(pm, .75)) %>% ungroup()

# write filtered output, forward and reverse chain
fo = file.path(diri, '10.bed')
write_tsv(ti, fo, col_names = F)
#
tsize = file.path(dirg, tgt, '15_intervals/01.chrom.sizes')
qsize = file.path(dirg, qry, '15_intervals/01.chrom.sizes')
fo1 = sprintf("%s/10.%s_%s.chain", diri, tgt, qry)
cmd=sprintf("chain fromBed %s %s %s > %s", fo, tsize, qsize, fo1)
system(cmd)
fo2 = sprintf("%s/10.%s_%s.chain", diri, qry, tgt)
cmd=sprintf("chainSwap %s %s", fo1, fo2)
system(cmd)

# write filtered SNP bed
fn = file.path(diri, '01.snp.bed')
tn = read_tsv(fn, col_names = F) %>%
    filter(X8 %in% tis$cid)
tn1 = tn %>% select(X4, X5, X6, X7, X8) %>%
    arrange(X4, X5)
tn2 = tn %>% select(X1, X2, X3, X7, X8) %>%
    arrange(X1, X2)
fn0 = sprintf("%s/11.snp.bed", diri)
write_tsv(tn, fn0, col_names = F)
fn1 = sprintf("%s/11.snp.%s.bed", diri, qry)
write_tsv(tn1, fn1, col_names = F)
fn2 = sprintf("%s/11.snp.%s.bed", diri, tgt)
write_tsv(tn2, fn2, col_names = F)
#}}}

#{{{ aln stats & genic coverge stats
t.gr_aln = with(ti, GRanges(seqnames = tchrom, ranges = IRanges(tstart+1, tend))) 
t.gr_exon_aln = intersect(t.gr_aln, tcfg$gr_exon)
q.gr_aln = with(ti, GRanges(seqnames = qchrom, ranges = IRanges(qstart+1, qend))) 
q.gr_exon_aln = intersect(q.gr_aln, qcfg$gr_exon)
#
t.genome = sum(width(tcfg$gr_chrom))
t.gap = sum(width(tcfg$gr_gap))
t.nogap = t.genome - t.gap
t.exon = sum(width(tcfg$gr_exon))
t.aln = sum(width(t.gr_aln))
t.exon.aln = sum(width(t.gr_exon_aln)) 
q.genome = sum(width(qcfg$gr_chrom))
q.gap = sum(width(qcfg$gr_gap))
q.nogap = q.genome - q.gap
q.exon = sum(width(qcfg$gr_exon))
q.aln = sum(width(q.gr_aln))
q.exon.aln = sum(width(q.gr_exon_aln)) 
#
outputs = c(
sprintf("%s:\n", tgt),
sprintf("%10.0f: genome space\n", t.genome),
sprintf("%10d: (%4.1f%%) gap-free space\n", t.nogap, t.nogap/t.genome*100),
sprintf("%10d: (%4.1f%%) synteny space\n", t.aln, t.aln/t.genome*100),
sprintf("%10d: (%4.1f%%) exonic space\n", t.exon, t.exon/t.genome*100),
sprintf("%10d: (%4.1f%%) synteny exon space\n", t.exon.aln, t.exon.aln/t.exon*100),
sprintf("%s:\n", qry),
sprintf("%10.0f: genome space\n", q.genome),
sprintf("%10d: (%4.1f%%) gap-free space\n", q.nogap, q.nogap/q.genome*100),
sprintf("%10d: (%4.1f%%) synteny space\n", q.aln, q.aln/q.genome*100),
sprintf("%10d: (%4.1f%%) exonic space\n", q.exon, q.exon/q.genome*100),
sprintf("%10d: (%4.1f%%) synteny exon space\n", q.exon.aln, q.exon.aln/q.exon*100)
)
cat(outputs, sep = '')
#}}}


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
tc = ta %>% group_by(bcid) %>%
    summarise(tchrom=tchrom[1], tstart=min(tstart), tend=max(tend), srd=srd[1],
              qchrom=qchrom[1], qstart=min(qstart), qend=max(qend), 
              size = sum(size))
tc
describe(tc$size)

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
fp = sprintf("%s/%s_%s.pdf", dird, qry, tgt)
ggsave(p, filename = fp, width = 8, height = 8)
#}}}


