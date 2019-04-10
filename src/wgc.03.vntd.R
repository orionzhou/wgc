source('functions.R')

#{{{ window-based synteny bases & variant density
qrys = c('Mo17', "W22", "PH207", "PHB47")
tgt = 'B73'
comps = sprintf("%s_%s", qrys, tgt)

fw = '~/projects/polyTE/data/nucdiv/01.win.bed'
tw = read_tsv(fw, col_names = c('chrom','start','end'))

get_syn_vnt <- function(qry, tgt, tw) {
    #{{{
    grw = with(tw, GRanges(seqnames=chrom, ranges=IRanges(start+1, end)))
    diri = sprintf("%s/raw_output/%s_%s", dird, qry, tgt)
    fc = file.path(diri, '10.bed')
    tc = read_tsv(fc, col_names = c('tchrom','tstart','tend','srd',
                                    'qchrom','qstart','qend','cid','nsnp'))
    grc = with(tc, GRanges(seqnames=tchrom, ranges=IRanges(tstart+1, tend)))
    fv = file.path(diri, '10.vnt.bed')
    tv = read_tsv(fv, col_names = c('tchrom','tstart','tend','srd',
                                    'qchrom','qstart','qend','cid','vtype','ref','alt'))
    tv1 = tv %>% filter(vtype == 'snp') %>% select(-vtype)
    tv2 = tv %>% filter(vtype != 'snp') %>% select(-vtype)
    grv1 = with(tv1, GRanges(seqnames=tchrom, ranges=IRanges(tstart+1, tend)))
    grv2 = with(tv2, GRanges(seqnames=tchrom, ranges=IRanges(tstart+1, tend)))
    bp_syn = intersect_basepair(grw, grc)
    num_snp = intersect_count(grw, grv1)
    num_indel = intersect_count(grw, grv2)
    tp = tw %>% mutate(bp_syn = bp_syn, num_snp = num_snp, num_indel = num_indel)
    tp
    #}}}
}

tp = tibble(qry = qrys, tgt = tgt) %>%
    mutate(data = map2(qry, tgt, get_syn_vnt, tw = tw)) %>%
    unnest()
fo = file.path(dird, '11_vnt_density', '01.rda')
save(tp, file = fo)
#}}}

#{{{ synteny / variant density plot
dirw = file.path(dird, '11_vnt_density')
fi = file.path(dirw, '01.rda')
x = load(fi)
range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
vtypes = c('psyn','dsnp','didl')
vtype.labs = c('syntenic content','SNP density','InDel density')
tp = tp %>% 
    mutate(posi = (end + start) / 2000000,
           psyn = range01(bp_syn/(end-start), na.rm = T),
           dsnp = range01(num_snp/bp_syn, na.rm = T),
           didl = range01(num_indel/bp_syn, na.rm = T)) %>%
    select(qry, chrom, posi, psyn, dsnp, didl) %>%
    gather(vtype, val, -chrom, -posi, -qry) %>%
    mutate(vtype = factor(vtype, levels = vtypes))

qrys = c('Mo17', "W22", "PH207", "PHB47")
qry = qrys[1]
tp1 = tp %>% filter(qry == !!qry)
p = ggplot(tp1) +
    geom_line(aes(x = posi, y = val, color = vtype), size = .2, alpha = .8) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_npg(labels = vtype.labs) +
    facet_grid(chrom~.) +
    otheme(legend.pos = 'top.center.out', legend.dir = 'h',
           xgrid = T, ygrid = T, 
           xtext = T, xticks = T)
fp = sprintf("%s/05.%s.pdf", dirw, qry)
ggsave(p, filename = fp, width = 10, height = 8)
#}}}

#{{{ pete's 100-bp tiles synteny & variant stats
qrys = c('Mo17', "PH207")
tgt = 'B73'
comps = sprintf("%s_%s", qrys, tgt)

comp = comps[2]
dirw = '~/projects/wgc/data/91_pete_tile'
f1 = sprintf("%s/15.%s.syn.bed", dirw, comp)
f2 = sprintf("%s/16.%s.snp.bed", dirw, comp)
f3 = sprintf("%s/17.%s.idl.bed", dirw, comp)
t1 = read_tsv(f1, col_names=c('chrom','start','end','id','score','phase',
                              'c2','s2','e2','srd','c3','s3','e3','cid','nsnp',
                              'bp'))
t2 = read_tsv(f2, col_names=c('chrom','start','end','id','score','phase','n'))
t3 = read_tsv(f3, col_names=c('chrom','start','end','id','score','phase','n'))
#
tp = t1 %>% group_by(chrom,start,end,id) %>%
    summarise(bp_syn=sum(bp)) %>% ungroup() %>%
    mutate(size = end-start) %>%
    select(chrom,start,end,id, size, bp_syn) %>%
    left_join(t2, by=c('chrom','start','end','id')) %>%
    select(-score,-phase) %>% rename(n_snp=n) %>%
    left_join(t3, by=c('chrom','start','end','id')) %>%
    select(-score,-phase) %>% rename(n_indel=n)

fo = sprintf("%s/20.%s.bed", dirw, comp)
write_tsv(tp, fo, col_names=F)
# liftOver -bedPlus=3 old chain new unmap
#}}}

