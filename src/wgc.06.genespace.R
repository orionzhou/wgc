source('functions.R')
require(GENESPACE)
dirw = glue("{dird}/06_genespace")

make_exampleDataDir(writeDir = dirw)

spec <- c("human","chimp","rhesus")
spec <- c("Osativa_MSU","Osativa_9311","Orufipogon",'NH001')
genomes = str_replace(spec, '_', '')
gpar <- init_genespace(
  genomeIDs=genomes, speciesIDs=spec, versionIDs=spec, ploidy=rep(1,length(spec)),
  diamondMode = "fast", orthofinderMethod = "fast", wd = dirw,
  orthofinderInBlk = F, overwrite = F, verbose = T,
  nCores = 4, minPepLen = 50,
  path2mcscanx = "/datalus/weiyu/software/MCScanX",
  gffString = "gff", pepString = "pep",
  path2orthofinder = "orthofinder", path2diamond = "diamond",
  rawGenomeDir = glue("{dirw}/demo2/rawGenomes"))

parse_annotations(gsParam = gpar,
  gffEntryType = "mRNA", gffIdColumn = "ID", gffStripText = "",
  headerEntryIndex = 1, headerSep = " ", headerStripText = "")

gpar = set_syntenyParams(gsParam = gpar)
gpar$params$synteny
gpar = run_orthofinder(gsParam=gpar, overwrite=T)
gpar = synteny(gsParam = gpar)

fo = glue("{dirw}/test.pdf")
pdf(fo, width=8, height=6)
ripSouceDat <- plot_riparianHits(gpar, 
  refGenome = genomes[1],
  #invertTheseChrs = data.frame(genome = "rhesus", chr = 2),
  genomeIDs=genomes[-3], labelTheseGenomes=genomes[-3],
  gapProp=.001, refChrCols = c("#BC4F43", "#F67243"),
  blackBg=F, returnSourceData=T, verbose=F)
dev.off()
