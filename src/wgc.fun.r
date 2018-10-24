#{{{ load required libraries, define common variables
require(grid)
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(Hmisc)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
source('~/projects/genomes/src/ge.fun.r')
#
dirp = '~/projects/wgc'
dird = file.path(dirp, 'data')
dirc = '/scratch.global/zhoux379/wgc'
#}}}


