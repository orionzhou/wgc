require(tidyverse)
require(lubridate)
source("~/projects/maize.expression/src/me.fun.r")
dirw= '~/projects/wgc/Rmd'

dates = ymd('2018-11-3') + 0:6
ta = tibble(Date = dates, Freq = c(5,5,0,1,0,1,1)) %>%
    mutate(DateLabel = sprintf("%s %s", month(Date, label=T), mday(Date)))

p = ggplot(ta) +
    geom_bar(aes(x=Date, y=Freq, fill=as.character(Date)), stat='identity', width=.8) +
    scale_x_discrete(breaks = ta$Date, labels = ta$DateLabel) +
    scale_y_continuous(name = 'Frequency of Angry') +
    scale_fill_npg() +
    otheme(legend.pos = 'none', 
    xtext=T, xtitle=T, ytext=T, ytitle=T, xgrid=F, ygrid=T,
    xtick = T, ytick = T)
fo = file.path(dirw, 'july1.pdf')
ggsave(p, file = fo, width = 5, height = 4)

tl = tibble(Date = rep(dates, times=c(5,5,1,1,1,1,1)),
            level = c(2,3,2,1,2, 2,3,1,4,4, 0, 4, 0, 3, 1))
tls = tl %>% group_by(Date) %>% summarise(level = mean(level)) %>%
    mutate(DateLabel = sprintf("%s %s", month(Date, label=T), mday(Date)))
irs_note = c(
"Not angry",
'Feel OK but a little bit angry',
'Feel not OK and a little bit angry',
'Feel somewhat angry and irritate insights',
'Feel very angry but can control actions',
'Feel extremely angry and cannot control actions'
)
irs = tibble(level=0:5, note=irs_note) %>%
    mutate(lab = sprintf("%d - %s", level, note))
p = ggplot(tls) +
    geom_point(aes(x=Date, y=level)) +
    geom_line(aes(x=Date, y=level)) +
    geom_point(data = tibble(x=ymd('2018-11-3'),y=0,level=0:5), aes(x=x, y=y, color=level), alpha = 0) +
    scale_x_date(breaks = tls$Date, labels = tls$DateLabel) +
    scale_color_viridis(name = "Individual Rating Scales (IRS)", breaks = irs$level, labels = irs$lab) +
    scale_y_continuous(name = 'Level of Angry', limits = c(0,5)) +
    otheme(legend.pos = 'right',, legend.dir = 'v',
    xtext=T, xtitle=T, ytext=T, ytitle=T, xgrid=F, ygrid=T,
    xtick = T, ytick = T) +
    theme(legend.title.align = 0,
          legend.title = element_text(size = 9),
          legend.key.size = unit(1, 'lines'), 
          legend.text = element_text(size = 8))
fo = file.path(dirw, 'july2.pdf')
ggsave(p, file = fo, width = 9, height = 5)

