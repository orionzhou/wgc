source('functions.R')
dirw = glue("{dird}/05_maize_merged")

#{{{ prepare NAM vcf to merge
fi = "~/projects/genome/nf/genomes.xlsx"
ti = read_xlsx(fi) %>% filter(str_detect(genome, '^Zmays')) %>%
    filter(!alias %in% c("B73","PH207",'B73v5','BxM','B104'))

ti2 = ti %>% mutate(fv = glue("{dird}/raw/{genome}-Zmays_B73/t.vcf.gz")) %>%
    mutate(fv = normalizePath(fv)) %>%
    mutate(fve = map_lgl(fv, file.exists))
ti2 %>% count(fve)

fo = glue("{dirw}/01.file.list.txt")
write(ti2$fv, fo)
#}}}



