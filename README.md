# whole genome comparison code and results

This repository hosts codes and results for pairwise whole genome comparison between any two genome (moslty maize) assemblies. 

- [minimap2](https://github.com/lh3/minimap2) was used to align one assembly to another;
- [blat chain/net tools](http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto) were used to process alignment resuls and build synteny chains/nets;
- [bcftools](https://samtools.github.io/bcftools/) and [GATK4](https://software.broadinstitute.org/gatk/gatk4) were used to call variants;
- [snpEff](http://snpeff.sourceforge.net) was used to evaluate variant effects on sytenic genes

See [this page](/Rmd/wgc.md) for a sumamry of comparing four maize de novo assemblies ([W22](), [Mo17](), [PH207]() and [PHB47]()) to the [B73 reference](https://www.maizegdb.org/assembly).
