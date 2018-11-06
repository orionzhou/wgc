Whole genome comparison of maize assemblies
================
November 05, 2018





































  - Pipeline components:
      - [minimap2](https://github.com/lh3/minimap2) was used to align
        one assembly to another;
      - [blat chain/net
        tools](http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto)
        were used to process alignment resuls and build synteny
        chains/nets;
      - [bcftools](https://samtools.github.io/bcftools/) and
        [GATK4](https://software.broadinstitute.org/gatk/gatk4) were
        used to call variants;
      - [snpEff](http://snpeff.sourceforge.net) was used to evaluate
        variant effects on sytenic genes
  - whole genome synteny plots
      - [Mo17 vs B73](/data/07_dotplot/Mo17_B73.pdf)
      - [W22 vs B73](/data/07_dotplot/W22_B73.pdf)
      - [PH207 vs B73](/data/07_dotplot/PH207_B73.pdf)
      - [PHB47 vs B73](/data/07_dotplot/PHB47_B73.pdf)
  - Table 1. Whole genome alignment statistics.

![](t1.png)<!-- -->

  - Table 2. Summary of variants called by synteny comparison.

![](t2.png)<!-- -->

  - Table 3. Summary of variant effects on syntenic genes.

![](t3.png)<!-- -->
