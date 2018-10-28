# whole genome comparison code and results

This repository hosts codes and results for pairwise whole genome comparison between any two genome (moslty maize) assemblies. 

- [minimap2]() was used to align one assembly to another;
- [blat chain/net tools]() were used to process alignment resuls and build synteny chains/nets;
- [bcftools]() and [GATK4]() were used to call variants;
- [snpEff]() was used to evaluate variant effects on sytenic genes

See [this page](/Rmd/wga.md) for a sumamry of comparing four maize de novo assemblies (W22, Mo17, PH207 and PHB47) to the [B73 reference]().
