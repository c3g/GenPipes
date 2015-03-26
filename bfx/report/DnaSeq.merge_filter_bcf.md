### Single Nucleotide Variations

#### SNV Calling

Single Nucleotide Variants (SNPs and INDELs) are called using [SAMtools and BCFtools][SAMtools]\ [@samtools]. SAMtools (mpileup) collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format. It does not call variants. BCFtools applies the prior and performs the actual calling. Every time a mapped read shows a mismatch from the reference genome, it applies Bayesian statistics to figure out whether the mismatch is caused by a real SNP. It incorporates different types of information, such as the number of different reads that share a mismatch from the reference, the sequence quality data, and the expected sequencing error rates, and it essentially figures out whether it is more likely that the observed mismatch is due simply to a sequencing error, or because of a true SNV.

The SNV calls are reprocessed by BCFtools (varfilter), which does an additional filtering of the variants and transforms the output into the VCF format. The final vcf files are filtered for long 'N' INDELs which are sometimes introduced and cause excessive memory usage in downstream analyses.

#### SNV Merging

All samples from the cohort are merged together into one global VCF file.

Table: Species-specific Global VCF Columns Description ([download full table](HumanVCFformatDescriptor.tsv))

