### Variants Metrics

The variant metrics and statistics have been generated using the [SnpEff] software [@snpeff].

#### Summary SNV Metrics

General summary statistics are provided for the entire set of variants without considering the sample information.

Table: Summary Statistics for all Variants ([download full table](SNV.SummaryTable.tsv))

Metrics|Value
-----|-----
$snv_summary_table$

* Number of variants before filter: total number of variants in the cohort. Please note that when measuring the effect counts, each snp could be evaluated several times depending on the number of transcripts affected by the SNP
* Number of variants filtered out: number of filtered variants (MAPQ < 15)
* %: Number of variants filtered out / Number of variants before filter
* Number of not variants: number of variants where the reference allele is the same than the alternate allele
* %: Number of not variants / Number of variants before filter
* Number of variants processed: number of analysed after removing the filter and not variant ones
* Number of known variants: number of variants with a non-empty ID field
* %: Numbe of known variants / Number of variants processed
* Transitions: number of SNP interchanges of purines (A<->G) or of pyrimidines (C<->T)
* Transversions: number of SNP interchanges of purine for pyrimidine bases (A|G<->C|T)
* Ts Tv ratio: Transitions / Transversions
* missense: number of variants in coding region for which the variation generate a codon modification
* silent: number of variants in coding region for which the variation does not generate a codon modification
* missense silent ratio: missense / silent
* high impact: number of variants which an effect annotated in one of these categories: SPLICE_SITE_ACCEPTOR, SPLICE_SITE_DONOR, START_LOST, EXON_DELETED, FRAME_SHIFT, STOP_GAINED, STOP_LOST or RARE_AMINO_ACID
* low impact: number of variant which an effect annotated in one of these categories: SYNONYMOUS_START, NON_SYNONYMOUS_START, START_GAINED, SYNONYMOUS_CODING or SYNONYMOUS_STOP
* moderate impact: number of variants which an effect annotated in one of these categories: NON_SYNONYMOUS_CODING, CODON_CHANGE, CODON_INSERTION, CODON_CHANGE_PLUS_CODON_INSERTION, CODON_DELETION, CODON_CHANGE_PLUS_CODON_DELETION, UTR_5_DELETED or UTR_3_DELETED
* modifier impact: number of variants which an effect annotated in one of the other categories

#### Cohort SNV Metrics

This section gives an overview of the different metrics estimated using the overall set of samples. These metrics describe either type, effects, localisation or quality of SNV changes.

![Variant Count Distribution as Function of Mapping Quality ([download data table](SNV.SNVQuality.tsv))](SNV.SNVQuality.jpeg)

---

![Variant Count Distribution as Function of Read Coverage ([download data table](SNV.SNVCoverage.tsv))](SNV.SNVCoverage.jpeg)

---

![INDEL Count Distribution as Function of their Length ([download data table](SNV.IndelLength.tsv))](SNV.IndelLength.jpeg)

---

![Variant Spatial Distribution over Chromosomes using Gene Structure as Reference ([download data table](SNV.CountRegions.tsv))](SNV.CountRegions.jpeg)

[Additionally, the genomic variant number per Mb representation can found here](SNV.chromosomeChange.zip)

---

![Variant Effect Distribution ([download data table](SNV.CountEffects.tsv))](SNV.CountEffects.jpeg)

---

![Base Changes induced by SNV ([download data table](SNV.BaseChange.tsv))](SNV.BaseChange.jpeg)

---

![Codon Changes induced by SNV ([download data table](SNV.codonChange.tsv))](SNV.codonChange.jpeg)

---

![Amino Acid Changes induced by SNV ([download data table](SNV.AminoAcidChange.tsv))](SNV.AminoAcidChange.jpeg)


#### Sample SNV Metrics

This section gives an overview of the different metrics estimated using each sample individually. These metrics describe the density of changes in the chromosiome and the general type of change found in each sample.

![Variant Change Rate (i.e. Mean Distance in bp between two Variants) for each Chromosome and each Sample ([download data table](SNV.changeRate.tsv))](SNV.changeRate.jpeg)

---

![Transitions - Transvertions rate for each Sample ([download data table](SNV.TsTv.tsv))](SNV.TsTv.jpeg)
