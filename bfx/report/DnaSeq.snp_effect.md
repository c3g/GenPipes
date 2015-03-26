#### SNV Annotations

An in-house database identifies regions in which reads are confidently mapped to the reference genome. Generally, low mappability corresponds nicely to RepeatMasker regions and decreases substantially with read length increase. A region is identified as HC = coverage too high, LC = low coverage, MQ = too low mean mapQ (< 20) and ND = unmappable region (no data) at the position.

The VCF files are also annotated for:

1. dbSNP using the [SnpSift] software [@snpsift]
2. Variant effects (predicted variant effects on gene, such as amino acid changes) using the [SnpEff] software [@snpeff]
