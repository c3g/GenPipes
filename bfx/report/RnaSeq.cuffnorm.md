### Transcriptome Assembly with Cufflinks

In order to provide the most complete analysis results the entire set of Cufflinks files are provided in the following archive:

[Complete Cuff suite analysis files are available here](cuffAnalysis.zip)

A reference-based transcript assembly was performed, which allows the detection of known and novel transcripts isoforms. Transcript assembly is accomplished using [Cufflinks]\ [@cufflinks].

Cufflinks constructs a parsimonious set of transcripts that "explains" the reads observed in an RNA-Seq experiment. The assembly algorithm explicitly handles paired-end reads by treating the alignment for a given pair as a single object in the covering relation. Cufflinks tries to find the most parsimonious set of transcripts by performing a minimum-cost maximum matching. Reference annotation based assembly seeks to build upon available information about the transcriptome of an organism to find novel genes and isoforms. When a reference GTF is provided, the reference transcripts are tiled with faux-reads that will aid in the assembly of novel isoforms. These faux-reads are combined with the sequencing reads and are input into the regular Cufflinks assembler. The assembled transfrags are then compared to the reference transcripts to determine if they are different enough to be considered novel. Individual results of these assemblies are available in the `cufflinks/` folder in the archive.

[Cuffmerge] was then used to merge all assemblies gtf to a single one (`cufflinks/AllSamples/merged.gtf`).

The resulting merged assembly gtf file was used as a reference to estimate the abundance of each transcript and to subsequently perform a differential analysis for each design as provided in the design file using [Cuffdiff]. The entire set of result files for the Cuffdiff analysis is provided in the `cuffdiff/` folder of the archive.

As Cuffdiff generates normalized data using only the sample implicated in the design comparison, we additionally run the [Cuffnorm] tool to generate a normalized data set that includes all the samples. The resulting tables of normalized expression values are provided in the `cuffnorm/` folder of the archive.
