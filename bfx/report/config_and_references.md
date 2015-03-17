## Analysis Configuration Parameters

All analysis parameters are described in this [configuration file](config.ini).

## References
<span />

[BWA]: http://bio-bwa.sourceforge.net/
[Cuffdiff]: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/
[Cufflinks]: http://cole-trapnell-lab.github.io/cufflinks/
[Cuffmerge]: http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/
[Cuffnorm]: http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/
[DESeq]: http://bioconductor.org/packages/release/bioc/html/DESeq.html
[edgeR]: http://bioconductor.org/packages/release/bioc/html/edgeR.html
[GATK]: https://www.broadinstitute.org/gatk/
[goseq]: http://www.bioconductor.org/packages/release/bioc/html/goseq.html
[HTSeq]: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
[IGV]: https://www.broadinstitute.org/igv/
[Picard]: http://broadinstitute.github.io/picard/
[SAMtools]: http://www.htslib.org/
[SnpEff]: http://snpeff.sourceforge.net/
[SnpSift]: http://snpeff.sourceforge.net/SnpSift.html
[STAR]: https://github.com/alexdobin/STAR
[Trimmomatic]: http://www.usadellab.org/cms/index.php?page=trimmomatic
[UCSC]: http://genome.ucsc.edu/

---
references:
- id: bwa
  title: 'Fast and accurate long-read alignment with Burrows–Wheeler transform'
  author:
    - family: Li
      given: H.
    - family: Durbin
      given: R.
  container-title: Bioinformatics
  volume: 26
  URL: 'http://bio-bwa.sourceforge.net/'
  DOI: 10.1093/bioinformatics/btp698
  issue: 5
  page: 589-595
  type: article-journal
  issued:
    year: 2010

- id: cuffdiff
  title: 'Differential analysis of gene regulation at transcript resolution with RNA-seq'
  author:
    - family: Trapnell
      given: C.
    - family: Hendrickson
      given: D.G.
    - family: Sauvageau
      given: M.
    - family: Goff
      given: L.
    - family: Rinn
      given: J.L.
    - family: Pachter
      given: L.
  container-title: Nature Biotechnology
  volume: 31
  URL: 'http://www.nature.com/nbt/journal/v31/n1/full/nbt.2450.html'
  DOI: 10.1038/nbt.2450
  issue: 1
  page: 46-53
  type: article-journal
  issued:
    year: 2012

- id: cufflinks
  title: 'Identification of novel transcripts in annotated genomes using RNA-Seq'
  author:
    - family: Roberts
      given: A.
    - family: Pimentel
      given: H.
    - family: Trapnell
      given: C.
    - family: Pachter
      given: L.
  container-title: Bioinformatics
  volume: 27
  URL: 'http://bioinformatics.oxfordjournals.org/content/27/17/2325'
  DOI: 10.1093/bioinformatics/btr355
  issue: 17
  page: 2325-2329
  type: article-journal
  issued:
    year: 2011

- id: deseq
  title: 'Differential expression analysis for sequence count data'
  author:
    - family: Anders
      given: S.
    - family: Huber
      given: W.
  container-title: Genome Biology
  volume: 11
  URL: 'http://genomebiology.com/2010/11/10/R106'
  DOI: 10.1186/gb-2010-11-10-r106
  issue: 10
  page: R106
  type: article-journal
  issued:
    year: 2010

- id: edger
  title: 'edgeR: a Bioconductor package for differential expression analysis of digital gene expression data'
  author:
    - family: Robinson
      given: M.D.
    - family: McCarthy
      given: D.J.
    - family: Smyth
      given: G.K.
  container-title: Bioinformatics
  volume: 26
  URL: 'http://bioinformatics.oxfordjournals.org/content/26/1/139'
  DOI: 10.1093/bioinformatics/btp616
  issue: 1
  page: 139-140
  type: article-journal
  issued:
    year: 2010

- id: gatk
  title: 'A framework for variation discovery and genotyping using next-generation DNA sequencing data'
  author:
    - family: DePristo et al.
  container-title: Nature Genetics
  volume: 43
  URL: 'http://www.nature.com/ng/journal/v43/n5/full/ng.806.html'
  DOI: 10.1038/ng.806
  issue: 5
  page: 491-498
  type: article-journal
  issued:
    year: 2011

- id: goseq
  title: 'Gene ontology analysis for RNA-seq: accounting for selection bias'
  author:
    - family: Young
      given: M.D.
    - family: Wakefield
      given: M.J.
    - family: Smyth
      given: G.K.
    - family: Oshlack
      given: A.
  container-title: Genome Biology
  volume: 11
  URL: 'http://genomebiology.com/2010/11/2/R14'
  DOI: 10.1186/gb-2010-11-2-r14
  issue: 2
  page: R14
  type: article-journal
  issued:
    year: 2010

- id: quantifying_rnaseq
  title: 'Mapping and quantifying mammalian transcriptomes by RNA-Seq'
  author:
    - family: Mortazavi
      given: A.
    - family: Williams
      given: B.A.
    - family: McCue
      given: K.
    - family: Schaeffer
      given: L.
    - family: Wold
      given: B.
  container-title: Nature Methods
  volume: 5
  URL: 'http://www.nature.com/nmeth/journal/v5/n7/full/nmeth.1226.html'
  DOI: 10.1038/nmeth.1226
  issue: 7
  page: 621-628
  type: article-journal
  issued:
    year: 2008

- id: samtools
  title: 'The Sequence Alignment/Map format and SAMtools'
  author:
    - family: Li
      given: H.
    - family: Handsaker
      given: B.
    - family: Wysoker
      given: A.
    - family: Fennell
      given: T.
    - family: Ruan
      given: J.
    - family: Homer
      given: N.
    - family: Marth
      given: G.
    - family: Abecasis
      given: G.
    - family: Durbin
      given: R.
    - family: 1000 Genome Project Data Processing Subgroup
  container-title: Bioinformatics
  volume: 25
  URL: 'http://bioinformatics.oxfordjournals.org/content/25/16/2078'
  DOI: 10.1093/bioinformatics/btp352
  issue: 16
  page: 2078-2079
  type: article-journal
  issued:
    year: 2009

- id: snpeff
  title: 'A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3'
  author:
    - family: Cingolani
      given: P.
    - family: Platts
      given: A.
    - family: Wang
      given: Le L.
    - family: Coon
      given: M.
    - family: Nguyen
      given: T.
    - family: Wang
      given: L.
    - family: Land
      given: S.J.
    - family: Lu
      given: X.
    - family: Ruden
      given: D.M.
  container-title: Fly
  volume: 6
  URL: 'http://www.tandfonline.com/doi/abs/10.4161/fly.19695'
  DOI: 10.4161/fly.19695
  issue: 2
  page: 80-92
  type: article-journal
  issued:
    year: 2012

- id: snpsift
  title: 'Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift'
  author:
    - family: Cingolani
      given: P.
    - family: Patel
      given: V.M.
    - family: Coon
      given: M.
    - family: Nguyen
      given: T.
    - family: Land
      given: S.J.
    - family: Ruden
      given: D.M.
    - family: Lu
      given: X.
  container-title: Frontiers in Genetics
  volume: 3
  URL: 'http://journal.frontiersin.org/article/10.3389/fgene.2012.00035/abstract'
  DOI: 10.3389/fgene.2012.00035
  type: article-journal
  issued:
    year: 2012

- id: star
  title: 'STAR: ultrafast universal RNA-seq aligner'
  author:
    - family: Dobin
      given: A.
    - family: Davis
      given: C.A.
    - family: Schlesinger
      given: F.
    - family: Drenkow
      given: J.
    - family: Zaleski
      given: C.
    - family: Jha
      given: S.
    - family: Batut
      given: P.
    - family: Chaisson
      given: M.
    - family: Gingeras
      given: T.R.
  container-title: Bioinformatics
  volume: 29
  URL: 'http://bioinformatics.oxfordjournals.org/content/29/1/15'
  DOI: 10.1093/bioinformatics/bts635
  issue: 1
  page: 15-21
  type: article-journal
  issued:
    year: 2012

- id: trimmomatic
  title: 'Trimmomatic: a flexible trimmer for Illumina sequence data'
  author:
    - family: Bolger
      given: A.M.
    - family: Lohse
      given: M.
    - family: Usadel
      given: B.
  container-title: Bioinformatics
  volume: 30
  URL: 'http://www.usadellab.org/cms/index.php?page=trimmomatic'
  DOI: 10.1093/bioinformatics/btu170
  issue: 15
  page: 2114-2120
  type: article-journal
  issued:
    year: 2014
...
