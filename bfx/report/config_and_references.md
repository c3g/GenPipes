## Analysis Configuration Parameters

All analysis parameters are described in this [configuration file](config.ini).

## References
<span />

[BLASR]: https://github.com/PacificBiosciences/blasr
[BLAST]: http://www.ncbi.nlm.nih.gov/books/NBK52640/
[BWA]: http://bio-bwa.sourceforge.net/
[Celera Assembler]: http://wgs-assembler.sourceforge.net/
[Cuffdiff]: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/
[Cufflinks]: http://cole-trapnell-lab.github.io/cufflinks/
[Cuffmerge]: http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/
[Cuffnorm]: http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/
[DESeq]: http://bioconductor.org/packages/release/bioc/html/DESeq.html
[edgeR]: http://bioconductor.org/packages/release/bioc/html/edgeR.html
[GATK]: https://www.broadinstitute.org/gatk/
[goseq]: http://www.bioconductor.org/packages/release/bioc/html/goseq.html
[HGAP]: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP
[HOMER]: http://homer.salk.edu/homer/
[HTSeq]: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
[IGV]: https://www.broadinstitute.org/igv/
[MACS]: https://github.com/taoliu/MACS
[PacBio]: http://www.pacificbiosciences.com/news_and_events/publications/
[Picard]: http://broadinstitute.github.io/picard/
[RSEM]: http://deweylab.biostat.wisc.edu/rsem/
[Sambamba]: https://lomereiter.github.io/sambamba/
[SAMtools]: http://www.htslib.org/
[SnpEff]: http://snpeff.sourceforge.net/
[SnpSift]: http://snpeff.sourceforge.net/SnpSift.html
[STAR]: https://github.com/alexdobin/STAR
[Trimmomatic]: http://www.usadellab.org/cms/index.php?page=trimmomatic
[Trinity]: https://trinityrnaseq.github.io/
[UCSC]: http://genome.ucsc.edu/
[FLASh]: http://ccb.jhu.edu/software/FLASH/
[Uchime]: http://drive5.com/usearch/usearch_docs.html
[Qiime]: http://qiime.org/
[Usearch]: http://drive5.com/usearch/usearch_docs.html
[Uclust]: http://drive5.com/usearch/manual/uclust_algo.html
[PyNAST]: http://biocore.github.io/pynast/
[FastTree]: http://www.microbesonline.org/fasttree/
[Krona]: http://sourceforge.net/p/krona/home/krona/
[VSEARCH]: https://github.com/torognes/vsearch
[CSS]: http://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html
[GOLD]: http://microbiomeutil.sourceforge.net/
[UNITE]: https://unite.ut.ee/repository.php
[BISMARK]: https://www.bioinformatics.babraham.ac.uk/projects/bismark/
[spp]: https://github.com/kundajelab/phantompeakqualtools

---
references:
- id: blasr
  title: 'Mapping single molecule sequencing reads using basic local alignment with successive refinement (BLASR): application and theory'
  author:
    - family: Chaisson
      given: M.J.
    - family: Tesler
      given: G.
  container-title: BMC Bioinformatics
  volume: 13
  URL: 'http://www.biomedcentral.com/1471-2105/13/238'
  DOI: 10.1186/1471-2105-13-238
  issue: 238
  type: article-journal
  issued:
    year: 2012

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

- id: celera
  title: 'A Whole-Genome Assembly of Drosophila'
  author:
    - family: Myers
      given: E.W.
    - family: . Sutton
      given: G.G.
    - family: Delcher
      given: A.L.
    - family: Dew
      given: I.M.
    - family: Fasulo
      given: D.P.
    - family: Flanigan
      given: M.J.
    - family: Kravitz
      given: S.A.
    - family: Mobarry
      given: C.M.
    - family: Reinert
      given: K.H.J.
    - family: Remington
      given: K.A.
    - family: Anson
      given: E.L.
    - family: Bolanos
      given: R.A.
    - family: Chou
      given: H.-H.
    - family: Jordan
      given: C.M.
    - family: Halpern
      given: A.L.
    - family: Lonardi
      given: S.
    - family: Beasley
      given: E.M.
    - family: Brandon
      given: R.C.
    - family: Chen
      given: L.
    - family: Dunn
      given: P.J.
    - family: Lai
      given: Z.
    - family: Liang
      given: Y.
    - family: Nusskern
      given: D.R.
    - family: Zhan
      given: M.
    - family: Zhang
      given: Q.
    - family: Zheng
      given: X.
    - family: Rubin
      given: G.M.
    - family: Adams
      given: M.D.
    - family: Venter
      given: J.C.
  container-title: Science
  volume: 287
  URL: 'http://www.sciencemag.org/content/287/5461/2196'
  DOI:  10.1126/science.287.5461.2196
  issue: 5461
  page: 2196-2204
  type: article-journal
  issued:
    year: 2000

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

- id: diginorm
  title: 'A Reference-Free Algorithm for Computational Normalization of Shotgun Sequencing Data'
  author:
    - family: Titus Brown
      given: C.
    - family: Howe
      given: A.
    - family: Zhang
      given: Q.
    - family: Pyrkosz
      given: A.B.
    - family: Brom
      given: T.H.
  container-title: arXiv.org
  URL: 'http://arxiv.org/abs/1203.4802'
  issued:
    year: 2012

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

- id: hgap
  title: 'Nonhybrid, finished microbial genome assemblies from long-read SMRT sequencing data'
  author:
    - family: Chin
      given: C.-S.
    - family: Alexander
      given: D.H.
    - family: Marks
      given: P.
    - family: Klammer
      given: A.A.
    - family: Drake
      given: J.
    - family: Heiner
      given: C.
    - family: Clum
      given: A.
    - family: Copeland
      given: A.
    - family: Huddleston
      given: J.
    - family: Eichler
      given: E.E.
    - family: Turner
      given: S.W.
    - family: Korlach
      given: J.
  container-title: Nature Methods
  volume: 10
  URL: 'http://www.nature.com/nmeth/journal/v10/n6/full/nmeth.2474.html'
  DOI: 10.1038/nmeth.2474
  issue: 6
  page: 563-569
  type: article-journal
  issued:
    year: 2013

- id: homer
  title: 'Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities'
  author:
    - family: Heinz
      given: S.
    - family: Benner
      given: C.
    - family: Spann
      given: N.
    - family: Bertolino
      given: E.
    - family: Lin
      given: Y.C.
    - family: Laslo
      given: P.
    - family: Cheng
      given: J.X.
    - family: Murre
      given: C.
    - family: Singh
      given: H.
    - family: Glass
      given: C.K.
  container-title: Molecular Cell
  volume: 38
  URL: 'http://www.sciencedirect.com/science/article/pii/S1097276510003667'
  DOI: 10.1016/j.molcel.2010.05.004
  issue: 4
  page: 576-589
  type: article-journal
  issued:
    year: 2010

- id: macs
  title: 'Model-based Analysis of ChIP-Seq (MACS)'
  author:
    - family: Zhang
      given: Y.
    - family: Liu
      given: T.
    - family: Meyer
      given: C.A.
    - family: Eeckhoute
      given: J.
    - family: Johnson
      given: D.S.
    - family: Bernstein
      given: B.E.
    - family: Nusbaum
      given: C.
    - family: Myers
      given: R.M.
    - family: Brown
      given: M.
    - family: Li
      given: W.
    - family: Liu
      given: X.S.
  container-title: Genome Biology
  volume: 9
  URL: 'http://genomebiology.com/2008/9/9/R137'
  DOI: 10.1186/gb-2008-9-9-r137
  issue: 9
  page: R137
  type: article-journal
  issued:
    year: 2008

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

- id: rsem
  title: 'RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome'
  author:
    - family: Li
      given: B.
    - family: Dewey
      given: C.N.
  container-title: BMC Bioinformatics
  volume: 12
  URL: 'http://www.biomedcentral.com/1471-2105/12/323'
  DOI: 10.1186/1471-2105-12-323
  issue: 323
  type: article-journal
  issued:
    year: 2011

- id: sambamba
  title: 'Sambamba: fast processing of NGS alignment formats'
  author:
    - family: Tarasov
      given: Artem
    - family: Vilella
      given: Albert J.
    - family: Cuppen
      given: Edwin
    - family: Nijman
      given: Isaac J.
    - family: Prins
      given: Pjotr
  container-title: Bioinformatics
  volume: 31
  URL: 'http://dx.doi.org/10.1093/bioinformatics/btv098'
  DOI: 10.1093/bioinformatics/btv098
  issue: 12
  page: 2032-2034
  type: article-journal
  issued:
    year: 2015

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

- id: trinity_protocol
  title: 'De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis'
  author:
    - family: Haas,
      given: B.J.
    - family: Papanicolaou,
      given: A.
    - family: Yassour,
      given: M.
    - family: Grabherr,
      given: M.
    - family: Blood,
      given: P.D.
    - family: Bowden,
      given: J.
    - family: Couger,
      given: M.B.
    - family: Eccles,
      given: D.
    - family: Li,
      given: B.
    - family: Lieber,
      given: M.
    - family: MacManes,
      given: M.D.
    - family: Ott,
      given: M.
    - family: Orvis,
      given: J.
    - family: Pochet,
      given: N.
    - family: Strozzi,
      given: F.
    - family: Weeks,
      given: N.
    - family: Westerman,
      given: R.
    - family: William,
      given: T.
    - family: Dewey,
      given: C.N.
    - family: Henschel,
      given: R.
    - family: LeDuc
      given: R.D.
    - family: Friedman
      given: N.
    - family: Regev
      given: A.
  container-title: Nature Protocols
  volume: 8
  URL: 'http://www.nature.com/nprot/journal/v8/n8/full/nprot.2013.084.html'
  DOI: 10.1038/nprot.2013.084
  issue: 8
  page: 1494-1512
  type: article-journal
  issued:
    year: 2013

- id: trinity
  title: 'Full-length transcriptome assembly from RNA-Seq data without a reference genome'
  author:
    - family: Grabherr,
      given: M.G.
    - family: Haas,
      given: B.J.
    - family: Yassour,
      given: M.
    - family: Levin,
      given: J.Z.
    - family: Thompson,
      given: D.A.
    - family: Amit,
      given: I.
    - family: Adiconis,
      given: X.
    - family: Fan,
      given: L.
    - family: Raychowdhury,
      given: R.
    - family: Zeng,
      given: Q.
    - family: Chen,
      given: Z.
    - family: Mauceli,
      given: E.
    - family: Hacohen,
      given: N.
    - family: Gnirke,
      given: A.
    - family: Rhind,
      given: N.
    - family: di Palma,
      given: F.
    - family: Birren,
      given: B.W.
    - family: Nusbaum,
      given: C.
    - family: Lindblad-Toh,
      given: K.
    - family: Friedman
      given: N.
    - family: Regev
      given: A.
  container-title: Nature Biotechnology
  volume: 29
  URL: 'http://www.nature.com/nbt/journal/v29/n7/abs/nbt.1883.html'
  DOI: 10.1038/nbt.1883
  issue: 7
  page: 644-652
  type: article-journal
  issued:
    year: 2011

- id: flash
  title: 'FLASH: Fast length adjustment of short reads to improve genome assemblies'
  author:
    - family: Magoc
      given: T.
    - family: Salzberg
      given: S.       
  container-title: Bioinformatics
  volume: 27
  URL: 'http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507.short'
  DOI: 10.1093/bioinformatics/btr507
  issue: 21
  page: 2957-2963
  type: article-journal
  issued:
    year: 2011      
        
- id: uchime
  title: 'UCHIME improves sensitivity and speed of chimera detection'
  author:
    - family: Edgar
      given: R.C.
    - family: Haas
      given: B.J.
    - family: Clemente
      given: JC.
    - family: Quince
      given: C.
    - family: Knight
      given: R.          
  container-title: Bioinformatics
  volume: 27
  URL: 'http://bioinformatics.oxfordjournals.org/content/early/2011/06/23/bioinformatics.btr381.short?rss=1'
  DOI: 10.1093/bioinformatics/btr381
  issue: 16
  page: 2194-2200
  type: article-journal
  issued:
    year: 2011     
    
- id: qiime
  title: 'QIIME allows analysis of high-throughput community sequencing data'
  author:
    - family: Caporaso
      given: J.G.
    - family: Kuczynski
      given: J.
    - family: Stombaugh
      given: J.
    - family: Bittinger
      given: K.
    - family: Bushman
      given: F.D.
    - family: Costello
      given: E.K.
    - family: Fierer
      given: N.
    - family: Gonzalez Pena
      given: A.
    - family: Goodrich
      given: J.K.
    - family: Gordon
      given: J.I.
    - family: Huttley
      given: G.A.
    - family: Kelley
      given: S.T.
    - family: Knights
      given: D.
    - family: Koenig
      given: J.E.
    - family: Ley
      given: R.E.
    - family: Lozupone
      given: C.A.
    - family: McDonald
      given: D.
    - family: Muegge
      given: B.D.
    - family: Pirrung
      given: M.
    - family: Reeder
      given: J.
    - family: Sevinsky
      given: J.R.
    - family: Turnbaugh
      given: P.J.
    - family: Walters
      given: W.A.
    - family: Widmann
      given: J.
    - family: Yatsunenko
      given: T.
    - family: Zaneveld
      given: J.
    - family: Knight
      given: R. 
  container-title: Nature Methods
  volume: 7
  URL: 'http://www.nature.com/nmeth/journal/v7/n5/full/nmeth.f.303.html'
  DOI: 10.1038/nmeth.f.303
  issue: 5
  page: 335-336
  type: article-journal
  issued:
    year: 2010 

- id: usearch
  title: 'Search and clustering orders of magnitude faster than BLAST'
  author:
    - family: Edgar
      given: R.C.         
  container-title: Bioinformatics
  volume: 26
  URL: 'http://bioinformatics.oxfordjournals.org/content/26/19/2460'
  DOI: 10.1093/bioinformatics/btq461
  issue: 19
  page: 2460-2461
  type: article-journal
  issued:
    year: 2010  
    
- id: uclust
  title: 'Search and clustering orders of magnitude faster than BLAST'
  author:
    - family: Edgar
      given: R.C.         
  container-title: Bioinformatics
  volume: 26
  URL: 'http://bioinformatics.oxfordjournals.org/content/26/19/2460'
  DOI: 10.1093/bioinformatics/btq461
  issue: 19
  page: 2460-2461
  type: article-journal
  issued:
    year: 2010  
    
- id: pynast
  title: 'PyNAST: a flexible tool for aligning sequences to a template alignment'
  author:
    - family: Caporaso
      given: J.G.
    - family: Bittinger
      given: K.
    - family: Bushman
      given: F.D.
    - family: DeSantis
      given: T.Z.
    - family: Andersen
      given: G.L.      
    - family: Knight
      given: R.           
  container-title: Bioinformatics
  volume: 26
  URL: 'http://bioinformatics.oxfordjournals.org/content/26/2/266'
  DOI: 10.1093/bioinformatics/btp636
  issue: 2
  page: 266-267
  type: article-journal
  issued:
    year: 2011  

- id: fasttree
  title: 'FastTree: computing large minimum evolution trees with profiles instead of a distance matrix'
  author:
    - family: Price
      given: M.N.
    - family: Dehal
      given: P.S.
    - family: Arkin
      given: A.P.       
  container-title: Mol Biol Evol
  volume: 26
  URL: 'http://mbe.oxfordjournals.org/content/26/7/1641'
  DOI: 10.1093/molbev/msp077
  issue: 7
  page: 1641-1650
  type: article-journal
  issued:
    year: 2009

- id: krona
  title: 'Interactive metagenomic visualization in a Web browser'
  author:
    - family: Ondov
      given: B.D.
    - family: Bergman
      given: N.H.
    - family: Phillippy
      given: A.M.       
  container-title: BMC Bioinformatics
  volume: 12
  URL: 'http://www.biomedcentral.com/1471-2105/12/385'
  DOI: 10.1186/1471-2105-12-385
  issue: 1
  page: 385
  type: article-journal
  issued:
    year: 2011

- id: vsearch
  title: 'VSEARCH'
  author:
    - family: Flouri
      given: T.
    - family: Zeeshan
      given: U.
    - family: Mahé
      given: F.       
    - family: Nichols
      given: B.   
    - family: Quince
      given: C. 
    - family: Rognes
      given: T. 
  URL: 'https://github.com/torognes/vsearch'
  DOI: 10.5281/zenodo.15524           
  issued:
    year: 2015

- id: css
  title: 'Differential abundance analysis for microbial marker-gene surveys'
  author:
    - family: Paulson
      given: J.N.
    - family: Stine
      given: O.C.
    - family: Bravo
      given: H.C.   
    - family: Pop
      given: M.     
  container-title: Nat Methods
  volume: 12
  URL: 'http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2658.html'
  DOI: 10.1038/nmeth.2658
  page: 1200
  type: article-journal
  issued:
    year: 2013   

- id: gold
  title: 'Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons'
  author:
    - family: Haas
      given: B.J.
    - family: Gevers
      given: D.
    - family: Earl
      given: A.   
    - family: Feldgarden
      given: M.   
    - family: Ward
      given: D.V.   
    - family: Giannokous
      given: G.
    - family: Ciulla
      given: D.
    - family: Tabbaa
      given: D.   
    - family: Highlander
      given: S.K.   
    - family: Sodergren
      given: E.  
    - family: Methe
      given: B.   
    - family: Desantis
      given: T.Z.   
    - family: Petrosino
      given: J.F.
    - family: Knight
      given: R.
    - family: Birren
      given: B.W.   
  container-title: Genome Res
  volume: 21
  URL: 'http://genome.cshlp.org/content/21/3/494.long'
  DOI: 10.1101/gr.112730.110
  page: 494
  type: article-journal
  issued:
    year: 2011     
    
- id: unite
  title: 'A Comprehensive, Automatically Updated Fungal ITS Sequence Dataset for Reference-Based Chimera Control in Environmental Sequencing Efforts'
  author:
    - family: Nilsson
      given: R.H.
    - family: Tedersoo
      given: L.
    - family: Ryberg
      given: M.   
    - family: Kristiansson
      given: E.   
    - family: Hartmann
      given: M.   
    - family: Unterseher
      given: M.
    - family: Porter
      given: T.M.
    - family: Bengtsson-Palme
      given: J.   
    - family: Walker
      given: D.M.   
    - family: de Sousa
      given: F.  
    - family: Gamper
      given: H.A.   
    - family: Larsson
      given: E.   
    - family: Larsson
      given: K.H.
    - family: Kõljalg
      given: U.
    - family: Edgar
      given: R.C.   
    - family: Abarenkov
      given: K.  
  container-title: Microbes Environ
  volume: 30
  URL: 'https://www.jstage.jst.go.jp/article/jsme2/30/2/30_ME14121/_article'
  DOI: 10.1264/jsme2.ME14121
  page: 145
  type: article-journal
  issued:
    year: 2015     

- id: bismark
  title: 'Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.'
  author:
    - family: Krueger
      given: Felix
    - family: Andrews
      given: Simon R.
  container-title: Bioinformatics
  volume: 27
  URL: 'https://www.ncbi.nlm.nih.gov/pubmed/21493656'
  DOI: 10.1093/bioinformatics/btr167
  page: 1571
  type: article-journal
  issued:
    year: 2011

- id: spp
  title: 'ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia.'
  author:
    - family: Kheradpour
      given: Pouya
    - family: Pauli
      given: Florencia
    - family: Batzoglou
      given: Serafim
    - family: Bernstein
      given: Bradley E
    - family: Bickel
      given: Peter
    - family: Brown
      given: James B
    - family: Cayting
      given: Philip
    - family: Chen
      given: Yiwen
    - family: DeSalvo
      given: Gilberto
    - family: Epstein
      given: Charles
    - family: Fisher-Aylor
      given: Katherine I
    - family: Euskirchen
      given: Ghia
    - family: Gerstein
      given: Mark
    - family: Gertz
      given: Jason
    - family: Hartemink
      given: Alexander J
    - family: Hoffman
      given: Michael M
    - family: Iyer
      given: Vishwanath R
    - family: Jung
      given: Youngsook L
    - family: Karmakar
      given: Subhradip
    - family: Kellis
      given: Manolis
    - family: Kharchenko
      given: Peter V
    - family: Li
      given: Qunhua
    - family: Liu
      given: Tao
    - family: Liu
      given: X Shirley
    - family: Ma
      given: Lijia
    - family: Milosavljevic
      given: Aleksandar
    - family: Myers
      given: Richard M
    - family: Park
      given: Peter J
    - family: Pazin
      given: Michael J
    - family: Perry
      given: Marc D
    - family: Raha
      given: Debasish
    - family: Reddy
      given: Timothy E
    - family: Rozowsky
      given: Joel
    - family: Shoresh
      given: Noam
    - family: Sidow
      given: Arend
    - family: Slattery
      given: Matthew
    - family: Stamatoyannopoulos
      given: John A
    - family: Tolstorukov
      given: Michael Y
    - family: White
      given: Kevin P
    - family: Xi
      given: Simon
    - family: Farnham
      given: Peggy J
    - family: Lieb
      given: Jason D
    - family: Wold
      given: Barbara J
    - family: Snyder
      given: Michael
  container-title: Genome research
  volume: 22
  URL: 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/'
  DOI: 10.1101/gr.136184.111
  page: 1813-31
  type: article-journal
  issued:
    year: 2012
                                   
...

