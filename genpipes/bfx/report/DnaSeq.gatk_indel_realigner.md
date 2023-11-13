### Read Alignment Improvement

Unfortunately, there is no prefect aligner to date. The alignment step is really sensitive to the aligner parameters as well as technical and biological variations. To increase the quality and speed of subsequent variant calling, we performed a series of alignment improvement procedures. These procedures consist of realigning the surrounding short insertion or deletion, fixing possible read mate discrepancy due to realignment, marking duplicated reads and recalibrating read base quality.

#### Realigning Short Insertions and Deletions (INDELs)

INDELs in reads (especially near the ends) can trick the mappers into mis-aligning with mismatches. These artifacts resulting from mismatches can harm base quality recalibration and variant detection. Realignment around INDELs helps to improve the accuracy of several downstream steps.

Insertion and deletion realignment is performed on regions where multiple base mismatches are preferred over INDELs by the aligner since it appears to be less costly for the algorithm. Such regions will introduce false positive variant calls which may be filtered out by realigning those regions properly. Mainly realignment will occurs in 3 different region types:

1. Known sites (e.g. dbSNP, 1000 Genomes)
2. INDELs seen in original alignments (in CIGARs)
3. Sites where evidence suggests a hidden INDEL

Realignment is done using [GATK]\ [@gatk].
