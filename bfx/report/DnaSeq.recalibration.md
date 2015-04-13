#### Base Quality Recalibration

This step allows recalibrating base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome. Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so provides not only more accurate quality scores but also more widely dispersed ones. The recalibration is done not only for the well-known base quality scores but also for base insertion and base deletion quality scores. These are per-base quantities which estimate the probability that the next base in the read was mis-incorporated or mis-deleted (due to slippage, for example).

The recalibration process is accomplished by analyzing the covariation among several features of a base. For example:

* Reported quality score
* The position within the read
* The preceding and current nucleotide (sequencing chemistry effect) observed by the sequencing machine

These covariates are then subsequently applied through a piecewise tabular correction to recalibrate the quality scores of all reads in a BAM file. The recalibration process is performed using [GATK]\ [@gatk].
