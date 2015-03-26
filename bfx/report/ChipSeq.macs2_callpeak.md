### Peak Calls

Peaks are called using [MACS] software [@macs]. The table below shows the general peak calling strategies for narrow and wide peaks. The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run. The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100. The default mfold parameter of MACS is [10,30].

Table: Peak Calling Strategies using MACS

Design Type | Narrow Peaks | Broad Peaks
------------|--------------|------------
Treatment without Control | --mfold=mfold --fix-bimodal --nolambda | --mfold=mfold --nomodel --nolambda
Treatment with Control | --mfold=mfold --nomodel | --mfold=mfold --nomodel

The following links will direct you to the MACS peak calls for each design:

