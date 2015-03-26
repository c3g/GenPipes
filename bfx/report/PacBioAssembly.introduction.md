Introduction
------------

This document contains the description of the current MUGQIC PacBio Assembly analysis. The information presented here reflects the current state of the analysis as of $date$.

Assembly Method
---------------

Contigs assembly was done using what is referred as the [HGAP] workflow [@hgap]. Briefly, raw subreads where generated from raw .bax.h5 PacBio data files. A subread length cutoff value was extracted from subreads and used into the preassembly ([BLASR]\ [@blasr]) step which consists of aligning short subreads on long subreads. Since errors in PacBio reads is random, the alignment of multiple short reads on longer reads allows to correct sequencing error on long reads. These long corrected reads are then used as seeds into assembly ([Celera Assembler]\ [@celera]) which gives contigs. These contigs are then 'polished' by aligning raw reads on contigs (BLASR) that are then processed through a variant calling algorithm (Quiver) that generates high quality consensus sequences using local realignments and PacBio quality scores. For more information, visit the [PacBio] publications web page.

Analysis and Results
--------------------
