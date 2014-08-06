#!/usr/bin/env python

# Python Standard Modules
import logging
import math
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *
from bio.sequence_dictionary import *

from bio import bvatools
from bio import bwa
from bio import gatk
from bio import gq_seq_utils
from bio import igvtools
from bio import metrics
from bio import picard
from bio import samtools
from bio import snpeff
from bio import tools
from bio import vcftools
from pipelines.illumina import illumina

log = logging.getLogger(__name__)

class Puure(dnaseq.DnaSeq):


    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.bwa_mem_picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.picard_mark_duplicates,
            self.metrics,
        ]

Puure().submit_jobs()
