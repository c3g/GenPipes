#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from bio.readset import *

from bio import smrtanalysis
from pipelines import common

log = logging.getLogger(__name__)

class PacBioAssembly(common.MUGQICPipeline):

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_pacbio_readset_file(self.args.readsets.name)
        return self._readsets

    """
    Filtering. This step uses smrtpipe.py (From the SmrtAnalysis package) and will filter reads and subreads  based on their length and QVs.
    1- fofnToSmrtpipeInput.py. 
    2- modify RS_Filtering.xml files according to reads filtering values entered in .ini file.
    3- smrtpipe.py with filtering protocol
    4- prinseq-lite.pl to write fasta file based on fastq file.
    Informative run metrics such as loading efficiency, readlengths and base quality are generated in this step as well.
    """
    def smrtanalysis_filtering(self):
        jobs = []

        for sample in self.samples:
            fofn = os.path.join("fofns", sample.name + ".fofn")
            bax_files = [bax_file for readset in sample.readsets for bax_file in readset.bax_files]
            filtering_directory = os.path.join(sample.name, "filtering")

            jobs.append(concat_jobs([
                Job(command="mkdir -p fofns"),
                Job(bax_files, [fofn], command="""\
\`cat > {fofn} << END
{bax_files}
END
\`""".format(bax_files="\n".join(bax_files), fofn=fofn)),
                Job(command="mkdir -p " + filtering_directory),
                smrtanalysis.filtering(
                    fofn,
                    os.path.join(filtering_directory, "input.xml"),
                    os.path.join(sample.name, "filtering.xml"),
                    filtering_directory,
                    os.path.join(filtering_directory, "smrtpipe.log")
                )
            ], name="filtering." + sample.name))

        return jobs

    @property
    def steps(self):
        return [
            self.smrtanalysis_filtering
        ]

if __name__ == '__main__': 
    PacBioAssembly().submit_jobs()
