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

from bio import tools
from pipelines import common

log = logging.getLogger(__name__)

class PacBioAssembly(common.MUGQICPipeline):

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_pacbio_readset_file(self.args.readsets.name)
        return self._readsets

    def filtering(self):
        return [Job(command="\n".join(["BAX: " + "\t".join(readset.bax_files) for readset in self.readsets]), name="filtering")]

    @property
    def steps(self):
        return [
            self.filtering
        ]

if __name__ == '__main__': 
    PacBioAssembly().submit_jobs()
