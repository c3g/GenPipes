#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import logging
import math
import os
import re
import sys

# GenPipes Modules
from ...core.config import global_conf
from ...core.job import Job, concat_jobs, pipe_jobs
from ...bfx.sample_tumor_pairs import parse_tumor_pair_file
from ...bfx.sequence_dictionary import split_by_size, parse_sequence_dictionary_file
import ...utils.utils

import gzip
from sys import stderr
from ..dnaseq import dnaseq

# utilizes
from ...bfx import sambamba
from ...bfx import bcftools
from ...bfx import tools
from ...bfx import metric_tools
from ...bfx import bvatools
from ...bfx import vt
from ...bfx import snpeff
from ...bfx import vawk
from ...bfx import deliverables
from ...bfx import bash_cmd as bash

# metrics
from ...bfx import conpair
from ...bfx import qualimap
from ...bfx import adapters
from ...bfx import fastqc
from ...bfx import multiqc

# variants
from ...bfx import htslib
from ...bfx import samtools
from ...bfx import varscan
from ...bfx import gatk
from ...bfx import gatk4
from ...bfx import vardict
from ...bfx import strelka2
from ...bfx import bcbio_variation_recall
from ...bfx import gemini

# sv
from ...bfx import delly
from ...bfx import manta
from ...bfx import lumpy
from ...bfx import svtyper
from ...bfx import wham
from ...bfx import metasv
from ...bfx import cnvkit
from ...bfx import scones
from ...bfx import sequenza
from ...bfx import amber
from ...bfx import cobalt
from ...bfx import purple
from ...bfx import svaba
from ...bfx import annotations

log = logging.getLogger(__name__)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = TumorPair.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = TumorPair.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

if __name__ == '__main__':
    main()
