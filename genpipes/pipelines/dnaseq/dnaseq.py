#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs, pipe_jobs
from .. import common
from ...bfx.sequence_dictionary import parse_sequence_dictionary_file, split_by_size
from ... import utils

from ...bfx import adapters
from ...bfx import bvatools
from ...bfx import bwa
from ...bfx import bwakit
from ...bfx import gatk4
from ...bfx import igvtools
from ...bfx import metrics
from ...bfx import samtools
from ...bfx import snpeff
from ...bfx import tools
from ...bfx import vcftools
from ...bfx import skewer
from ...bfx import sambamba
from ...bfx import picard2
from ...bfx import vt
from ...bfx import htslib
from ...bfx import gemini
from ...bfx import qualimap
from ...bfx import fastqc
from ...bfx import ngscheckmate
from ...bfx import peddy
from ...bfx import verify_bam_id
from ...bfx import multiqc
from ...bfx import deliverables

from ...bfx import bash_cmd as bash

##Structural variants
from ...bfx import delly
from ...bfx import manta
from ...bfx import lumpy
from ...bfx import wham
from ...bfx import cnvkit
from ...bfx import breakseq2
from ...bfx import metasv
from ...bfx import svtyper
from ...bfx import bcftools
from ...bfx import metric_tools
from ...bfx import vawk
from ...bfx import svaba

log = logging.getLogger(__name__)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = DnaSeq.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = DnaSeq.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

    # Pipeline config
    config_files = parsed_args.config
    # Genpipes Config

    # Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    report = parsed_args.report
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    protocol = parsed_args.protocol
    design_file = parsed_args.design_file

    pipeline = DnaSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                      clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                      protocol=protocol, design_file=design_file, no_json=no_json, container=container)

    pipeline.submit_jobs()


if __name__ == '__main__':
    main()
