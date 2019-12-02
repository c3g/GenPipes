#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os
import sys
import logging

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, SanitycheckError, _raise
from core.job import Job
from bfx.readset import parse_nanopore_readset_file
from pipelines import common
from bfx import minimap2
from bfx import svim
from bfx import pycoqc

log = logging.getLogger(__name__)

class Nanopore(common.MUGQICPipeline):
    """
    Nanopore Pipeline
    ==============

    Experimental Nanopore pipeline for QC, alignment and SV calling.

    Allows for two protocols:
        - fastq : for nanopore runs where basecalling has already been done and fastq files are available (DEFAULT).
        - fast5 : for running the pipeline from raw output (requires GPU access).

    """

    def __init__(self, protocol='fastq'):
        self._protocol = protocol
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
        self.argparser.add_argument("-t", "--type", help="Type of starting file (default fastq)",
                                    choices=["fastq", "fast5"], default="fastq")
        super(Nanopore, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_nanopore_readset_file(self.args.readsets.name)
            else:
                self.argparser.error("argument -r/--readsets is required!")
        return self._readsets

    @property
    def samples(self):
        if not hasattr(self, "_samples"):
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
        return self._samples


    def guppy(self):
        """
        Use the Guppy basecaller to perform basecalling on all raw fast5 files.
        Uses the 'flip-flop' basecalling model by default.
        """
        jobs = []

        return jobs


    def blastqc(self):
        """
        Uses BLAST to perform a basic QC test by aligning 1000bp of randomly selected
        reads to the NCBI nt database in order to detect potential contamination.
        """
        jobs = []

        for readset in self.readsets:

            blast_directory = os.path.join("blastQC", readset.name)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            job = Job(
                [reads_fastq_dir],
                [os.path.join(blast_directory, readset.name + "blastHit_20MF_species.txt")],
                [["blastqc", "module_mugqic_tools"],
                 ["blastqc", "module_blast"]],
                command="""\
mkdir -p {output_directory} && \\ 
cd {output_directory} && \\ 
cat {reads_fastq_dir}/*.fastq >> full_input.tmp.fastq && \\
Nseq=$(cat full_input.tmp.fastq | awk ' {{ if (substr($0,0,1) == "+") {{ print $0}} }}' | wc -l) \\
thrC=$(echo " scale=6; 1000 / $Nseq" | bc) \\
if [ $thrC == 0 ]; then \\
    thrC=0.000001 \\
fi \\
fastqPickRandom.pl --threshold 0$thrC --input1 full_input.tmp.fastq --out1 subsample_input.fastq && \\
rm full_input.tmp.fastq && \\
python trim_nanopore.py -i subsample_input.fastq -o subsample_input.trim.fastq -s 1000 && \\
fastq2FastaQual.pl subsample_input.trim.fastq subsample_input.trim.fasta subsample_input.trim.qual && \\
blastn -query subsample_input.trim.fasta -db nt -out subsample_input.trim.blastres \\
-perc_identity 80 -num_descriptions 1 -num_alignments 1 && \\
grep ">" subsample_input.trim.blastres | awk ' {{ print $2 "_" $3}} ' | sort | uniq -c | sort -n -r | head -20 > \\
{readset_name}.blastHit_20MF_species.txt 
                """.format(
                    output_directory=blast_directory,
                    reads_fastq_dir=reads_fastq_dir,
                    readset_name=readset.name
                )
            )
            job.name = "blastqc." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs


    def minimap2_align(self):
        """
        Uses minimap2 to align the Fastq reads that passed the minimum QC threshold to
        the provided reference genome. By default, it aligns to GRCh38.
        """
        jobs = []

        sort_bam = config.param('minimap2', 'sort_output_bam')

        for readset in self.readsets:

            align_directory = os.path.join("alignments", readset.name)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            job = minimap2.minimap2_ont(reads_fastq_dir, align_directory, sort_bam)
            job.name = "minimap2_align." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs


    def pycoqc(self):
        """
        Use pycoQC to produce an interactive quality report based on the summary file and
        alignment outputs.
        """
        jobs = []

        for readset in self.readsets:

            pycoqc_directory = os.path.join("pycoQC", readset.name)

            if readset.summary_file:
                in_summary = readset.summary_file
            else:
                _raise(SanitycheckError("Error: summary file not available for readset \"" + readset.name + "\"!"))

            align_directory = os.path.join("alignments", readset.name)
            in_bam = os.path.join(align_directory, "Aligned.sortedByCoord.out.bam")

            job = pycoqc.pycoqc(readset.name, in_summary, pycoqc_directory, in_bam)
            job.name = "pycoQC." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs


    def svim(self):
        """
        Use SVIM to perform SV calling on each sample.
        """
        jobs = []

        for readset in self.readsets:

            align_directory = os.path.join("alignments", readset.name)
            in_bam = os.path.join(align_directory, "Aligned.sortedByCoord.out.bam")

            svim_directory = os.path.join("svim", readset.name)

            job = svim.svim_ont(in_bam, svim_directory)
            job.name = "svim." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs


    @property
    def steps(self):
        return [
            [self.blastqc,
             self.minimap2_align,
             self.pycoqc,
             self.svim
             ],
            [self.guppy,
             self.blastqc,
             self.minimap2_align,
             self.pycoqc,
             self.svim]
        ]


if __name__ == '__main__':
    Nanopore(protocol=['fastq', 'fast5'])
