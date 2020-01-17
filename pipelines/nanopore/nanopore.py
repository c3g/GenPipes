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
import collections
import re

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, SanitycheckError, _raise
from core.job import Job, concat_jobs
from bfx.readset import parse_nanopore_readset_file
from pipelines import common
from bfx import minimap2
from bfx import svim
from bfx import pycoqc
from bfx import gatk4

log = logging.getLogger(__name__)

class Nanopore(common.MUGQICPipeline):
    """
    Nanopore Pipeline
    ==============

    Experimental Nanopore pipeline for QC, alignment and SV calling.
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
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
                 ["blastqc", "module_blast"],
                 ["blastqc", "module_python"]],
                command="""\
mkdir -p {output_directory} && \\ 
cat {reads_fastq_dir}/*.fastq >> {output_directory}/full_input.tmp.fastq && \\
Nseq=$(cat {output_directory}/full_input.tmp.fastq | awk ' {{ if (substr($0,0,1) == "+") {{ print $0}} }}' | wc -l) && \\
thrC=$(echo " scale=6; 1000 / $Nseq" | bc) && \\
if [ $thrC == 0 ]; then thrC=0.000001; fi  && \\
fastqPickRandom.pl --threshold 0$thrC --input1 {output_directory}/full_input.tmp.fastq --out1 {output_directory}/subsample_input.fastq && \\
rm {output_directory}/full_input.tmp.fastq && \\
trim_nanopore.py -i {output_directory}/subsample_input.fastq -o {output_directory}/subsample_input.trim.fastq -s 1000 && \\
fastq2FastaQual.pl {output_directory}/subsample_input.trim.fastq {output_directory}/subsample_input.trim.fasta {output_directory}/subsample_input.trim.qual && \\
blastn -query {output_directory}/subsample_input.trim.fasta -db nt -out {output_directory}/subsample_input.trim.blastres -perc_identity 80 -num_descriptions 1 -num_alignments 1 && \\
grep ">" {output_directory}/subsample_input.trim.blastres | awk ' {{ print $2 "_" $3}} ' | sort | uniq -c | sort -n -r | head -20 > {output_directory}/{readset_name}.blastHit_20MF_species.txt 
                """.format(
                    output_directory=blast_directory,
                    reads_fastq_dir=reads_fastq_dir,
                    readset_name=readset.name
                ),
                removable_files=[os.path.join(blast_directory, "subsample_input.trim.blastres"),
                                 os.path.join(blast_directory, "subsample_input.trim.fasta"),
                                 os.path.join(blast_directory, "subsample_input.trim.fastq"),
                                 os.path.join(blast_directory, "subsample_input.trim.qual"),
                                 os.path.join(blast_directory, "subsample_input.fastq")]
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

        for readset in self.readsets:

            alignment_directory = os.path.join("alignment", readset.sample.name, readset.name)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            read_group = "'@RG" + \
                         "\\tID:" + readset.name + \
                         "\\tSM:" + readset.sample.name + \
                         "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                         ("\\tPU:run" + readset.run if readset.run else "") + \
                         "\\tPL:Nanopore" + \
                         "'"

            job = minimap2.minimap2_ont(readset.name, reads_fastq_dir, alignment_directory, read_group)
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

            align_directory = os.path.join("alignment", readset.sample.name, readset.name)
            in_bam = os.path.join(align_directory, readset.name + ".sorted.bam")

            job = pycoqc.pycoqc(readset.name, in_summary, pycoqc_directory, in_bam)
            job.name = "pycoqc." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample.
        Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:
        Aligned and sorted BAM output files from previous minimap2_align step
        """
        jobs = []

        for sample in self.samples:

            alignment_directory = os.path.join("alignment", sample.name)

            # Find input readset BAMs first from previous minimap2_align job,
            readset_bams = self.select_input_files([
                [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
                 for readset in sample.readsets]
            ])

            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            mkdir_job = Job(command="mkdir -p " + os.path.dirname(sample_bam), samples=[sample])

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:

                readset_bam = readset_bams[0]

                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, alignment_directory)

                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                target_readset_index = re.sub("\.bam$", ".bai", target_readset_bam)
                sample_index = re.sub("\.bam$", ".bai", sample_bam)

                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam],
                        [sample_bam],
                        command="ln -s -f " + target_readset_bam + " " + sample_bam,
                        removable_files=[sample_bam]),
                    Job([readset_index],
                        [sample_index],
                        command="ln -s -f " + target_readset_index + " " + sample_index + " && sleep 180",
                        removable_files=[sample_index])
                ], name="symlink_readset_sample_bam." + sample.name)
                job.samples = [sample]

            elif len(sample.readsets) > 1:

                job = concat_jobs([
                    mkdir_job,
                    gatk4.merge_sam_files(readset_bams, sample_bam)
                ])
                job.samples = [sample]
                job.name = "picard_merge_sam_files." + sample.name

            jobs.append(job)

        return jobs

    def svim(self):
        """
        Use SVIM to perform SV calling on each sample.
        """
        jobs = []

        for sample in self.samples:

            align_directory = os.path.join("alignment", sample.name)
            in_bam = os.path.join(align_directory, sample.name + ".sorted.bam")

            svim_directory = os.path.join("svim", sample.name)

            job = svim.svim_ont(in_bam, svim_directory)
            job.name = "svim." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.blastqc,
            self.minimap2_align,
            self.pycoqc,
            self.picard_merge_sam_files,
            self.svim
        ]


if __name__ == '__main__':
    Nanopore()
