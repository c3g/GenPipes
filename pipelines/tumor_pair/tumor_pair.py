<<<<<<< HEAD
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
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, pipe_jobs
from bfx.sample_tumor_pairs import parse_tumor_pair_file
from bfx.sequence_dictionary import split_by_size, parse_sequence_dictionary_file
import utils.utils

import gzip
from sys import stderr
from pipelines.dnaseq import dnaseq

#utilizes
from bfx import sambamba
from bfx import bcftools
from bfx import tools
from bfx import bed_file
from bfx import metric_tools
from bfx import bvatools
from bfx import vt
from bfx import snpeff
from bfx import vawk
from bfx import deliverables
from bfx import bash_cmd as bash

#metrics
from bfx import conpair
from bfx import multiqc

#variants
from bfx import htslib
from bfx import samtools
from bfx import varscan
from bfx import gatk
from bfx import gatk4
from bfx import vardict
from bfx import strelka2
from bfx import bcbio_variation_recall
from bfx import gemini

#sv
from bfx import delly
from bfx import manta
from bfx import lumpy
from bfx import svtyper
from bfx import wham
from bfx import metasv
from bfx import cnvkit
from bfx import scones
from bfx import sequenza
from bfx import shapeit
from bfx import scnaphase
from bfx import svaba
from bfx import annotations

log = logging.getLogger(__name__)

class TumorPair(dnaseq.DnaSeqRaw):
    """
    Tumor Pair Pipeline
    =================

    The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline with the exception of the
    indel realignment (IR) step. In the tumor pipeline the IR step utilizes both the normal and tumor bam to further reduce
    false positives (FPs) in and around indels. The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step. 
    At this point, a paired caller is used to call SNVs and Indels from the pairs given as input. Additional, muliple cancer callers 
    are utilized using an ensemble approach and SNVs and Indels seen in at least 2 different callers are retained for further 
    investigation.

    Example command:
    python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv
    
    -c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin

    -r readset: derived from GQ lims or made yourself. See : https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-readset-file

    -p pairs : format - patient_name,normal_sample_name,tumor_sample_name 
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)
        self.argparser.add_argument("--profyle", help="adjust deliverables to PROFYLE folder conventions (Default: False)", action="store_true")
        self.argparser.add_argument("-t", "--type", help="Tumor pair analysis type",choices = ["fastpass", "ensemble", "sv"], default="ensemble")
        super(TumorPair, self).__init__(protocol)


    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(
                self.args.pairs.name,
                self.samples,
                self.args.profyle
            )
        return self._tumor_pairs

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(
                config.param('DEFAULT', 'genome_dictionary', type='filepath'),
                variant=True
            )
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary_variant()):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary_variant()]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))
            windows = []

            for sequence in self.sequence_dictionary_variant():
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))

        return windows

    def is_gz_file(self, name):
        if not os.path.isfile(name):
            return True
        #if os.stat(name).st_size == 0:
        #    return False

        with gzip.open(name, 'rb') as f:
            try:
                file_content = f.read(1)
                return len(file_content) > 0
            except:
                return False

    def build_ped_file(self, directory, tumor_pair):
        ped_file = os.path.join(directory, tumor_pair.name + ".ped")
        ped_job = Job(
            command="""\
`cat > {ped_file} << END
#Family_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\tPhenotype\tEthnicity
1\t{normal}\t-9\t-9\t0\t1\t-9
1\t{tumor}\t-9\t-9\t0\t2\t-9
END`""".format(
            ped_file=ped_file,
            normal=tumor_pair.normal.name,
            tumor=tumor_pair.tumor.name,
            )
        )

        return ped_job

    def sym_link_fastq_pair(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Normal"] = [self.select_input_files([
                [readset.fastq1], [os.path.join("raw_reads", readset.sample.name, readset.name + ".pair1.fastq.gz")]])
                                for readset in tumor_pair.readsets[tumor_pair.normal.name]][0]
            inputs["Normal"].append([self.select_input_files([
                [readset.fastq2], [os.path.join("raw_reads", readset.sample.name, readset.name + ".pair2.fastq.gz")]])
                                for readset in tumor_pair.readsets[tumor_pair.normal.name]][0][0])

            inputs["Tumor"] = [self.select_input_files([
                [readset.fastq1], [os.path.join("raw_reads", readset.sample.name, readset.name + ".pair1.fastq.gz")]])
                                for readset in tumor_pair.readsets[tumor_pair.tumor.name]][0]
            inputs["Tumor"].append([self.select_input_files([
                [readset.fastq2], [os.path.join("raw_reads", readset.sample.name, readset.name + ".pair2.fastq.gz")]])
                                for readset in tumor_pair.readsets[tumor_pair.tumor.name]][0][0])
            
            for key,input in inputs.iteritems():
                for readset in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            readset,
                            readset + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            readset,
                            tumor_pair,
                            self.output_dir,
                            type="raw_reads",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            readset + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="raw_reads",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_fastq.pairs." + tumor_pair.name + "." + key))

        return jobs

    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.

        Note: modified to use both normal and tumor bams to reduce FPs around indels

        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            normal_alignment_directory = os.path.join("alignment", tumor_pair.normal.name)
            tumor_alignment_directory = os.path.join("alignment", tumor_pair.tumor.name)
            pair_directory = os.path.join("alignment", "realign", tumor_pair.name)

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.abspath(os.path.join(pair_directory, "all.intervals"))
                bam_postfix = ".realigned.all.bam"
                
                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.all.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory,
                                                 tumor_pair.normal.name + ".sorted.realigned.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                
                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.all.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory,
                                                tumor_pair.tumor.name + ".sorted.realigned.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    bash.mkdir(
                        pair_directory,
                        remove=True
                    ),
                    bash.chgdir(
                        pair_directory
                    ),
                    gatk.realigner_target_creator(
                        input_normal,
                        realign_intervals,
                        output_dir=self.output_dir,
                        input2=input_tumor
                    ),
                    gatk.indel_realigner(
                        input_normal,
                        input2=input_tumor,
                        output_dir=self.output_dir,
                        output_norm_dep=[normal_bam,normal_index],
                        output_tum_dep=[tumor_bam,tumor_index],
                        target_intervals=realign_intervals,
                        optional=bam_postfix
                    ),
                    # Move sample realign
                    bash.ln(
                        normal_bam,
                        normal_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        normal_index,
                        normal_output_index,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_bam,
                        tumor_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_index,
                        tumor_output_index,
                        self.output_dir
                    ),
                ], name="gatk_indel_realigner." + tumor_pair.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
                normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
                tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")
                
                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(pair_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory,
                                                     tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory,
                                                    tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            pair_directory,
                            remove=True
                        ),
                        bash.mkdir(
                            normal_realign_directory,
                            remove=True
                        ),
                        bash.mkdir(
                            tumor_realign_directory,
                            remove=True
                        ),
                        bash.chgdir(
                            pair_directory
                        ),
                        gatk.realigner_target_creator(
                            input_normal,
                            realign_intervals,
                            output_dir=self.output_dir,
                            input2=input_tumor,
                            intervals=intervals
                        ),
                        gatk.indel_realigner(
                            input_normal,
                            input2=input_tumor,
                            output_dir=self.output_dir,
                            output_norm_dep=[normal_bam,normal_index],
                            output_tum_dep=[tumor_bam,tumor_index],
                            target_intervals=realign_intervals,
                            intervals=intervals,
                            optional=bam_postfix
                        ),
                        bash.ln(
                            normal_bam,
                            normal_output_bam,
                            self.output_dir
                        ),
                        bash.ln(
                            normal_index,
                            normal_output_index,
                            self.output_dir
                        ),
                        bash.ln(
                            tumor_bam,
                            tumor_output_bam,
                            self.output_dir
                        ),
                        bash.ln(
                            tumor_index,
                            tumor_output_index,
                            self.output_dir
                        ),
                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(tumor_realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory,
                                                 tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory,
                                                tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    bash.mkdir(
                        pair_directory,
                        remove=True
                    ),
                    bash.mkdir(
                        normal_realign_directory,
                        remove=True
                    ),
                    bash.mkdir(
                        tumor_realign_directory,
                        remove=True
                    ),
                    bash.chgdir(
                        pair_directory
                    ),
                    gatk.realigner_target_creator(
                        input_normal,
                        realign_intervals,
                        output_dir=self.output_dir,
                        input2=input_tumor,
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                    gatk.indel_realigner(
                        input_normal,
                        input2=input_tumor,
                        output_dir=self.output_dir,
                        output_norm_dep=[normal_bam, normal_index],
                        output_tum_dep=[tumor_bam, tumor_index],
                        target_intervals=realign_intervals,
                        exclude_intervals=unique_sequences_per_job_others,
                        optional=bam_postfix
                    ),
                    bash.ln(
                        normal_bam,
                        normal_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        normal_index,
                        normal_output_index,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_bam,
                        tumor_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_index,
                        tumor_output_index,
                        self.output_dir
                    ),
                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))

        return jobs

    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            output = alignment_file_prefix + "sorted.dup.bam"
            [input] = self.select_input_files([
                [alignment_file_prefix + "sorted.matefixed.bam"],
                [alignment_file_prefix + "sorted.realigned.bam"],
                [alignment_file_prefix + "sorted.bam"]
            ])

            job = sambamba.markdup(input, output, os.path.join("alignment", sample.name))
            job.name = "sambamba_mark_duplicates." + sample.name
            job.samples=[sample]
            jobs.append(job)

        return jobs

    def sym_link_final_bam(self):
        """
        Create sym link of final bam for delivery of data to clients
        :return:
        """
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Normal"] = [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal")]
            inputs["Tumor"] =  [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal")]

            for key, input in inputs.iteritems():
                for sample_bam in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample_bam + ".bam",
                            sample_bam + ".bam.md5",
                            self.output_dir
                        ),
                        deliverables.md5sum(
                            sample_bam + ".bai",
                            sample_bam + ".bai.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample_bam + ".bam",
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_bam + ".bai",
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_bam + ".bam.md5",
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_bam + ".bai.md5",
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_final_bam.pairs." + tumor_pair.name + "." + key))

        return jobs

    def conpair_concordance_contamination(self):
        """
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            metrics_directory = os.path.join("metrics")
            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")
            pileup_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".gatkPileup")
            pileup_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".gatkPileup")

            concordance_out = os.path.join(metrics_directory, tumor_pair.name + ".concordance.tsv")
            contamination_out = os.path.join(metrics_directory, tumor_pair.name + ".contamination.tsv")

            jobs.append(concat_jobs([
                conpair.pileup(
                    input_normal,
                    pileup_normal
                ),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.normal.name))

            jobs.append(concat_jobs([
                conpair.pileup(
                    input_tumor,
                    pileup_tumor
                ),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.tumor.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    metrics_directory,
                    remove=False
                ),
                conpair.concordance(
                    pileup_normal,
                    pileup_tumor,
                    concordance_out
                ),
                conpair.contamination(
                    pileup_normal,
                    pileup_tumor,
                    contamination_out
                )
            ], name="conpair_concordance_contamination." + tumor_pair.name))

        return jobs

    def rawmpileup_panel(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bedfile = config.param('rawmpileup_panel', 'panel')

            for sequence in self.sequence_dictionary_variant():
                if sequence['type'] is 'primary':
                    pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                    jobs.append(concat_jobs([
                        bash.mkdir(
                            varscan_directory,
                            remove=True
                        ),
                        samtools.mpileup(
                            [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam"),
                             os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                            pair_output,
                            config.param('rawmpileup_panel', 'mpileup_other_options'),
                            region=sequence['name'],
                            regionFile=bedfile
                        ),
                        ], name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def paired_varscan2_panel(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                if sequence['type'] is 'primary':
                    input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                    output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                    output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                    output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                    output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")

                    jobs.append(concat_jobs([
                        bash.mkdir(
                            varscan_directory,
                            remove=True
                        ),
                        varscan.somatic(
                            input_pair,
                            output,
                            config.param('varscan2_somatic_panel', 'other_options'),
                            output_vcf_dep=output_vcf_gz,
                            output_snp_dep=output_snp,
                            output_indel_dep=output_indel
                        ),
                        htslib.bgzip_tabix(
                            output_snp,
                            os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz")
                        ),
                        htslib.bgzip_tabix(
                            output_indel,
                            os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")
                        ),
                        pipe_jobs([
                            bcftools.concat(
                                [os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"),
                                 os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")],
                                None
                            ),
                            Job(
                                [None],
                                [None],
                                command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' "
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vcf_gz
                            ),
                        ]),
                    ], name="varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2_panel(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")
                          for sequence in self.sequence_dictionary_variant() if sequence['type'] is 'primary']

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            somatic_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")
            germline_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz")

            for input_vcf in all_inputs:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete panel varscan2 vcf: %s\n" % input_vcf)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(
                        all_inputs,
                        None),
                    tools.fix_varscan_output(
                        None,
                        None
                    ),
                    Job(
                        [None],
                        [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                    ),
                    Job([None],
                        [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                    Job(
                        [None],
                        [None],
                        command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                    ),
                    htslib.bgzip_tabix(
                        None,
                        all_output
                    ),
                ]),
                bcftools.view(
                    all_output,
                    somatic_output,
                    config.param('merge_varscan2', 'somatic_filter_options')
                ),
                htslib.tabix(
                    somatic_output,
                    config.param('merge_varscan2', 'tabix_options', required=False)
                ),
                bcftools.view(
                    all_output,
                    germline_output,
                    config.param('merge_varscan2', 'germline_filter_options')
                ),
                htslib.tabix(
                    germline_output,
                    config.param('merge_varscan2', 'tabix_options', required=False)
                ),
            ], name="merge_varscan2." + tumor_pair.name))

        return jobs

    def preprocess_vcf_panel(self):
        """
        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into gemini
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")

            prefix = os.path.join(pair_directory, tumor_pair.name)
            output_somatic = prefix + ".varscan2.somatic.vt.vcf.gz"

            output_germline = prefix + ".varscan2.germline.vt.vcf.gz"

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        prefix + ".varscan2.somatic.vcf.gz" ,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".prep.vt.vcf.gz"
                    ),
                ]),
                tools.preprocess_varscan(
                    prefix + ".prep.vt.vcf.gz",
                    output_somatic
                ),
            ], name="preprocess_vcf_panel.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        prefix + ".varscan2.germline.vcf.gz" ,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".germline.prep.vt.vcf.gz"
                    ),
                ]),
                tools.preprocess_varscan(
                    prefix + ".germline.prep.vt.vcf.gz",
                    output_germline
                ),
            ], name="preprocess_vcf_panel.germline." + tumor_pair.name))

        return jobs

    def snp_effect_panel(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")
            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")

            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf")
            output_germline_gz = os.path.join(pair_directory,
                                              tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(
                    input_somatic,
                    output_somatic,
                    cancer_sample_file=cancer_pair_filename,
                    options=config.param('compute_cancer_effects_somatic', 'options')
                ),
                htslib.bgzip_tabix(
                    output_somatic,
                    output_somatic_gz
                ),
            ], name = "compute_cancer_effects_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(
                    input_germline,
                    output_germline,
                    cancer_sample_file=cancer_pair_filename,
                    options=config.param('compute_cancer_effects_germline', 'options')
                ),
                htslib.bgzip_tabix(
                    output_germline,
                    output_germline_gz
                ),
            ], name = "compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def gemini_annotations_panel(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            temp_dir = os.path.join(os.getcwd(), pair_directory)
            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(
                    gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                    gemini_prefix + ".somatic.gemini.db", temp_dir
                )
            ], name="gemini_annotations.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(
                    gemini_prefix + ".varscan2.germline.vt.snpeff.vcf.gz",
                    gemini_prefix + ".germline.gemini.db",
                    temp_dir
                )
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs

    def sym_link_panel(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] =  [os.path.join("pairedVariants", tumor_pair.name, "panel", tumor_pair.name)]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(
                            sample + ".varscan2.vcf.gz",
                            tumor_pair, self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".varscan2.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".varscan2.somatic.vt.snpeff.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".varscan2.somatic.vt.snpeff.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle),
                        deliverables.sym_link_pair(
                            sample + ".varscan2.germline.vt.snpeff.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".varscan2.germline.vt.snpeff.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".somatic.gemini.db",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".germline.gemini.db",
                            tumor_pair, self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_panel." + tumor_pair.name + "." + key))

        return jobs

    def run_pair_multiqc(self):

        jobs = []

        metrics_directory = os.path.join("metrics", "dna")
        input_dep = []
        for tumor_pair in self.tumor_pairs.itervalues():
            normal_directory = os.path.join(metrics_directory, tumor_pair.normal.name)
            input_normal_oxog = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name + ".oxog_metrics.txt")
            input_normal_qcbias = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name +".qcbias_metrics.txt")
            input_normal_all_picard = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name + ".all.metrics.quality_distribution.pdf")
            input_normal_qualimap = os.path.join(normal_directory, "qualimap", tumor_pair.normal.name, "genome_results.txt")

            [input_normal_fastqc] = self.select_input_files([
                [os.path.join(normal_directory, "fastqc", tumor_pair.normal.name + ".sorted.dup_fastqc.zip")],
                [os.path.join(normal_directory, "fastqc", tumor_pair.normal.name + "_fastqc.zip")],
            ])
            #input_normal_flagstat = os.path.join(normal_directory, "flagstat", tumor_pair.normal.name + ".flagstat")

            tumor_directory = os.path.join(metrics_directory, tumor_pair.tumor.name)
            input_tumor_oxog = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".oxog_metrics.txt")
            input_tumor_qcbias = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".qcbias_metrics.txt")
            input_tumor_all_picard = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".all.metrics.quality_distribution.pdf")
            input_tumor_qualimap = os.path.join(tumor_directory, "qualimap", tumor_pair.tumor.name, "genome_results.txt")

            [input_tumor_fastqc] = self.select_input_files([
                [os.path.join(tumor_directory, "fastqc", tumor_pair.tumor.name + ".sorted.dup_fastqc.zip")],
                [os.path.join(tumor_directory, "fastqc", tumor_pair.tumor.name + "_fastqc.zip")],
            ])
            #input_tumor_flagstat = os.path.join(tumor_directory, "flagstat", tumor_pair.tumor.name + ".flagstat")

            input = [
                os.path.join(metrics_directory, tumor_pair.normal.name),
                os.path.join(metrics_directory,tumor_pair.tumor.name)
            ]

            input_dep += [
                input_normal_oxog,
                input_normal_qcbias,
                input_normal_all_picard,
                input_normal_qualimap,
                input_normal_fastqc,
                input_tumor_oxog,
                input_tumor_qcbias,
                input_tumor_all_picard,
                input_tumor_qualimap,
                input_tumor_fastqc
            ]

            output = os.path.join(metrics_directory, tumor_pair.name + ".multiqc")

            jobs.append(
                concat_jobs([
                    multiqc.run(
                        input_dep,
                        output
                        )
            ], name="multiqc." + tumor_pair.name))

        return jobs

    def sym_link_report(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] = [os.path.join("metrics", "dna", tumor_pair.name + ".multiqc.html")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(
                            sample, tumor_pair,
                            self.output_dir,
                            type="metrics",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_fastq.report." + tumor_pair.name + "." + key))

        return jobs

    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed

            input_normal = self.select_input_files([[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                                                    [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                                                    [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files([[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                                                   [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                                                   [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            nb_jobs = config.param('rawmpileup', 'nb_jobs', type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
            
            if nb_jobs == 1:
                pair_output = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
                jobs.append(concat_jobs([
                    bash.mkdir(
                        varscan_directory,
                        remove=True
                    ),
                    samtools.mpileup(
                        [input_normal[0], input_tumor[0]],
                        pair_output,
                        config.param('rawmpileup', 'mpileup_other_options'),
                        regionFile=bed_file
                    ),
                ], name="rawmpileup." + tumor_pair.name))

            else:
                
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] is 'primary':
                        pair_output = os.path.join(varscan_directory,
                                                   tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        jobs.append(concat_jobs([
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            samtools.mpileup(
                                [input_normal[0], input_tumor[0]],
                                pair_output,
                                config.param('rawmpileup', 'mpileup_other_options'),
                                region=sequence['name'],
                                regionFile=bed_file
                            ),
                        ], name="rawmpileup." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data. 
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collison with Samtools output during ensemble                     
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
            output = os.path.join(varscan_directory, tumor_pair.name)

            nb_jobs = config.param('rawmpileup', 'nb_jobs', type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
    
                output_snp = os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")
    
                jobs.append(concat_jobs([
                    bash.mkdir(
                        varscan_directory,
                        remove=True
                    ),
                    varscan.somatic(
                        input_pair,
                        output,
                        config.param('varscan2_somatic', 'other_options'),
                        output_vcf_dep=output_vcf,
                        output_snp_dep=output_snp,
                        output_indel_dep=output_indel
                    ),
                    htslib.bgzip_tabix(
                        output_snp,
                        os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz")
                    ),
                    htslib.bgzip_tabix(
                        output_indel,
                        os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                    ),
                    pipe_jobs([
                        bcftools.concat(
                            [os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz"),
                             os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")],
                            None
                        ),
                        Job(
                            [None],
                            [output_vcf],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                    + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > "
                                    + output_vcf
                        ),
                    ]),
                    htslib.bgzip_tabix(
                        output_vcf,
                        output_vcf_gz
                    ),
                ], name="varscan2_somatic." + tumor_pair.name))

            else:

                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] is 'primary':
                        input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                        output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                        output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                        output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")

                        jobs.append(concat_jobs([
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                config.param('varscan2_somatic', 'other_options'),
                                output_vcf_dep=output_vcf,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")
                            ),
                            pipe_jobs([
                                bcftools.concat(
                                    [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"),
                                     os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")],
                                    None
                                ),
                                Job(
                                    [None],
                                    [output_vcf],
                                    command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                            + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > "
                                            + output_vcf
                                ),
                            ]),
                            htslib.bgzip_tabix(
                                output_vcf,
                                output_vcf_gz
                            ),
                        ], name="varscan2_somatic." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup', 'nb_jobs', type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            all_inputs = []
            if nb_jobs == 1:
                all_inputs = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")
                
            else:
                all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")
                              for sequence in self.sequence_dictionary_variant() if sequence['type'] is 'primary']

            for input_vcf in all_inputs:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete varscan2 vcf: %s\n" % input_vcf)

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")

            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.view(
                            all_inputs,
                            None
                        ),
                        tools.fix_varscan_output(
                            None,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        #vt.sort("-", all_output, "-m full"),
                        htslib.bgzip_tabix(
                            None,
                            all_output
                        ),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            all_output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            all_output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            all_output_vt,
                            None,
                            config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            somtic_output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            all_output_vt,
                            None,
                            config.param('varscan2_readcount_fpfilter', 'germline_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            germline_output_vt
                        ),
                    ]),
            	], name="merge_varscan2." + tumor_pair.name))

            else:
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            all_inputs,
                            None
                        ),
                        tools.fix_varscan_output(
                            None,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        #vt.sort("-", all_output, "-m full"),
                        htslib.bgzip_tabix(
                            None,
                            all_output
                        ),
                ]),
                #htslib.tabix(all_output),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        all_output,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        all_output_vt
                    ),
                ]),
                pipe_jobs([
                    bcftools.view(
                        all_output_vt,
                        None,
                        config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')
                    ),
                    htslib.bgzip_tabix(
                        None,
                        somtic_output_vt
                    ),
                ]),
                pipe_jobs([
                    bcftools.view(
                        all_output_vt,
                        None,
                        config.param('varscan2_readcount_fpfilter', 'germline_filter_options')
                    ),
                    htslib.bgzip_tabix(
                        None,
                        germline_output_vt
                    ),
                ]),
            ], name="merge_varscan2." + tumor_pair.name))

        return jobs

    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        created_interval_lists = []
        
        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            
            input_normal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)
    
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(
                        None,
                        coverage_bed,
                        interval_list
                    )
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    bash.mkdir(
                        mutect_directory,
                        remove=True
                    ),
                    gatk4.mutect2(
                        input_normal[0],
                        tumor_pair.normal.name,
                        input_tumor[0],
                        tumor_pair.tumor.name,
                        os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"),
                        interval_list=interval_list
                    )
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    outprefix = tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            mutect_directory,
                            remove=True
                        ),
                        gatk4.mutect2(
                            input_normal[0],
                            tumor_pair.normal.name,
                            input_tumor[0],
                            tumor_pair.tumor.name,
                            os.path.join(mutect_directory, outprefix + ".vcf.gz"),
                            intervals=sequences,
                            interval_list=interval_list
                        )
                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    bash.mkdir(
                        mutect_directory,
                        remove=True
                    ),
                    gatk4.mutect2(
                        input_normal[0],
                        tumor_pair.normal.name,
                        input_tumor[0],
                        tumor_pair.tumor.name,
                        os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"),
                        exclude_intervals=unique_sequences_per_job_others,
                        interval_list=interval_list
                    )
                ], name="gatk_mutect2." + tumor_pair.name + ".others"))

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants        
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_flt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.flt.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                input_vcf = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                jobs.append(concat_jobs([
                    Job(
                        [input_vcf],
                        [output_gz],
                        command="ln -s -f " + os.path.abspath(input_vcf) + " "
                                + os.path.abspath(output_gz), samples=[tumor_pair.normal, tumor_pair.tumor]
                    ),
		            #gatk4.filter_mutect_calls(output_gz, output_flt),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output_gz,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                    + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic_vt
                        ),
                    ]),
                ], name="symlink_mutect_vcf." + tumor_pair.name))

            elif nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    inputs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                for input_vcf in inputs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete mutect2 vcf: %s\n" % input_vcf)

                if config.param('gatk_mutect2', 'module_gatk').split("/")[2] > "4":
                    output_stats = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz.stats")
                    stats = []
                    for idx, sequences in enumerate(unique_sequences_per_job):
                        stats.append(
                            os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz.stats"))
                    stats.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz.stats"))

                    jobs.append(concat_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        gatk4.cat_variants(
                            inputs,
                            output_gz
                        ),
                        gatk4.merge_stats(
                            stats,
                            output_stats
                        ),
                        gatk4.filter_mutect_calls(
                            output_gz,
                            output_flt
                        ),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_flt,
                                None
                            ),
                            Job(
                                [None],
                                [None],
                                command=" grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vt_gz
                            ),
                        ]),
                        pipe_jobs([
                            bcftools.view(
                                output_vt_gz,
                                None,
                                config.param('merge_filter_mutect2', 'filter_options')
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="merge_filter_mutect2." + tumor_pair.name))
                    
                else:
                    jobs.append(concat_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        pipe_jobs([
                            bcftools.concat(
                                inputs,
                                None,
                                config.param('merge_filter_mutect2', 'bcftools_options')
                            ),
                            Job(
                                [None],
                                [None],
                                command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                        + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                            ),
    
                            htslib.bgzip_tabix(
                                None,
                                output_gz
                            ),
                        ]),
                        #gatk4.filter_mutect_calls(output_gz, output_flt),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_gz,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vt_gz
                            ),
                        ]),
                        pipe_jobs([
                            bcftools.view(
                                output_vt_gz,
                                None,
                                config.param('merge_filter_mutect2', 'filter_options')
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="merge_filter_mutect2." + tumor_pair.name))

        return jobs

    def strelka2_paired_somatic(self):
        """

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            strelka2_directory = os.path.abspath(os.path.join(pair_directory, "rawStrelka2"))
            output_prefix = os.path.abspath(os.path.join(pair_directory, tumor_pair.name))

            input_normal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            mantaIndels = None
            if os.path.isfile(os.path.join("SVariants", tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")):
                mantaIndels = os.path.join("SVariants", tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if os.path.isdir(strelka2_directory):
                jobs.append(concat_jobs([
                    bash.rm(
                        strelka2_directory
                    )
                ], name="rm_strelka2_directory." + tumor_pair.name))

            if os.path.isdir(strelka2_directory):
                jobs.append(concat_jobs([
                    bash.rm(
                        strelka2_directory
                    )
                ], name="rm_strelka2_directory." + tumor_pair.name))

            if coverage_bed:
                bed_file = coverage_bed + ".gz"
                jobs.append(concat_jobs([
                    Job(
                        [coverage_bed],
                        [coverage_bed + ".sort"],
                        command="sort -V -k1,1 -k2,2n -k3,3n " + coverage_bed + " | sed 's#chr##g' > "
                                + coverage_bed + ".sort ; sleep 15"
                    ),
                    htslib.bgzip(
                        coverage_bed + ".sort",
                        coverage_bed + ".gz"
                    ),
                    htslib.tabix(
                        coverage_bed + ".gz",
                        "-p bed"
                    ),
                 ],name="bed_index." + tumor_pair.name))

            output_dep = [os.path.join(strelka2_directory, "results/variants/somatic.snvs.vcf.gz"),
                          os.path.join(strelka2_directory, "results/variants/somatic.indels.vcf.gz")]

            jobs.append(concat_jobs([
                strelka2.somatic_config(
                    input_normal[0],
                    input_tumor[0],
                    strelka2_directory,
                    bed_file,
                    mantaIndels
                ),
                strelka2.run(
                    strelka2_directory,
                    output_dep=output_dep
                ),
            ], name="strelka2_paired_somatic.call." + tumor_pair.name))
            
            jobs.append(concat_jobs([
                pipe_jobs([
                    bcftools.concat(
                        output_dep,
                        None
                    ),
                    Job(
                        [None],
                        [None],
                        command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name
                                + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_prefix + ".strelka2.vcf.gz"
                    ),
                ]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        output_prefix + ".strelka2.vcf.gz",
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_prefix + ".strelka2.vt.vcf.gz"
                    ),
                ]),
                tools.fix_genotypes_strelka(
                    output_prefix + ".strelka2.vt.vcf.gz",
                    output_prefix + ".strelka2.somatic.gt.vcf.gz",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name
                ),
                bcftools.view(
                    output_prefix + ".strelka2.somatic.gt.vcf.gz",
                    output_prefix + ".strelka2.somatic.vt.vcf.gz",
                    config.param('strelka2_paired_somatic', 'filter_options')
                ),
            ], name="strelka2_paired_somatic.filter." + tumor_pair.name))

        return jobs

    def samtools_paired(self):
        """
        Samtools caller for SNVs and Indels using verison 0.1.19.
        """

        jobs = []

        nb_jobs = config.param('samtools_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            
            input_normal = self.select_input_files([[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                                                    [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                                                    [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files([[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                                                   [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                                                   [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            bed_file = ""
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    bash.mkdir(
                        samtools_directory,
                        remove=True
                    ),
                    pipe_jobs([
                        bcftools.mpileup(
                            [input_normal[0], input_tumor[0]],
                            None,
                            options=config.param('samtools_paired', 'mpileup_other_options'),
                            regionFile=bed_file
                        ),
                        bcftools.call(
                            "",
                            os.path.join(samtools_directory, tumor_pair.name + ".bcf"),
                            options=config.param('samtools_paired', 'bcftools_calls_options')
                        ),
                    ]),
                    bcftools.index(
                        os.path.join(
                            samtools_directory,
                            tumor_pair.name + ".bcf")
                    ),
                ], name="samtools_paired." + tumor_pair.name))

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] is 'primary':
                        jobs.append(concat_jobs([
                            bash.mkdir(
                                samtools_directory,
                                remove=True
                            ),
                            pipe_jobs([
                                bcftools.mpileup(
                                    [input_normal[0], input_tumor[0]],
                                    None,
                                    options=config.param('samtools_paired', 'mpileup_other_options'),
                                    regions=sequence['name']
                                ),
                                bcftools.call(
                                    "",
                                    os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf"),
                                    config.param('samtools_paired', 'bcftools_calls_options')
                                ),
                            ]),
                            bcftools.index(
                                os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf")
                            ),
                        ], name="samtools_paired." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_filter_paired_samtools(self):
        """
        bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
        The output of bcftools is fed to varfilter, which does an additional filtering of the variants
        and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
        for all samples in the experiment.
        Additional somatic filters are performed to reduce the number of FPs: 
        1. vcflibs vcfsamplediff tags each variant with <tag>={germline,somatic,loh} to specify the type 
        of variant given the genotype difference between the two samples.
        2. bcftools filter is used to retain only variants with CLR>=15 and have STATUS=somatic from 
        vcfsamplediff
        3. bcftools filter is used to retain only variants that have STATUS=germline or STATUS=loh from
        vcfsamplediff
        """

        jobs = []
        nb_jobs = config.param('samtools_paired', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory, tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".samtools.vcf.gz")
            output_vcf_vt = os.path.join(pair_directory, tumor_pair.name + ".samtools.vt.vcf.gz")
            output_somatics = os.path.join(pair_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".samtools.germline.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory, tumor_pair.name + ".bcf")
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.view(
                            inputs,
                            None
                        ),
                        #vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job([None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                            ),
                        Job([None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                            ),
                        Job([None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                            ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf
                        ),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output_vcf,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf_vt
                        ),
                    ]),
                    pipe_jobs([
                        vawk.paired_somatic(
                            output_vcf_vt,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatics
                        ),
                    ]),
                    pipe_jobs([
                        vawk.paired_germline(
                            output_vcf_vt,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline
                        ),
                    ]),
                ], name="merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf") for sequence in self.sequence_dictionary_variant() if sequence['type'] is 'primary']

                for input_vcf in inputs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete samtools vcf: %s\n" % input_vcf)

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    bcftools.concat(
                        inputs,
                        output,
                        config.param('merge_filter_paired_samtools', 'concat_options')
                    ),
                    pipe_jobs([
                        bcftools.view(
                            output,
                            None
                        ),
                        #vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf
                        ),
                    ]),
                    vt.decompose_and_normalize_mnps(output_vcf, output_vcf_vt),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output_vcf,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf_vt
                        ),
                    ]),
                    pipe_jobs([
                        vawk.paired_somatic(
                            output_vcf_vt,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatics
                        ),
                    ]),
                    pipe_jobs([
                        vawk.paired_germline(
                            output_vcf_vt,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline
                        ),
                    ]),
                ], name="merge_filter_paired_samtools." + tumor_pair.name))

        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is AUPAC nomenclature 
        """

        ##TO DO - the BED system needs to be revisted !! 
        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')

        bed_file_list = []
        if use_bed:
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(
                bed
            )
            last_bed_file = 'vardict.tmp.' + str(nb_jobs - 1) + '.bed'
            if not os.path.exists(last_bed_file):
                bed_file_list = bed_file.split_by_size(
                    bed_intervals,
                    interval_size,
                    nb_jobs,
                    output="./vardict.tmp")
            else:
                for idx in range(nb_jobs):
                    bed_file_list.append(os.path.join("vardict.tmp." + str(idx) + ".bed"))

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            
            input_normal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            if use_bed:
                idx = 0
                for bf in bed_file_list:
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            vardict_directory,
                            remove=True
                        ),
                        pipe_jobs([
                            vardict.paired_java(
                                input_normal[0],
                                input_tumor[0],
                                tumor_pair.name,
                                None,
                                bf
                            ),
                            vardict.testsomatic(
                                None,
                                None
                            ),
                            vardict.var2vcf(
                                None,
                                tumor_pair.normal.name,
                                tumor_pair.tumor.name,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                os.path.abspath(output)
                            ),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
                if nb_jobs == 1:
                    bedjob = vardict.dict2beds(genome_dictionary, beds)
                    output = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            vardict_directory,
                            remove=True
                        ),
                        bedjob,
                        pipe_jobs([
                            vardict.paired_java(
                                input_normal[0],
                                input_tumor[0],
                                tumor_pair.name,
                                None,
                                beds.pop()
                            ),
                            vardict.testsomatic(
                                None,
                                None
                            ),
                            vardict.var2vcf(
                                None,
                                tumor_pair.normal.name,
                                tumor_pair.tumor.name,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                os.path.abspath(output)
                            ),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + ".0"))
                    
                else:
                    bedjob = vardict.dict2beds(
                        genome_dictionary,
                        beds
                    )
                    jobs.append(concat_jobs([bedjob], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                        jobs.append(concat_jobs([
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            pipe_jobs([
                                vardict.paired_java(
                                    input_normal[0],
                                    input_tumor[0],
                                    tumor_pair.name,
                                    None,
                                    beds[idx]
                                ),
                                vardict.testsomatic(
                                    None,
                                    None
                                ),
                                vardict.var2vcf(
                                    None,
                                    tumor_pair.normal.name,
                                    tumor_pair.tumor.name,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output
                                ),
                            ]),
                        ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))
        return jobs

    def merge_filter_paired_vardict(self):
        """
        The fully merged vcf is filtered using following steps:
        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
        2. Somatics identified in step 1 must have PASS filter
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_tmp = os.path.abspath(os.path.join(pair_directory, tumor_pair.name + ".vardict.tmp.vcf.gz"))
            output = os.path.join(pair_directory, tumor_pair.name + ".vardict.vcf.gz")
            output_vt = os.path.join(pair_directory, tumor_pair.name + ".vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job(
                        [os.path.abspath(inputs)],
                        [output_tmp],
                        command="ln -s -f " + os.path.abspath(inputs) + " " + output_tmp,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    ),
                    pipe_jobs([
                        Job(
                            [output_tmp],
                            [None],
                            command="zcat " + output_tmp
                                    + " | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"),
                        htslib.bgzip_tabix(
                            None,
                            output
                        )
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline_loh
                        ),
                    ]),
                ], name="symlink_vardict_vcf." + tumor_pair.name))
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(
                        os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))

                for input_vcf in inputVCFs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete vardict vcf: %s\n" % input_vcf)

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            inputVCFs,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output
                        ),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline_loh
                        ),
                    ]),
                ], name="merge_filter_paired_vardict." + tumor_pair.name))

        return jobs

    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")
            input_strelka2 = os.path.abspath(os.path.join(input_directory, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz"))
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            inputs_somatic = [input_mutect2, input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_somatic:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                bash.mkdir(
                    paired_ensemble_directory,
                    remove=True
                ),
                bcbio_variation_recall.ensemble(
                    inputs_somatic,
                    output_ensemble,
                    config.param('bcbio_ensemble_somatic', 'options')
                ),
            ], name="bcbio_ensemble_somatic." + tumor_pair.name))

        return jobs

    def ensemble_germline_loh(self):
        """
        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.germline.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")
            inputs_germline = [input_vardict, input_varscan2, input_samtools]

            for input_vcf in inputs_germline:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory,
                                           tumor_pair.name + ".ensemble.germline.vt.vcf.gz")

            if os.path.isdir(os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt-work")):
                rm_job = bash.rm(
                    os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt-work")
                )
                jobs.append(rm_job)

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                bash.mkdir(
                    paired_ensemble_directory,
                    remove=True
                ),
                bcbio_variation_recall.ensemble(
                    inputs_germline,
                    output_ensemble,
                    config.param('bcbio_ensemble_germline', 'options')
                ),
            ], name="bcbio_ensemble_germline." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            annot_directory = os.path.join("pairedVariants", "ensemble", tumor_pair.name, "rawAnnotation")
            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
    
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_somatic_variants,
                        output_somatic_variants
                    ),
                ], name="gatk_variant_annotator.somatic." + tumor_pair.name))
                
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) + ".vcf.gz")

                    jobs.append(concat_jobs([
                        bash.mkdir(
                            annot_directory,
                            remove=True
                        ),
                        gatk.variant_annotator(
                            input_normal,
                            input_tumor,
                            input_somatic_variants,
                            output_somatic_variants,
                            intervals=sequences
                        ),
                    ], name="gatk_variant_annotator.somatic." + str(idx) + "." + tumor_pair.name))

                output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz")

                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_somatic_variants,
                        output_somatic_variants,
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                ], name="gatk_variant_annotator.somatic.others." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            annot_directory = os.path.join("pairedVariants", "ensemble", tumor_pair.name, "rawAnnotation")
            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.vcf.gz")
    
            if nb_jobs == 1:
                output_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
        
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_germline_variants,
                        output_germline_variants
                    ),
                ], name="gatk_variant_annotator.germline." + tumor_pair.name))
    
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")
            
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            annot_directory,
                            remove=True
                        ),
                        gatk.variant_annotator(
                            input_normal,
                            input_tumor,
                            input_germline_variants,
                            output_germline_variants,
                            intervals=sequences
                        ),
                    ], name="gatk_variant_annotator.germline." + str(idx) + "." + tumor_pair.name))
        
                output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz")
        
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_germline_variants,
                        output_germline_variants,
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                ], name="gatk_variant_annotator.germline.others." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_somatic = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)
                vcfs_to_merge = [os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) +".vcf.gz")
                                  for idx in xrange(len(unique_sequences_per_job))]
                
                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz"))
                
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            vcfs_to_merge,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                ], name="merge_gatk_variant_annotator.somatic." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_germline = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
            
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)
                vcfs_to_merge = [os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation", tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")
                                 for idx in xrange(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz"))
        
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            vcfs_to_merge,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline
                        ),
                    ]),
                ], name="merge_gatk_variant_annotator.germline." + tumor_pair.name))

        return jobs

    def somatic_signature(self):
        """
		Extract somatic signature composition of each sample based on Alexandrov signature reference
		Analysis is done using the SomaticSignatures and the deconstructSigs R packages
		"""
    
        jobs = []
    
        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        somatic_signature_directory = os.path.join("somaticSignature")
        # input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name , tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]

        for tumor_pair in self.tumor_pairs.itervalues():
            jobs.append(concat_jobs([
                bash.mkdir(
                    somatic_signature_directory + "/tmp_vcf",
                    remove=True
                ),
                Job(
                    [os.path.join(ensemble_directory, tumor_pair.name,
                                  tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")],
                    [os.path.join(somatic_signature_directory, "tmp_vcf", tumor_pair.name + "list_tmp.tsv")],
                    command="echo -e '" + tumor_pair.tumor.name + "\t" + os.path.join(ensemble_directory, tumor_pair.name,
                                  tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") + "' > " + os.path.join(
                        somatic_signature_directory, "tmp_vcf", tumor_pair.name + "list_tmp.tsv"))
            ], name="somatic_signature.gzip." + tumor_pair.name)
            )
    
        input_list = [os.path.join(somatic_signature_directory, "tmp_vcf", tumor_pair.name + "list_tmp.tsv") for
                      tumor_pair in self.tumor_pairs.itervalues()]
    
        jobs.append(concat_jobs([
            Job(
                input_list,
                [os.path.join(somatic_signature_directory, "Samplelist.tsv")],
                command="echo -e 'tumor\tpath' > " + os.path.join(somatic_signature_directory, "Samplelist.tsv")
                        + " && cat " + " ".join(input_list) + " >> " + os.path.join(somatic_signature_directory, "Samplelist.tsv"),
                removable_files=input_list
            ),
            tools.r_somatic_signature(
                os.path.join(somatic_signature_directory, "Samplelist.tsv"),
                somatic_signature_directory
            )
        ], name="somatic_signature.allPairs"))
    
        return jobs

    def compute_cancer_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        if not os.path.exists(ensemble_directory):
            os.makedirs(ensemble_directory)

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            if not os.path.exists(paired_directory):
                os.makedirs(paired_directory)

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                snpeff.compute_effects(
                    input_somatic,
                    output_somatic,
                    cancer_sample_file=cancer_pair_filename,
                                       options=config.param('compute_cancer_effects_somatic', 'options')
                ),
                htslib.bgzip_tabix(
                    output_somatic,
                    output_somatic + ".gz"
                ),
            ], name="compute_cancer_effects_somatic." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory,
                                           tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                snpeff.compute_effects(
                    input_germline,
                    output_germline,
                    options=config.param('compute_cancer_effects_germline', 'options')
                ),
                htslib.bgzip_tabix(
                    output_germline,
                    output_germline + ".gz"
                ),
            ], name="compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def ensemble_somatic_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """
    
        jobs = []
    

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.dbnsfp.vcf")
            
            jobs.append(concat_jobs([
                snpeff.snpsift_dbnsfp(
                    input_vcf,
                    output_vcf
                ),
                htslib.bgzip_tabix(
                    output_vcf,
                    output_vcf + ".gz"
                ),
            ], name="dbnsfp_annotation.somatic." + tumor_pair.name))
        # job.samples = self.samples
    
        return jobs

    def ensemble_germline_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """
    
        jobs = []
    
        ensemble_directory = os.path.join("pairedVariants", "ensemble")
    
        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory,
                                      tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.dbnsfp.vcf")
        
            jobs.append(concat_jobs([
                snpeff.snpsift_dbnsfp(
                    input_vcf,
                    output_vcf
                ),
                htslib.bgzip_tabix(
                    output_vcf,
                    output_vcf + ".gz"
                ),
            ], name="dbnsfp_annotation.germline." + tumor_pair.name))
        # job.samples = self.samples
    
        return jobs

    def sample_gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)


            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                gemini.gemini_annotations(
                    gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                    gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                    temp_dir
                )
            ], name="gemini_annotations.somatic." + tumor_pair.name))

        return jobs

    def sample_gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                gemini.gemini_annotations(
                    gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                    gemini_prefix + ".germline.gemini." + gemini_version + ".db",
                    temp_dir
                )
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs

    def sym_link_ensemble(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] =  [os.path.join("pairedVariants", "ensemble", tumor_pair.name, tumor_pair.name)]

            for key,input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample + ".ensemble.somatic.vt.annot.vcf.gz",
                            sample + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.somatic.vt.annot.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.somatic.vt.annot.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.md5sum(
                            sample + ".ensemble.germline.vt.annot.vcf.gz",
                            sample + ".ensemble.germline.vt.annot.vcf.gz.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.germline.vt.annot.vcf.gz.md5",
                            tumor_pair, self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.germline.vt.annot.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample + ".ensemble.germline.vt.annot.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_ensemble." + tumor_pair.name + "." + key))

        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [
            os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                Job(
                    [input_merged_vcfs[0]],
                    [output],
                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output
                )
            ], name="gatk_combine_variants.somatic.allPairs"))

        else:

            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                gatk.combine_variants(
                    input_merged_vcfs,
                    output
                )
            ], name="gatk_combine_variants.somatic.allPairs"))

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name,
                                          tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz") for tumor_pair in
                             self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                Job(
                    [input_merged_vcfs[0]],
                    [output],
                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output
                )
            ], name="gatk_combine_variants.germline.allPairs"))

        else:

            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                gatk.combine_variants(
                    input_merged_vcfs,
                    output
                )
            ], name="gatk_combine_variants.germline.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            vt.decompose_and_normalize_mnps(
                input,
                output
            )
        ], name="decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.annot.vcf.gz")
        output_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input_vcf, output_vcf)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            vt.decompose_and_normalize_mnps(
                input_vcf,
                output_vcf
            )
        ], name="decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        for tumor_pair in self.tumor_pairs.itervalues():
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            snpeff.compute_effects(
                input,
                output,
                cancer_sample_file=cancer_pair_filename,
                options=config.param('compute_cancer_effects_somatic', 'options')
            ),
            htslib.bgzip_tabix(
                output,
                output_gz
            ),
        ], name="compute_effects.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf.gz")

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            snpeff.compute_effects(
                input,
                output,
                options=config.param('compute_cancer_effects_germline', 'options')
            ),
            htslib.bgzip_tabix(
                output,
                output_gz
            ),
        ], name="compute_effects.germline.allPair"))

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            gemini.gemini_annotations(
                gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                gemini_prefix + ".somatic.gemini.db",
                temp_dir
            )
        ], name="gemini_annotations.somatic.allPairs"))

        return jobs

    def gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        jobs.append(concat_jobs([
            bash.mkdir(
                ensemble_directory,
                remove=True
            ),
            gemini.gemini_annotations(
                gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                gemini_prefix + ".germline.gemini.db",
                temp_dir
            )
        ], name="gemini_annotations.germline.allPairs"))

        return jobs

    def sequenza(self):
        """
        Sequenza is a novel set of tools providing a fast python script to genotype cancer samples,
        and an R package to estimate cancer cellularity, ploidy, genome wide copy number profile and infer
        for mutated alleles.

        """
        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            sequenza_directory = os.path.join(pair_directory, "sequenza")
            rawSequenza_directory = os.path.join(sequenza_directory, "rawSequenza")
            
            inputNormal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])

            inputTumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed

            for sequence in self.sequence_dictionary_variant():
                if sequence['type'] is 'primary':
                    normal_mpileup = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                    tumor_mpileup = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")
                    normal_gz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.normal.name + "." + sequence['name'] + ".mpileup.gz")
                    tumor_gz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup.gz")
                    out_seqz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + "." + sequence['name'] + ".seqz.gz")
                    binned_seqz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".binned.seqz." + sequence['name'] + ".gz")

                    if os.path.isfile(normal_mpileup) and os.path.isfile(tumor_mpileup):
                        jobs.append(concat_jobs([
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            Job(
                                [normal_mpileup],
                                [normal_gz],
                                command="gzip -cf " + normal_mpileup + " > " + normal_gz
                            ),
                            Job(
                                [tumor_mpileup],
                                [tumor_gz],
                                command="gzip -cf " + tumor_mpileup + " > " + tumor_gz
                            ),
                            pipe_jobs([
                                sequenza.seqz(
                                    normal_gz,
                                    tumor_gz,
                                    config.param('sequenza', 'gc_file'),
                                    None
                                ),
                                Job(
                                    [None],
                                    [out_seqz],
                                    command="gzip -cf > " + out_seqz
                                )
                            ]),
                            pipe_jobs([
                                sequenza.bin(
                                    out_seqz,
                                    None
                                ),
                                Job(
                                    [None],
                                    [binned_seqz],
                                    command="gzip -c > " + binned_seqz
                                ),
                            ]),
                        ], name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name))

                    else:

                        jobs.append(concat_jobs([
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            pipe_jobs([
                                samtools.mpileup(
                                    [inputNormal[0]],
                                    None,
                                    config.param('sequenza', 'mpileup_options'),
                                    sequence['name'],
                                    bed_file
                                ),
                                Job(
                                    [None],
                                    [normal_gz],
                                    command="gzip -cf > " + normal_gz
                                ),
                            ]),
                            pipe_jobs([
                                samtools.mpileup(
                                    [inputTumor[0]],
                                    None,
                                    config.param('sequenza', 'mpileup_options'),
                                    sequence['name'],
                                    bed_file
                                ),
                                Job(
                                    [None],
                                    [tumor_gz],
                                    command="gzip -cf > " + tumor_gz
                                ),
                            ]),
                        ], name="mpileup_sequenza." + sequence['name'] + "." + tumor_pair.name))

                        jobs.append(concat_jobs([
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            pipe_jobs([
                                sequenza.seqz(
                                    normal_gz,
                                    tumor_gz,
                                    config.param('sequenza', 'gc_file'),
                                    None
                                ),
                                Job(
                                    [None],
                                    [out_seqz],
                                    command="gzip -c > " + out_seqz
                                ),
                            ]),
                            pipe_jobs([
                                sequenza.bin(
                                    out_seqz,
                                    None
                                ),
                                Job(
                                    [None],
                                    [binned_seqz],
                                    command="gzip -c > " + binned_seqz
                                ),
                            ]),
                        ], name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name))

            seqz_outputs = [os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".binned.seqz." + sequence['name'] + ".gz")
                            for sequence in self.sequence_dictionary_variant() if sequence['type'] is 'primary']
            #seqz_input = seqz_outputs[0]
            merged_seqz = os.path.join(sequenza_directory, tumor_pair.name + ".binned.merged.seqz.gz")

            jobs.append(concat_jobs([
                bash.mkdir(
                    rawSequenza_directory,
                    remove=True
                ),
                Job(seqz_outputs,
                    [merged_seqz],
                    command="zcat " + " \\\n".join(seqz_outputs)
                            + " \\\n | gawk 'FNR==1 && NR==1{print;}{ if($1!=\"chromosome\" && $1!=\"MT\" && $1!=\"chrMT\" && $1!=\"chrM\") {print $0} }' | \\\n   gzip -cf > "
                            + merged_seqz
                    ),
            ], name="sequenza.merge_binned_seqz." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    rawSequenza_directory,
                    remove=True
                ),
                sequenza.main(
                    merged_seqz,
                    sequenza_directory,
                    tumor_pair.name
                ),
                #sequenza.filter(os.path.join(sequenza_directory, tumor_pair.name + "_segments.txt"), tumor_pair.name, os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt")),
                #sequenza.annotate(os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt"), os.path.join(sequenza_directory, tumor_pair.name + ".annotated"),
                #                  os.path.join(sequenza_directory, tumor_pair.name + ".tmp"))
            ], name="sequenza." + tumor_pair.name))

        return jobs

    def sym_link_sequenza(self):
        jobs = []

        inputs = dict()

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory, "sequenza", tumor_pair.name + "_chromosome_view.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_genome_view.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CN_bars.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CP_contours.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_ploidy_celularity.tsv")]
                               #os.path.join(pair_directory, "sequenza", tumor_pair.name + ".annotated.TumS.filteredSV.annotate.txt")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(
                            sample, tumor_pair,
                            self.output_dir,
                            type="sv/cnv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link.sequenza." + tumor_pair.name + "." + key))

        return jobs

    def sCNAphase(self):
        """


        """
        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            raw_directory = os.path.join(pair_directory, "sCNAphase", "rawsCNAphase")
            scnaphase_directory = os.path.join(pair_directory, "sCNAphase")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            output_normal_vcf = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vcf.gz")
            output_normal_vcf_vt = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vt.vcf.gz")
            output_normal_vcf_vt_flt = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vt.flt.vcf.gz")
            output_tumor_vcf = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vcf.gz")
            output_tumor_vcf_vt = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vt.vcf.gz")
            output_tumor_vcf_vt_flt = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vt.flt.vcf.gz")

            nb_jobs = config.param('samtools_single', 'nb_jobs', type='posint')

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    bash.mkdir(
                        raw_directory,
                        remove=True
                    ),
                    pipe_jobs([
                        samtools.mpileup(
                            [inputNormal],
                            None,
                            config.param('samtools_single', 'mpileup_other_options'),
                            ini_section="samtools_single"
                        ),
                        samtools.bcftools_view(
                            "-",
                            os.path.join(raw_directory, tumor_pair.normal.name + ".bcf"),
                            config.param('samtools_single', 'bcftools_view_options'),
                            ini_section="samtools_single"
                        ),
                    ]),
                ], name="samtools_single." + tumor_pair.normal.name))

                jobs.append(concat_jobs([
                    bash.mkdir(
                        raw_directory,
                        remove=True
                    ),
                    pipe_jobs([
                        samtools.mpileup(
                            [inputTumor],
                            None,
                            config.param('samtools_single', 'mpileup_other_options'),
                            ini_section="samtools_single"
                        ),
                        samtools.bcftools_view(
                            "-",
                            os.path.join(raw_directory, tumor_pair.tumor.name + ".bcf"),
                            config.param('samtools_single', 'bcftools_view_options'),
                            ini_section="samtools_single"
                        ),
                    ]),
                ], name="samtools_single." + tumor_pair.tumor.name))

            else:
                for region in self.generate_approximate_windows(nb_jobs):  # for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            raw_directory,
                            remove=True
                        ),
                        pipe_jobs([
                            samtools.mpileup(
                                [inputNormal],
                                None,
                                config.param('samtools_single', 'mpileup_other_options'),
                                region,
                                ini_section="samtools_single"
                            ),
                            samtools.bcftools_view(
                                "-",
                                os.path.join(raw_directory, tumor_pair.normal.name + "." + region + ".bcf"),
                                config.param('samtools_single', 'bcftools_view_options'),
                                ini_section="samtools_single"
                            ),
                        ]),
                    ], name="samtools_single." + tumor_pair.normal.name + "." + region))

                    jobs.append(concat_jobs([
                        bash.mkdir(
                            raw_directory,
                            remove=True
                        ),
                        pipe_jobs([
                            samtools.mpileup(
                                [inputTumor],
                                None,
                                config.param('samtools_single', 'mpileup_other_options'),
                                region,
                                ini_section="samtools_single"),
                            samtools.bcftools_view(
                                "-",
                                os.path.join(raw_directory, tumor_pair.tumor.name + "." + region + ".bcf"),
                                config.param('samtools_single', 'bcftools_view_options'),
                                ini_section="samtools_single"
                            ),
                        ]),
                    ], name="samtools_single." + tumor_pair.tumor.name + "." + region))

                inputsNormal = [os.path.join(raw_directory, tumor_pair.normal.name + "." + region + ".bcf") for region in
                                self.generate_approximate_windows(nb_jobs)]
                
                jobs.append(concat_jobs([
                    pipe_jobs([
                        samtools.bcftools_cat(
                            inputsNormal,
                            None,
                            ini_section="samtools_single"
                        ),
                        samtools.bcftools_view(
                            "-",
                            None,
                            ini_section="samtools_single"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_normal_vcf
                        ),
                    ]),
                    vt.decompose_and_normalize_mnps(
                        output_normal_vcf,
                        output_normal_vcf_vt
                    ),
                    pipe_jobs([
                        Job(
                            [output_normal_vcf_vt],
                            [None],
                            command="zgrep -Pv '\\tN\\t' " + output_normal_vcf_vt
                                    + " | grep -v 'INDEL' | grep -v '\.\/' | grep -v '\/\.' "
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_normal_vcf_vt_flt
                        ),
                    ]),
                ], name="merge_samtools_single." + tumor_pair.normal.name))

                inputsTumor = [os.path.join(raw_directory, tumor_pair.tumor.name + "." + region + ".bcf") for
                               region in self.generate_approximate_windows(nb_jobs)]
                jobs.append(concat_jobs([
                    pipe_jobs([
                        samtools.bcftools_cat(
                            inputsTumor,
                            None,
                            ini_section="samtools_single"
                        ),
                        samtools.bcftools_view(
                            "-",
                            None,
                            ini_section="samtools_single"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_tumor_vcf
                        ),
                    ]),
                    vt.decompose_and_normalize_mnps(
                        output_tumor_vcf,
                        output_tumor_vcf_vt
                    ),
                    pipe_jobs([
                        Job(
                            [output_tumor_vcf_vt],
                            [None],
                            command="zgrep -Pv '\\tN\\t' " + output_tumor_vcf_vt
                                    + " | grep -v 'INDEL' | grep -v '\.\/' | grep -v '\/\.' "
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_tumor_vcf_vt_flt),
                    ]),
                ], name="merge_samtools_single." + tumor_pair.tumor.name))

                shapeit_nprefix = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".")
                shapeit_tprefix = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".")

                for chr in range(1, 23):
                    jobs.append(concat_jobs([
                        htslib.tabix_split(
                            output_normal_vcf_vt_flt,
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"),
                            str(chr)
                        ),
                        htslib.tabix_split(
                            output_tumor_vcf_vt_flt,
                            os.path.join(shapeit_tprefix + "chr" + str(chr) + ".vcf"),
                            str(chr)
                        ),
                    ], name="tabix_split." + tumor_pair.name + "." + str(chr)))

                for chr in range(1, 23):
                    jobs.append(concat_jobs([
                        shapeit.check(
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"),
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".alignments"),
                            str(chr)
                        ),
                    ], name="shapeit.check." + tumor_pair.normal.name + "." + str(chr)))

                    jobs.append(concat_jobs([
                        shapeit.phase(
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"),
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".alignments.snp.strand.exclude"),
                            os.path.join(shapeit_nprefix + "chr" + str(chr)),
                            os.path.join(shapeit_nprefix + "chr" + str(chr) + ".phase"),
                            str(chr)
                        ),
                    ], name="shapeit.phase." + tumor_pair.normal.name + "." + str(chr)))

            jobs.append(concat_jobs([
                bash.mkdir(
                    raw_directory,
                    remove=True
                ),
                Job(
                    command="cd " + scnaphase_directory
                ),
                scnaphase.run(
                    tumor_pair.name,
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name
                ),
            ], name="scnaphase." + tumor_pair.name))

        return jobs

    def delly_call_filter(self):
        """
        Delly2 is an integrated structural variant prediction method that can
        discover, genotype and visualize deletions, tandem duplications, inversions and translocations
        at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends
        and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.
        Structural variants can be visualized using Delly-maze and Delly-suave.
        input: normal and tumor final bams
        Returns:bcf file

        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():

            pair_directory = os.path.join("SVariants", tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")

            filename = os.path.join(delly_directory, tumor_pair.name + '.tsv')
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))
           
            cancer_pair = open(filename, 'w')
            cancer_pair.write(tumor_pair.tumor.name + "\ttumor\n")
            cancer_pair.write(tumor_pair.normal.name + "\tcontrol\n")

            inputNormal = os.path.join("alignment", tumor_pair.normal.name,
                                       tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            inputs = [inputTumor, inputNormal]

            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            for sv_type in SV_types:
                output_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf")
                output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".somatic.flt.vcf.gz")

                jobs.append(concat_jobs([
                    bash.mkdir(
                        delly_directory,
                        remove=True
                    ),
                    delly.call(
                        inputs,
                        output_bcf,
                        sv_type
                    ),
                    pipe_jobs([
                        bcftools.view(
                            output_bcf,
                            None,
                            config.param('delly_call_filter_somatic', 'bcftools_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf
                        ),
                    ]),
                ], name="delly_call_filter." + str(sv_type) + "." + tumor_pair.name))

        return jobs

    def delly_sv_annotation(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():

            pair_directory = os.path.join("SVariants", tumor_pair.name)
            final_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")
            output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly.merge.sort.vcf.gz")
            output_flt_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.merge.sort.flt.vcf.gz")
            
            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            inputBCF = []
            for sv_type in SV_types:
                inputBCF.append(os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf"))

            jobs.append(concat_jobs([
                pipe_jobs([
                    bcftools.concat(
                        inputBCF,
                        None,
                        "-O v"
                    ),
                    vt.sort(
                        "-",
                        "-",
                        "-m full"
                    ),
                    htslib.bgzip(
                        None,
                        output_vcf
                    ),
                ]),
                pipe_jobs([
                    bcftools.view(
                        output_vcf,
                        None,
                        "-f PASS"
                    ),
                    htslib.bgzip(
                        None,
                        output_flt_vcf
                    ),
                ]),
            ], name="sv_annotation.delly.merge_sort_filter." + tumor_pair.name))

            jobs.append(concat_jobs([
                vawk.paired_somatic(
                    output_flt_vcf,
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    final_directory + ".delly.somatic.vcf"
                ),
                snpeff.compute_effects(
                    final_directory + ".delly.somatic.vcf",
                    final_directory + ".delly.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    final_directory + ".delly.somatic.snpeff.vcf",
                    final_directory + ".delly.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    final_directory + ".delly.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "DELLY",
                    final_directory + ".delly.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.delly.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                vawk.paired_germline(
                    output_flt_vcf,
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    final_directory + ".delly.germline.vcf"
                ),
                snpeff.compute_effects(
                    final_directory + ".delly.germline.vcf",
                    final_directory + ".delly.germline.snpeff.vcf"
                ),
                annotations.structural_variants(
                    final_directory + ".delly.germline.snpeff.vcf",
                    final_directory + ".delly.germline.snpeff.annot.vcf"
                ),
                vawk.sv(
                    final_directory + ".delly.germline.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "DELLY",
                    final_directory + ".delly.germline.prioritize.tsv"
                ),
            ], name="sv_annotation.delly.germline." + tumor_pair.name))
            
        return jobs
    
    def sym_link_delly(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".delly.somatic.snpeff.annot.vcf",
                               pair_directory + ".delly.somatic.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_delly.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".delly.germline.snpeff.annot.vcf",
                               pair_directory + ".delly.germline.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_delly.germline." + tumor_pair.name + "." + key))

        return jobs
    
        
    def manta_sv_calls(self):
        """
        Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
        analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
        Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
        single efficient workflow.
        Returns:Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
         in VCF 4.1 format.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            manta_directory = os.path.abspath(os.path.join(pair_directory, "rawManta"))
            output_prefix = os.path.abspath(os.path.join(pair_directory, tumor_pair.name))

            mkdir_job = Job(command="mkdir -p " + manta_directory, removable_files=[manta_directory], samples = [tumor_pair.normal, tumor_pair.tumor])

            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            manta_somatic_output = os.path.join(manta_directory, "results/variants/somaticSV.vcf.gz")
            manta_germline_output = os.path.join(manta_directory, "results/variants/diploidSV.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed + ".gz"
                jobs.append(concat_jobs([
                    Job(
                        [coverage_bed],
                        [coverage_bed + ".sort"],
                        command="sort -V -k1,1 -k2,2n -k3,3n " + coverage_bed + " | sed 's#chr##g' > "
                                + coverage_bed + ".sort"
                    ),
                    htslib.bgzip(
                        coverage_bed + ".sort",
                        coverage_bed + ".gz"
                    ),
                    htslib.tabix(
                        coverage_bed + ".gz",
                        "-p bed"
                    ),
                 ],name="bed_index." + tumor_pair.name))

            output_dep = [manta_somatic_output, manta_somatic_output + ".tbi", manta_germline_output, manta_germline_output + ".tbi"]

            jobs.append(concat_jobs([
                bash.mkdir(
                    manta_directory,
                    remove=True
                ),
                manta.manta_config(
                    inputNormal,
                    inputTumor,
                    manta_directory,
                    bed_file
                ),
                manta.manta_run(
                    manta_directory,
                    output_dep=output_dep
                ),
                bash.ln(
                    manta_somatic_output,
                    output_prefix + ".manta.somatic.vcf.gz",
                    self.output_dir,
                ),
                bash.ln(
                    manta_somatic_output + ".tbi",
                    output_prefix + ".manta.somatic.vcf.gz.tbi",
                    self.output_dir
                ),
                bash.ln(
                    manta_germline_output,
                    output_prefix + ".manta.germline.vcf.gz",
                    self.output_dir,
                ),
                bash.ln(
                    manta_germline_output + ".tbi",
                    output_prefix + ".manta.germline.vcf.gz.tbi",
                    self.output_dir,
                ),
            ], name="manta_sv." + tumor_pair.name))

        return jobs

    def manta_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    pair_directory + ".manta.somatic.vcf.gz",
                    pair_directory + ".manta.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    pair_directory + ".manta.somatic.snpeff.vcf",
                    pair_directory + ".manta.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    pair_directory + ".manta.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "MANTA",
                    pair_directory + ".manta.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.manta_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    pair_directory + ".manta.germline.vcf.gz",
                    pair_directory + ".manta.germline.snpeff.vcf"
                ),
                annotations.structural_variants(
                    pair_directory + ".manta.germline.snpeff.vcf",
                    pair_directory + ".manta.germline.snpeff.annot.vcf"
                ),
                vawk.sv(
                    pair_directory + ".manta.germline.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "MANTA",
                    pair_directory + ".manta.germline.prioritize.tsv"
                )
            ], name="sv_annotation.manta_germline." + tumor_pair.name))

        return jobs

    def sym_link_manta(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".manta.somatic.snpeff.annot.vcf"),
                               pair_directory + ".manta.somatic.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_manta.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".manta.germline.snpeff.annot.vcf"),
                               pair_directory + ".manta.germline.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_manta.germline." + tumor_pair.name + "." + key))

        return jobs

    def lumpy_paired_sv(self):
        """
        A probabilistic framework for structural variant discovery.
        Lumpy traditional with paired ends and split reads on tumor normal pair.
        Returns:bams.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            lumpy_directory = os.path.join(pair_directory, "rawLumpy")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            discordants_normal = os.path.join(lumpy_directory, tumor_pair.normal.name + ".discordants.sorted.bam")
            discordants_tumor = os.path.join(lumpy_directory, tumor_pair.tumor.name + ".discordants.sorted.bam")

            splitters_tumor = os.path.join(lumpy_directory, tumor_pair.tumor.name + ".splitters.sorted.bam")
            splitters_normal = os.path.join(lumpy_directory, tumor_pair.normal.name + ".splitters.sorted.bam")

            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf")
            gzip_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf.gz")

            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")
            genotype_gzip = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf.gz")

            jobs.append(concat_jobs([
                bash.mkdir(
                    lumpy_directory,
                    remove=True
                ),
                pipe_jobs([
                    samtools.view(
                        inputNormal,
                        None,
                        "-b -F 1294"
                    ),
                    sambamba.sort(
                        "/dev/stdin",
                        discordants_normal,
                        lumpy_directory,
                        config.param('extract_discordant_reads', 'options')
                    ),
                ]),
                pipe_jobs([
                    samtools.view(
                        inputTumor,
                        None,
                        "-b -F 1294"
                    ),
                    sambamba.sort(
                        "/dev/stdin",
                        discordants_tumor,
                        lumpy_directory,
                        config.param('extract_discordant_reads', 'options')
                    ),
                ]),
            ], name="extract_discordant_reads." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    lumpy_directory,
                    remove=True
                ),
                pipe_jobs([
                    samtools.view(
                        inputNormal,
                        None,
                        "-h"
                    ),
                    Job(
                        [None],
                        [None],
                        [['lumpy_sv', 'module_lumpy']],
                        command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"
                    ),
                    samtools.view(
                        "-",
                        None,
                        " -Sb "
                    ),
                    sambamba.sort(
                        "/dev/stdin",
                        splitters_normal,
                        lumpy_directory,
                        config.param('extract_split_reads', 'options')
                    ),
                ]),
                pipe_jobs([
                    samtools.view(
                        inputTumor,
                        None,
                        "-h"
                    ),
                    Job(
                        [None],
                        [None],
                        [['lumpy_sv', 'module_lumpy']],
                        command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"
                    ),
                    samtools.view(
                        "-",
                        None,
                        " -Sb "
                    ),
                    sambamba.sort(
                        "/dev/stdin",
                        splitters_tumor,
                        lumpy_directory,
                        config.param('extract_split_reads', 'options')
                    ),
                ]),
            ], name="extract_split_reads." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    lumpy_directory,
                    remove=True
                ),
                lumpy.lumpyexpress_pair(
                    inputNormal,
                    inputTumor,
                    output_vcf,
                    spl_normal=splitters_normal,
                    spl_tumor=splitters_tumor,
                    dis_normal=discordants_normal,
                    dis_tumor=discordants_tumor
                ),
                htslib.bgzip(
                    output_vcf,
                    gzip_vcf
                ),
            ], name="lumpy_paired_sv_calls." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    Job(
                        [gzip_vcf],
                        [None],
                        command="zcat " + gzip_vcf + " | grep -v \"^##contig\""
                    ),
                    bcftools.annotate(
                        None,
                        None,
                        config.param('lumpy_paired_sv_calls', 'header_options')
                    ),
                    vt.sort(
                        "-",
                        os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"),
                        "-m full"
                    ),
                ]),
                svtyper.genotyper(
                    inputTumor,
                    inputNormal,
                    os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"),
                    genotype_vcf
                ),
                htslib.bgzip(
                    genotype_vcf,
                    genotype_gzip
                ),
            ], name="lumpy_paired_sv_calls.genotype." + tumor_pair.name))

        return jobs

    def lumpy_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            prefix = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            
            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.germline.vcf.gz")

            jobs.append(concat_jobs([
                pipe_jobs([
                    vawk.paired_somatic(
                        genotype_vcf,
                        tumor_pair.normal.name,
                        tumor_pair.tumor.name,
                        None
                    ),
                    htslib.bgzip(
                        None,
                        somatic_vcf
                    ),
                ]),
                pipe_jobs([
                    vawk.paired_germline(
                        genotype_vcf,
                        tumor_pair.normal.name,
                        tumor_pair.tumor.name,
                        None
                    ),
                    htslib.bgzip(
                        None,
                        germline_vcf
                    ),
                ]),
            ], name="sv_annotation.lumpy.genotypes." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    somatic_vcf,
                    prefix + ".lumpy.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    prefix + ".lumpy.somatic.snpeff.vcf",
                    prefix + ".lumpy.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    prefix + ".lumpy.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "LUMPY",
                    prefix + ".lumpy.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.lumpy.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    germline_vcf,
                    prefix + ".lumpy.germline.snpeff.vcf"
                ),
                annotations.structural_variants(
                    prefix + ".lumpy.germline.snpeff.vcf",
                    prefix + ".lumpy.germline.snpeff.annot.vcf"
                ),
                vawk.sv(
                    prefix + ".lumpy.germline.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "LUMPY",
                    prefix + ".lumpy.germline.prioritize.tsv"
                ),
            ], name="sv_annotation.lumpy.germline." + tumor_pair.name))

        return jobs

    def sym_link_lumpy(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".lumpy.somatic.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".lumpy.somatic.prioritize.tsv")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_lumpy.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".lumpy.germline.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".lumpy.germline.prioritize.tsv")]
        
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_lumpy.germline." + tumor_pair.name + "." + key))

        return jobs

    def wham_call_sv(self):
        """
        Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
        calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
        to employ SVs in association testing.
        Returns:vcf.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            wham_directory = os.path.join(pair_directory, "rawWham")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            output_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.vcf")
            merge_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.merged.vcf")
            genotyped_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.germline.vcf.gz")

            jobs.append(concat_jobs([
                bash.mkdir(
                    wham_directory,
                    remove=True
                ),
                wham.call_sv(
                    inputTumor,
                    inputNormal,
                    output_vcf
                ),
                pipe_jobs([
                    wham.merge(
                        output_vcf,
                        None
                    ),
                    Job(
                        [None],
                        [merge_vcf],
                        command="sed 's/NONE/" + tumor_pair.tumor.name + "/g' | sed -e 's#\"\"#\"#g' > " + merge_vcf
                    ),
                ]),
            ], name="wham_call_sv.call_merge." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    wham_directory,
                    remove=True
                ),
                pipe_jobs([
                    Job(
                        [merge_vcf],
                        [None],
                        command="cat " + merge_vcf + " | grep -v \"^##contig\""
                    ),
                    bcftools.annotate(
                        None,
                        None,
                        config.param('wham_call_sv', 'header_options')
                    ),
                    vt.sort(
                        "-",
                        os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"),
                        "-m full"
                    ),
                ]),
                pipe_jobs([
                    svtyper.genotyper(
                        inputTumor,
                        inputNormal,
                        os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"),
                        None
                    ),
                    Job(
                        [None],
                        [None],
                        command=" sed -e 's#\"\"#\"#g' "
                    ),
                    htslib.bgzip_tabix(
                        None,
                        genotyped_vcf
                    ),
                ]),
            ], name="wham_call_sv.genotype." + tumor_pair.name))

        return jobs

    def wham_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            prefix = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            genotyped_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")

            jobs.append(concat_jobs([
                pipe_jobs([
                    vawk.paired_somatic(
                        genotyped_vcf,
                        tumor_pair.normal.name,
                        tumor_pair.tumor.name,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".wham.somatic.vcf.gz"
                    ),
                ]),
                snpeff.compute_effects(
                    prefix + ".wham.somatic.vcf.gz",
                    prefix + ".wham.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    prefix + ".wham.somatic.snpeff.vcf",
                    prefix + ".wham.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    prefix + ".wham.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "WHAM",
                    prefix + ".wham.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.wham.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    vawk.paired_germline(
                        genotyped_vcf,
                        tumor_pair.normal.name,
                        tumor_pair.tumor.name,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".wham.germline.vcf.gz"
                    ),
                ]),
                snpeff.compute_effects(
                    prefix + ".wham.germline.vcf.gz",
                    prefix + ".wham.germline.snpeff.vcf"
                ),
                annotations.structural_variants(
                    prefix + ".wham.germline.snpeff.vcf",
                    prefix + ".wham.germline.snpeff.annot.vcf"
                ),
                vawk.sv(
                    prefix + ".wham.germline.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "WHAM", prefix + ".wham.germline.prioritize.tsv"
                ),
            ], name="sv_annotation.wham.germline." + tumor_pair.name))

        return jobs

    def sym_link_wham(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".wham.somatic.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".wham.somatic.prioritize.tsv")]
            
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_wham.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".wham.germline.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".wham.germline.prioritize.tsv")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_wham.germline." + tumor_pair.name + "." + key))

        return jobs

    def cnvkit_batch(self):
        """
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            tarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.targetcoverage.cnn")
            antitarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.antitargetcoverage.cnn")
            ref_cnn = os.path.join(cnvkit_dir, tumor_pair.name + ".reference.cnn")
            tumor_cns = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".cns")
            vcf_gz = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            metrics = os.path.join("SVariants", "cnvkit_reference")
            poolRef = os.path.join(metrics, "pooledReference.cnn")

            if os.path.isfile(poolRef):
                pool_ref_cnn = poolRef
                ref_cnn = None

            else:
                pool_ref_cnn = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            bed = None

            if coverage_bed:
                bed = coverage_bed

            vardict_vcf = os.path.join("pairedVariants", tumor_pair.name,
                                       tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            input_vcf = None
            normal = None
            tumor = None
            if os.path.isfile(vardict_vcf):
                input_vcf = vardict_vcf
                normal = tumor_pair.normal.name
                tumor = tumor_pair.tumor.name

            jobs.append(concat_jobs([
                bash.mkdir(
                    cnvkit_dir,
                    remove=True
                ),
                cnvkit.batch(
                    inputTumor,
                    inputNormal,
                    cnvkit_dir,
                    tar_dep=tarcov_cnn,
                    antitar_dep=antitarcov_cnn,
                    target_bed=bed,
                    reference=pool_ref_cnn,
                    output_cnn=ref_cnn
                ),
            ], name="cnvkit_batch." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    cnvkit_dir,
                    remove=True
                ),
                cnvkit.fix(
                    tarcov_cnn,
                    antitarcov_cnn,
                    os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                    reference=pool_ref_cnn,
                    ref_cnn=ref_cnn
                ),
                cnvkit.segment(
                    os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                    tumor_cns
                ),
            ], name="cnvkit_batch.correction." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    cnvkit_dir,
                    remove=True
                ),
                cnvkit.call(
                    tumor_cns,
                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns")
                ),
                pipe_jobs([
                    cnvkit.export(
                        os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                        None,
                        sample_id=tumor_pair.tumor.name
                    ),
                    htslib.bgzip_tabix(
                        None,
                        vcf_gz
                    ),
                ]),
            ], name="cnvkit_batch.call." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    cnvkit_dir,
                    remove=True
                ),
                cnvkit.metrics(
                    os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                    os.path.join(metrics, tumor_pair.name + ".metrics.tsv")
                ),
                cnvkit.scatter(
                    os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                    os.path.join(cnvkit_dir, tumor_pair.name + ".scatter.pdf"),
                    input_vcf,
                    normal,
                    tumor
                ),
                cnvkit.diagram(
                    os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                    os.path.join(cnvkit_dir, tumor_pair.name + ".diagram.pdf")
                ),
            ], name="cnvkit_batch.metrics." + tumor_pair.name))

        return jobs

    def cnvkit_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    pair_directory + ".cnvkit.vcf.gz",
                    pair_directory + ".cnvkit.snpeff.vcf"
                ),
                annotations.structural_variants(
                    pair_directory + ".cnvkit.snpeff.vcf",
                    pair_directory + ".cnvkit.snpeff.annot.vcf"
                ),
            ], name="sv_annotation.cnvkit." + tumor_pair.name))

        return jobs

    def sym_link_cnvkit(self):
        jobs = []
    
        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".cnvkit.snpeff.annot.vcf"]
        
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_cnvkit.somatic." + tumor_pair.name + "." + key))

        return jobs
     
    def ensemble_metasv_somatic(self):
        """

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            ensemble_directory = os.path.join("SVariants", "ensemble", tumor_pair.name)

            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            isize_file = os.path.join("metrics", "dna", tumor_pair.tumor.name, "picard_metrics", "picard_metrics.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join("pairedVariants", "ensemble", tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vcf.gz")
            gatk_pass = os.path.join("pairedVariants", "ensemble", tumor_pair.name, tumor_pair.name + ".ensemble.somatic.flt.pass.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.somatic.vcf.gz")
            manta_vcf = os.path.abspath(os.path.join(pair_directory, tumor_pair.name + ".manta.somatic.vcf.gz"))
            wham_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.somatic.vcf.gz")
            delly_vcf= os.path.join(pair_directory, tumor_pair.name + ".delly.somatic.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(
                    isize_file
                )

            else:
                isize_mean = 325
                isize_sd = 75
                
            gatk_pass = None
            if os.path.isfile(gatk_vcf):
                jobs.append(concat_jobs([
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    bcftools.view(
                        gatk_vcf,
                        gatk_pass,
                        config.param('metasv_ensemble', 'filter_somatic_options')
                    ),
                ], name="metasv_ensemble.ensemble_pass." + tumor_pair.name))
            
            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                metasv.ensemble(
                    lumpy_vcf,
                    manta_vcf,
                    cnvkit_vcf,
                    wham_vcf,
                    delly_vcf,
                    gatk_pass,
                    inputTumor,
                    tumor_pair.tumor.name,
                    os.path.join(ensemble_directory, "rawMetaSV_somatic"),
                    ensemble_directory,
                    isize_mean=str(isize_mean),
                    isize_sd=str(isize_sd),
                    output_vcf=os.path.join(ensemble_directory, "variants.vcf.gz")
                ),
            ], name="metasv_ensemble." + tumor_pair.name))

        return jobs

    def ensemble_metasv_germline(self):
        """

        """
        jobs = []
    
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            ensemble_directory = os.path.join("SVariants", "ensemble", tumor_pair.name)
        
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            isize_file = os.path.join("metrics", "dna", tumor_pair.tumor.name, "picard_metrics",
                                      "picard_metrics.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join("pairedVariants", "ensemble", tumor_pair.name,
                                    tumor_pair.name + ".ensemble.germline.vcf.gz")
            gatk_pass = os.path.join("pairedVariants", "ensemble", tumor_pair.name,
                                     tumor_pair.name + ".ensemble.germline.flt.pass.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.germline.vcf.gz")
            manta_vcf = os.path.abspath(os.path.join(pair_directory, tumor_pair.name + ".manta.germline.vcf.gz"))
            wham_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.germline.vcf.gz")
            delly_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.germline.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(
                    isize_file
                )
        
            else:
                isize_mean = 325
                isize_sd = 75
        
            gatk_pass = None
            if os.path.isfile(gatk_vcf):
                jobs.append(concat_jobs([
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    bcftools.view(
                        gatk_vcf,
                        gatk_pass,
                        config.param('metasv_ensemble', 'filter_germline_options')
                    ),
                ], name="metasv_ensemble.ensemble_pass." + tumor_pair.name))
        
            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                metasv.ensemble(
                    lumpy_vcf,
                    manta_vcf,
                    cnvkit_vcf,
                    wham_vcf,
                    delly_vcf,
                    gatk_pass,
                    inputTumor,
                    tumor_pair.tumor.name,
                    os.path.join(ensemble_directory, "rawMetaSV_germline"),
                    ensemble_directory,
                    isize_mean=str(isize_mean),
                    isize_sd=str(isize_sd),
                    output_vcf=os.path.join(ensemble_directory, "variants.vcf.gz")
                ),
            ], name="metasv_ensemble." + tumor_pair.name))
    
        return jobs

    def metasv_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            ensemble_directory = os.path.join("SVariants", "ensemble", tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    os.path.join(ensemble_directory, "variants.vcf.gz"),
                    os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf")
                ),
                annotations.structural_variants(
                    os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf"),
                    os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.annot.vcf")
                ),
            ], name="sv_annotation.metasv_ensemble." + tumor_pair.name))

        return jobs

    def sym_link_metasv(self):
        jobs = []
    
        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", "ensemble", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".metasv.snpeff.annot.vcf"]
        
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_metasv." + tumor_pair.name + "." + key))
                    
        return jobs

    def scones(self):
        """
        This step aims to estimate somatic Copy Number Variation using BVAtools and SCoNEs. BVAtools generate the bined Depth ratio values from the
        tumor and normal BAM files. SCoNEs is tool to deconvolution the logR signal of the tumor-normal coverage into a mixture of baysian sub-signal
        for each copy number state. The result is a set of several deconvolution using  0-7 sub-signal. As each tumor sample is unique the choice of
        the best final model (number of sub-signal) needs to be manually evaluated using the log ratio graphical representation.

        """
        window_size = config.param('scones', 'window', required=True)
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            sv_directory = os.path.join("SVariants", tumor_pair.name)
            scones_directory = os.path.join(sv_directory, "SCoNEs")
            inputNormal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])[0]
            inputTumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])[0]

            bined_count_file = os.path.join(scones_directory, tumor_pair.normal.name + ".bin" + window_size + ".tsv")
            bined_count_fix_file = os.path.join(scones_directory, tumor_pair.normal.name + ".bin" + window_size + ".fix.tsv")

            output_scones_basename = os.path.join(scones_directory,
                                                  tumor_pair.normal.name + ".bin" + window_size + "_SCoNEs")
            scones_best_model_basename = output_scones_basename + "_Model_" + config.param('scones', 'best_model',
                                                                                           required=True)
            scones_calls_file = scones_best_model_basename + "_CNVcalls.txt"
            scones_filtered_file = scones_best_model_basename + "_CNVcalls.filtered.tsv"
            scones_annotate_basename = scones_best_model_basename + "_CNVcalls.filtered.anotated"
            scones_annotate_tmp_basename = scones_best_model_basename + "_CNVcalls.filtered.tmp"

            jobs.append(concat_jobs([
                bash.mkdir(
                    scones_directory,
                    remove=True
                ),
                bvatools.bincounter(
                    bam=inputTumor,
                    refbam=inputNormal,
                    out=bined_count_fix_file,
                    window=window_size
                ),
                Job(
                    [bined_count_fix_file],
                    [bined_count_file],
                    command="cat <(head -1 " + bined_count_fix_file + ") <(grep -v \"_\" " + bined_count_fix_file
                            + " | grep -v \"EBV\" ) > " + bined_count_file
                ),
            ], name="bvatools_bincounter." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    scones_directory,
                    remove=True
                ),
                scones.scones_pair(
                    bined_file=bined_count_file,
                    output_basename=output_scones_basename,
                    window=window_size
                )
            ], name="scones_pair." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    scones_directory,
                    remove=True
                ),
                scones.scones_filter(
                    scones_calls=scones_calls_file,
                    pair_name=tumor_pair.name,
                    output=scones_filtered_file
                )
            ], name="scones_filter." + tumor_pair.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    scones_directory,
                    remove=True
                ),
                scones.scones_annotate(
                    scones_calls_filtered=scones_filtered_file,
                    output_basename=scones_annotate_basename,
                    tmp_basename=scones_annotate_tmp_basename
                )
            ], name="scones_annotate." + tumor_pair.name))

        return jobs

    def svaba_assemble(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            svaba_directory = os.path.join(pair_directory, "rawSvaba")
            abs_alignment = os.path.abspath("alignment")
            input_normal = os.path.join(abs_alignment, tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(abs_alignment, tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            somatic_input = tumor_pair.name + ".svaba.somatic.sv.vcf"
            somatic_output = os.path.join(os.path.abspath(pair_directory), tumor_pair.name + ".svaba.somatic.vcf")

            germline_input = tumor_pair.name + ".svaba.germline.sv.vcf"
            germline_output = os.path.join(os.path.abspath(pair_directory), tumor_pair.name + ".svaba.germline.vcf")

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.tumor.readsets[0]
            )

            bed = None

            if coverage_bed:
                bed = coverage_bed

            jobs.append(concat_jobs([
                bash.mkdir(
                    svaba_directory,
                    remove=True
                ),
                Job(
                    command="cd " + svaba_directory
                ),
                svaba.run(
                    input_tumor,
                    tumor_pair.name,
                    input_normal,
                    bed
                ),
                Job(
                    [somatic_input],
                    [somatic_output],
                    command="sed -e 's#" + os.path.abspath(input_normal) + "#" + tumor_pair.normal.name + "#g' " + somatic_input + " | "
                                                               "sed -e 's#" + os.path.abspath(input_tumor) + "#" + tumor_pair.tumor.name + "#g' > " + somatic_output),
                Job([germline_input], [germline_output], command="sed -e 's#" + os.path.abspath(input_normal) + "#" + tumor_pair.normal.name + "#g' " + germline_input + " | "
                                                               "sed -e 's#" + os.path.abspath(input_tumor) + "#" + tumor_pair.tumor.name + "#g' > " + germline_output)
            ], name="svaba_run." + tumor_pair.name))

        return jobs

    def svaba_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(
                    [os.path.abspath(pair_directory) + ".svaba.somatic.vcf"],
                    [os.path.abspath(pair_directory) + ".svaba.somatic.flt.vcf"],
                    command="cat <(grep \"^#\" " + os.path.abspath(pair_directory)
                            + ".svaba.somatic.vcf) <(grep -v \"^#\" " + os.path.abspath(pair_directory)
                            + ".svaba.somatic.vcf | cut -f1-9,13-14) > " + os.path.abspath(pair_directory) + ".svaba.somatic.flt.vcf"
                ),
                snpeff.compute_effects(
                    os.path.abspath(pair_directory) + ".svaba.somatic.flt.vcf",
                    pair_directory + ".svaba.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    pair_directory + ".svaba.somatic.snpeff.vcf",
                    pair_directory + ".svaba.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    pair_directory + ".svaba.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "SVABA",
                    pair_directory + ".svaba.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.svaba_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(
                    [os.path.abspath(pair_directory) + ".svaba.germline.vcf"],
                    [os.path.abspath(pair_directory) + ".svaba.germline.flt.vcf"],
                    command="cat <(grep \"^#\" " + os.path.abspath(pair_directory) + ".svaba.germline.vcf) <(grep -v \"^#\" "
                            + os.path.abspath(pair_directory) + ".svaba.germline.vcf | cut -f1-9,13-14) > "
                            + os.path.abspath(pair_directory) + ".svaba.germline.flt.vcf"
                ),
                snpeff.compute_effects(os.path.abspath(pair_directory) + ".svaba.germline.flt.vcf",
                                       pair_directory + ".svaba.germline.snpeff.vcf"
                                       ),
                annotations.structural_variants(
                    pair_directory + ".svaba.germline.snpeff.vcf",
                    pair_directory + ".svaba.germline.snpeff.annot.vcf"
                ),
                vawk.sv(
                    pair_directory + ".svaba.germline.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "SVABA",
                    pair_directory + ".svaba.germline.prioritize.tsv"
                ),
            ], name="sv_annotation.svaba_germline." + tumor_pair.name))

        return jobs

    def sym_link_svaba(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".svaba.somatic.snpeff.annot.vcf",
                               pair_directory + ".svaba.somatic.prioritize.tsv"]
                               
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_svaba.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [pair_directory + ".svaba.germline.sv.snpeff.annot.vcf",
                               pair_directory + ".svaba.germline.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample,
                            sample + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample,
                            tumor_pair,
                            self.output_dir,
                            type="sv",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_svaba.germline." + tumor_pair.name + "." + key))

        return jobs

    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.rawmpileup_panel,
                self.paired_varscan2_panel,
                self.merge_varscan2_panel,
                self.preprocess_vcf_panel,
                self.snp_effect_panel,
                self.gemini_annotations_panel,
                #self.set_somatic_and_actionable_mutations_panel,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc,
                self.run_pair_multiqc,
                self.sym_link_report,
                self.sym_link_fastq_pair,
                self.sym_link_panel,
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc,
                self.rawmpileup,
                self.paired_varscan2,
                self.merge_varscan2,
                self.paired_mutect2,
                self.merge_mutect2,
                self.samtools_paired,
                self.merge_filter_paired_samtools,
                self.vardict_paired,
                self.merge_filter_paired_vardict,
                self.strelka2_paired_somatic,
                self.ensemble_somatic,
                self.gatk_variant_annotator_somatic,
                self.merge_gatk_variant_annotator_somatic,
                #self.somatic_signature,
                self.compute_cancer_effects_somatic,
                self.ensemble_somatic_dbnsfp_annotation,
                self.sample_gemini_annotations_somatic,
                #self.set_somatic_and_actionable_mutations,
                self.ensemble_germline_loh,
                self.gatk_variant_annotator_germline,
                self.merge_gatk_variant_annotator_germline,
                self.compute_cancer_effects_germline,
                self.sample_gemini_annotations_germline,
                #self.combine_tumor_pairs_somatic,
                #self.decompose_and_normalize_mnps_somatic,
                #self.all_pairs_compute_effects_somatic,
                #self.gemini_annotations_somatic,
                #self.combine_tumor_pairs_germline,
                #self.decompose_and_normalize_mnps_germline,
                #self.all_pairs_compute_effects_germline,
                #self.gemini_annotations_germline,
                self.run_pair_multiqc,
                self.sym_link_fastq_pair,
                self.sym_link_final_bam,
                self.sym_link_report,
                self.sym_link_ensemble,
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.metrics_dna_picard_metrics,
                self.sequenza,
                #self.sCNAphase,
                self.delly_call_filter,
                self.delly_sv_annotation,
                self.manta_sv_calls,
                self.manta_sv_annotation,
                self.lumpy_paired_sv,
                self.lumpy_sv_annotation,
                self.wham_call_sv,
                self.wham_sv_annotation,
                self.cnvkit_batch,
                self.cnvkit_sv_annotation,
                self.scones,
                self.svaba_assemble,
                self.svaba_sv_annotation,
                self.ensemble_metasv_somatic,
                self.ensemble_metasv_germline,
                self.metasv_sv_annotation,
                self.sym_link_sequenza,
                self.sym_link_metasv,
                self.sym_link_delly,
                self.sym_link_manta,
                self.sym_link_lumpy,
                self.sym_link_wham,
                self.sym_link_cnvkit,
                #self.sym_link_svaba
            ]
        ]


if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        TumorPair(protocol=['fastpass','ensemble','sv'])
=======
#!/usr/bin/env python################################################################################# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre## This file is part of MUGQIC Pipelines.## MUGQIC Pipelines is free software: you can redistribute it and/or modify# it under the terms of the GNU Lesser General Public License as published by# the Free Software Foundation, either version 3 of the License, or# (at your option) any later version.## MUGQIC Pipelines is distributed in the hope that it will be useful,# but WITHOUT ANY WARRANTY; without even the implied warranty of# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the# GNU Lesser General Public License for more details.## You should have received a copy of the GNU Lesser General Public License# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.################################################################################# Python Standard Modulesimport loggingimport mathimport osimport reimport sys# Append mugqic_pipelines directory to Python library pathsys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))# MUGQIC Modulesfrom core.config import *from core.job import *from core.pipeline import *from bfx.sample_tumor_pairs import *from bfx.sequence_dictionary import *from bfx import sambambafrom bfx import samtoolsfrom bfx import vcflibfrom bfx import htslibfrom bfx import varscanfrom bfx import bamreadcountfrom bfx import bcftoolsfrom bfx import gatkfrom bfx import toolsfrom bfx import bed_filefrom bfx import vardictfrom bfx import bcbio_variation_recallfrom bfx import vtfrom bfx import snpefffrom bfx import geminifrom bfx import conpairfrom bfx import multiqcfrom pipelines.dnaseq import dnaseqlog = logging.getLogger(__name__)class TumorPair(dnaseq.DnaSeq):    """    Tumor Pair Pipeline    =================    The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline with the exception of the    indel realignment (IR) step. In the tumor pipeline the IR step utilizes both the normal and tumor bam to further reduce    false positives (FPs) in and around indels. The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step.     At this point, a paired caller is used to call SNVs and Indels from the pairs given as input. Additional, muliple cancer callers     are utilized using an ensemble approach and SNVs and Indels seen in at least 2 different callers are retained for further     investigation.    Example command:    python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv        -c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin    -r readset: derived from GQ lims or made yourself. See : https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-readset-file    -p pairs : format - patient_name,normal_sample_name,tumor_sample_name     """    def __init__(self):        # Add pipeline specific arguments        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)        self.argparser.add_argument("--profyle", help="adjust deliverables to PROFYLE folder conventions (Default: False)", action="store_true")        super(TumorPair, self).__init__()    @property    def tumor_pairs(self):        if not hasattr(self, "_tumor_pairs"):            self._tumor_pairs = parse_tumor_pair_file(self.args.pairs.name, self.samples)        return self._tumor_pairs    def sequence_dictionary_variant(self):        if not hasattr(self, "_sequence_dictionary_variant"):            self._sequence_dictionary_variant = parse_sequence_dictionary_file(                config.param('DEFAULT', 'genome_dictionary', type='filepath'), variant=True)        return self._sequence_dictionary_variant    def generate_approximate_windows(self, nb_jobs):        if nb_jobs <= len(self.sequence_dictionary_variant()):            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in                    self.sequence_dictionary_variant()]        else:            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])            approximate_window_size = int(                math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))            windows = []            for sequence in self.sequence_dictionary_variant():                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in                                   range(1, sequence['length'] + 1, approximate_window_size)]:                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))        return windows    def gatk_indel_realigner(self):        """        Insertion and deletion realignment is performed on regions where multiple base mismatches        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.        Such regions will introduce false positive variant calls which may be filtered out by realigning        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).        The reference genome is divided by a number regions given by the `nb_jobs` parameter.        Note: modified to use both normal and tumor bams to reduce FPs around indels        """        jobs = []        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')        if nb_jobs > 50:            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")        for tumor_pair in self.tumor_pairs.itervalues():            normal_alignment_directory = os.path.join("alignment", tumor_pair.normal.name)            normal_realign_directory = os.path.join(normal_alignment_directory, "realign")            tumor_alignment_directory = os.path.join("alignment", tumor_pair.tumor.name)            tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")            if nb_jobs == 1:                realign_intervals = os.path.join(tumor_realign_directory, "all.intervals")                bam_postfix = ".all.realigned.bam"                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.all.realigned.bam")                normal_index = re.sub("\.bam$", ".bai", normal_bam)                normal_output_bam = os.path.join(normal_alignment_directory,                                                 tumor_pair.normal.name + ".realigned.qsorted.bam")                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.all.realigned.bam")                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)                tumor_output_bam = os.path.join(tumor_alignment_directory,                                                tumor_pair.tumor.name + ".realigned.qsorted.bam")                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)                jobs.append(concat_jobs([                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor),                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,                                         output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,                                         optional=bam_postfix),                    # Move sample realign                     Job([input_normal], [normal_output_bam],                        command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),                    Job([input_tumor], [tumor_output_bam],                        command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)                ], name="gatk_indel_realigner." + tumor_pair.name))            else:                # The first sequences are the longest to process.                # Each of them must be processed in a separate job.                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,                                                                                          nb_jobs - 1)                # Create one separate job for each of the first sequences                for idx, sequences in enumerate(unique_sequences_per_job):                    realign_prefix = os.path.join(tumor_realign_directory, str(idx))                    realign_intervals = realign_prefix + ".intervals"                    intervals = sequences                    if str(idx) == 0:                        intervals.append("unmapped")                    bam_postfix = ".realigned." + str(idx) + ".bam"                    normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")                    normal_index = re.sub("\.bam$", ".bai", normal_bam)                    tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)                    normal_output_bam = os.path.join(normal_realign_directory,                                                     tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)                    tumor_output_bam = os.path.join(tumor_realign_directory,                                                    tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)                    jobs.append(concat_jobs([                        # Create output directory since it is not done by default by GATK tools                        Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),                        Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),                        gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor,                                                      intervals=intervals),                        gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,                                             output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,                                             intervals=intervals, optional=bam_postfix),                        Job([input_normal], [normal_output_bam],                            command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),                        Job([input_tumor], [tumor_output_bam],                            command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))                # Create one last job to process the last remaining sequences and 'others' sequences                realign_prefix = os.path.join(tumor_realign_directory, "others")                realign_intervals = realign_prefix + ".intervals"                bam_postfix = ".realigned.others.bam"                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned.others.bam")                normal_index = re.sub("\.bam$", ".bai", normal_bam)                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned.others.bam")                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)                normal_output_bam = os.path.join(normal_realign_directory,                                                 tumor_pair.normal.name + ".sorted.realigned.others.bam")                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)                tumor_output_bam = os.path.join(tumor_realign_directory,                                                tumor_pair.tumor.name + ".sorted.realigned.others.bam")                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)                jobs.append(concat_jobs([                    # Create output directory since it is not done by default by GATK tools                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor,                                                  exclude_intervals=unique_sequences_per_job_others),                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,                                         output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,                                         exclude_intervals=unique_sequences_per_job_others, optional=bam_postfix),                    Job([input_normal], [normal_output_bam],                        command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),                    Job([input_tumor], [tumor_output_bam],                        command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))        return jobs    def sambamba_mark_duplicates(self):        """        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).        """        jobs = []        for sample in self.samples:            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")            input = alignment_file_prefix + "realigned.sorted.bam"            output = alignment_file_prefix + "sorted.dup.bam"            job = sambamba.markdup(input, output, os.path.join("alignment", sample.name, sample.name))            job.name = "sambamba_mark_duplicates." + sample.name            jobs.append(job)        report_file = os.path.join("report", "DnaSeq.picard_mark_duplicates.md")        jobs.append(            Job(                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in self.samples],                [report_file],                command="""\mkdir -p report && \\cp \\  {report_template_dir}/{basename_report_file} \\  {report_file}""".format(                    report_template_dir=self.report_template_dir,                    basename_report_file=os.path.basename(report_file),                    report_file=report_file                ),                report_files=[report_file],                name="picard_mark_duplicates_report")        )        return jobs    def conpair_concordance_contamination(self):        """                """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            metrics_directory = os.path.join("metrics")            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")            pileup_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".gatkPileup")            pileup_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".gatkPileup")            concordance_out = os.path.join(metrics_directory, tumor_pair.name + ".concordance.tsv")            contamination_out = os.path.join(metrics_directory, tumor_pair.name + ".contamination.tsv")            jobs.append(concat_jobs([                conpair.pileup(input_normal, pileup_normal),            ], name="conpair_concordance_contamination.pileup." + tumor_pair.normal.name))            jobs.append(concat_jobs([                conpair.pileup(input_tumor, pileup_tumor),            ], name="conpair_concordance_contamination.pileup." + tumor_pair.tumor.name))            jobs.append(concat_jobs([                Job(command="mkdir -p " + metrics_directory),                conpair.concordance(pileup_normal, pileup_tumor, concordance_out),                conpair.contamination(pileup_normal, pileup_tumor, contamination_out)            ], name="conpair_concordance_contamination." + tumor_pair.name))        return jobs    def rawmpileup_panel(self):        """        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.        One packaged mpileup file is created per sample/chromosome.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            bedfile = config.param('rawmpileup_panel', 'panel')            for sequence in self.sequence_dictionary_variant():                pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")                jobs.append(concat_jobs([                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),                    samtools.mpileup(                        [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam"),                         os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],                        pair_output, config.param('rawmpileup_panel', 'mpileup_other_options'), region=sequence['name'],                        regionFile=bedfile),                ], name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name']))        return jobs    def paired_varscan2_panel(self):        """        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            for sequence in self.sequence_dictionary_variant():                input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])                output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")                output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")                output_vcf_gz = os.path.join(varscan_directory,                                             tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")                jobs.append(concat_jobs([                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),                    varscan.somatic(input_pair, output, config.param('varscan2_somatic_panel', 'other_options'),                                    output_vcf_dep=output_vcf_gz, output_snp_dep=output_snp,                                    output_indel_dep=output_indel),                    htslib.bgzip_tabix_vcf(output_snp, os.path.join(varscan_directory,                                                                    tumor_pair.name + ".snp." + sequence[                                                                        'name'] + ".vcf.gz")),                    htslib.bgzip_tabix_vcf(output_indel, os.path.join(varscan_directory,                                                                      tumor_pair.name + ".indel." + sequence[                                                                          'name'] + ".vcf.gz")),                    pipe_jobs([                        bcftools.concat(                            [os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"),                             os.path.join(varscan_directory,                                          tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")], None),                        Job([None], [None],                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' "),                        htslib.bgzip_tabix_vcf(None, output_vcf_gz),                    ]),                ], name="varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name']))        return jobs    def merge_varscan2_panel(self):        """        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")                          for sequence in self.sequence_dictionary_variant()]            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")            somatic_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")            germline_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vcf.gz")            jobs.append(concat_jobs([                pipe_jobs([                    bcftools.concat(all_inputs, None),                    tools.fix_varscan_output(None, None),                    Job([None], [None],                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                    Job([None], [None],                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                    htslib.bgzip_tabix_vcf(None, all_output),                ]),                pipe_jobs([                    bcftools.view(all_output, None, config.param('merge_varscan2', 'somatic_filter_options')),                    htslib.bgzip_tabix_vcf(None, somatic_output),                ]),                pipe_jobs([                    bcftools.view(all_output, None, config.param('merge_varscan2', 'germline_loh_filter_options')),                    htslib.bgzip_tabix_vcf(None, germline_output),                ]),            ], name="merge_varscan2." + tumor_pair.name))        return jobs    def preprocess_vcf_panel(self):        """        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and        vcf FORMAT modification for correct loading into gemini        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            prefix = os.path.join(pair_directory, tumor_pair.name)            output_somatic = prefix + ".varscan2.somatic.vt.vcf.gz"            output_germline = prefix + ".varscan2.germline_loh.vt.vcf.gz"            jobs.append(concat_jobs([                pipe_jobs([                    vt.decompose_and_normalize_mnps(prefix + ".varscan2.somatic.vcf.gz", None),                    htslib.bgzip_tabix_vcf(None, prefix + ".prep.vt.vcf.gz"),                ]),                tools.preprocess_varscan(prefix + ".prep.vt.vcf.gz", output_somatic),            ], name="preprocess_vcf_panel.somatic." + tumor_pair.name))            jobs.append(concat_jobs([                pipe_jobs([                    vt.decompose_and_normalize_mnps(prefix + ".varscan2.germline_loh.vcf.gz", None),                    htslib.bgzip_tabix_vcf(None, prefix + ".germline_loh.prep.vt.vcf.gz"),                ]),                tools.preprocess_varscan(prefix + ".germline_loh.prep.vt.vcf.gz", output_germline),            ], name="preprocess_vcf_panel.germline." + tumor_pair.name))        return jobs    def snp_effect_panel(self):        """        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            if not os.path.exists(varscan_directory):                os.makedirs(varscan_directory)            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf")            output_germline_gz = os.path.join(pair_directory,                                              tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf.gz")            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')            cancer_pair = open(cancer_pair_filename, 'w')            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")            jobs.append(concat_jobs([                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename,                                       options=config.param('compute_cancer_effects_somatic', 'options')),                htslib.bgzip_tabix_vcf(output_somatic, output_somatic_gz),            ], name="compute_cancer_effects_somatic." + tumor_pair.name))            jobs.append(concat_jobs([                snpeff.compute_effects(input_germline, output_germline, cancer_sample_file=cancer_pair_filename,                                       options=config.param('compute_cancer_effects_germline', 'options')),                htslib.bgzip_tabix_vcf(output_germline, output_germline_gz),            ], name="compute_cancer_effects_germline." + tumor_pair.name))        return jobs    def gemini_annotations_panel(self):        """        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html        """        jobs = []        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            if not os.path.exists(varscan_directory):                os.makedirs(varscan_directory)            temp_dir = os.path.join(os.getcwd(), pair_directory)            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)            jobs.append(concat_jobs([                gemini.gemini_annotations(gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",                                          gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)            ], name="gemini_annotations.somatic." + tumor_pair.name))            jobs.append(concat_jobs([                gemini.gemini_annotations(gemini_prefix + ".varscan2.germline_loh.vt.snpeff.vcf.gz",                                          gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)            ], name="gemini_annotations.germline." + tumor_pair.name))        return jobs    def run_pair_multiqc(self):        jobs = []        metrics_directory = os.path.join("metrics", "dna")        inputs = []        outputs = []        for tumor_pair in self.tumor_pairs.itervalues():            input_normal_oxog = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.oxog_metrics.txt")            input_normal_qcbias = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.qcbias_metrics.pdf")            input_normal_all_picard = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.all.metrics.quality_distribution.pdf")            input_normal_qualimap = os.path.join(metrics_directory, tumor_pair.normal.name, "qualimap", tumor_pair.normal.name, "genome_results.txt")            input_normal_fastqc = os.path.join(metrics_directory, tumor_pair.normal.name, "fastqc", tumor_pair.normal.name + ".sorted.dup_fastqc.zip")            input_normal_flagstat = os.path.join(metrics_directory, tumor_pair.normal.name, "flagstat", tumor_pair.normal.name + ".flagstat")            input_tumor_oxog = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.oxog_metrics.txt")            input_tumor_qcbias = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.qcbias_metrics.pdf")            input_tumor_all_picard = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.all.metrics.quality_distribution.pdf")            input_tumor_qualimap = os.path.join(metrics_directory, tumor_pair.tumor.name, "qualimap", tumor_pair.tumor.name, "genome_results.txt")            input_tumor_fastqc = os.path.join(metrics_directory, tumor_pair.tumor.name, "fastqc", tumor_pair.tumor.name + ".sorted.dup_fastqc.zip")            input_tumor_flagstat = os.path.join(metrics_directory, tumor_pair.tumor.name, "flagstat", tumor_pair.tumor.name + ".flagstat")            input_dep = [input_normal_oxog, input_normal_qcbias, input_normal_all_picard, input_normal_qualimap, input_normal_fastqc,                         input_normal_flagstat, input_tumor_oxog, input_tumor_qcbias, input_tumor_all_picard, input_tumor_qualimap, input_tumor_fastqc,                        input_tumor_flagstat]            input = [os.path.join(metrics_directory, tumor_pair.normal.name), os.path.join(metrics_directory,tumor_pair.tumor.name)]            output = os.path.join(metrics_directory, tumor_pair.name + ".multiqc.html")            inputs.append(input)            outputs.append(output)            jobs.append(concat_jobs([                multiqc.run(input, output, input_dep=input_dep, sample=tumor_pair.name),            ], name="multiqc." + tumor_pair.name))        jobs.append(concat_jobs([            multiqc.run(inputs, os.path.join(metrics_directory, "allPairs.multiqc.html") , input_dep=outputs, sample="allPairs"),        ], name="multiqc.allPairs"))        return jobs    def rawmpileup(self):        """        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.        One packaged mpileup file is created per sample/chromosome.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            for sequence in self.sequence_dictionary_variant():                pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")                jobs.append(concat_jobs([                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),                    samtools.mpileup([os.path.join("alignment", tumor_pair.normal.name,                                                   tumor_pair.normal.name + ".sorted.dup.recal.bam"),                                      os.path.join("alignment", tumor_pair.tumor.name,                                                   tumor_pair.tumor.name + ".sorted.dup.recal.bam")], pair_output,                                     config.param('rawmpileup', 'mpileup_other_options'), sequence['name']),                ], name="rawmpileup." + tumor_pair.name + "." + sequence['name']))        return jobs    def rawmpileup_cat(self):        """        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            mpileup_normal_file_prefix = os.path.join(varscan_directory, tumor_pair.normal.name + ".")            mpileup_normal_inputs = [mpileup_normal_file_prefix + sequence['name'] + ".mpileup" for sequence in                                     self.sequence_dictionary_variant()]            mpileup_tumor_file_prefix = os.path.join(varscan_directory, tumor_pair.tumor.name + ".")            mpileup_tumor_inputs = [mpileup_tumor_file_prefix + sequence['name'] + ".mpileup" for sequence in                                    self.sequence_dictionary_variant()]            normal_output = mpileup_normal_file_prefix + "mpileup"            tumor_output = mpileup_tumor_file_prefix + "mpileup"            jobs.append(concat_jobs([                Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),                Job(mpileup_normal_inputs, [normal_output],                    command="cat \\\n  " + " \\\n  ".join(mpileup_normal_inputs) + " \\\n  > " + normal_output),                Job(mpileup_tumor_inputs, [tumor_output],                    command="cat \\\n  " + " \\\n  ".join(mpileup_tumor_inputs) + " \\\n  > " + tumor_output)            ], name="rawmpileup_cat." + tumor_pair.name))        return jobs    def paired_varscan2(self):        """        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.         Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing        Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases        SSC INFO field remove to prevent collison with Samtools output during ensemble                             """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            for sequence in self.sequence_dictionary_variant():                input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])                output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")                output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")                output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")                output_vcf_gz = os.path.join(varscan_directory,                                             tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")                jobs.append(concat_jobs([                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),                    varscan.somatic(input_pair, output, config.param('varscan2_somatic', 'other_options'),                                    output_vcf_dep=output_vcf, output_snp_dep=output_snp,                                    output_indel_dep=output_indel),                    htslib.bgzip_tabix_vcf(output_snp, os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                        'name'] + ".snp.vcf.gz")),                    htslib.bgzip_tabix_vcf(output_indel, os.path.join(varscan_directory,                                                                      tumor_pair.name + "." + sequence[                                                                          'name'] + ".indel.vcf.gz")),                    pipe_jobs([                        bcftools.concat(                            [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"),                             os.path.join(varscan_directory,                                          tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")], None),                        Job([None], [output_vcf],                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > " + output_vcf),                    ]),                    htslib.bgzip_tabix_vcf(output_vcf, output_vcf_gz),                ], name="varscan2_somatic." + tumor_pair.name + "." + sequence['name']))        return jobs    def varscan2_fpfilter(self):        """        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing        Varscan2 filtering thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases        Somatic and germline calls are filter using bam-readcount info from tumor and normal bam, respectively.         Variants denoted as PASS and somatic (SS=1) or germline (SS=2) and loh (SS=3) are retained for further analysis.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")            for sequence in self.sequence_dictionary_variant():                input_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")                output_bed = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.bed")                output_normal_readcount = os.path.join(varscan_directory,                                                       tumor_pair.name + "." + sequence['name'] + ".normal.readcount")                output_tumor_readcount = os.path.join(varscan_directory,                                                      tumor_pair.name + "." + sequence['name'] + ".tumor.readcount")                output_fpfilter_somatic = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.fpfilter.somatic.vcf")                output_fpfilter_somatic_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.fpfilter.somatic.vcf.gz")                output_fpfilter_germline = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.fpfilter.germline_loh.vcf")                output_fpfilter_germline_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.fpfilter.germline_loh.vcf.gz")                output_vcf_somatic_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.somatic.vcf.gz")                output_vcf_germline_loh_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence[                    'name'] + ".varscan2.germline_loh.vcf.gz")                jobs.append(concat_jobs([                    Job(command="mkdir -p " + varscan_directory,                        removable_files=[varscan_directory, output_fpfilter_somatic, output_fpfilter_germline]),                    tools.vcf2bed(input_vcf, output_bed),                    bamreadcount.readcount(input_normal, output_bed, output_normal_readcount),                    bamreadcount.readcount(input_tumor, output_bed, output_tumor_readcount),                    varscan.fpfilter_somatic(input_vcf, output_tumor_readcount, output_fpfilter_somatic),                    htslib.bgzip_tabix_vcf(output_fpfilter_somatic, output_fpfilter_somatic_gz),                    pipe_jobs([                        bcftools.view(output_fpfilter_somatic_gz, None,                                      config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),                        vcflib.vcfstreamsort(None, None),                        htslib.bgzip_tabix_vcf(None, output_vcf_somatic_gz),                    ]),                    varscan.fpfilter_somatic(input_vcf, output_normal_readcount, output_fpfilter_germline),                    htslib.bgzip_tabix_vcf(output_fpfilter_germline, output_fpfilter_germline_gz),                    pipe_jobs([                        bcftools.view(output_fpfilter_germline_gz, None,                                      config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),                        vcflib.vcfstreamsort(None, None),                        htslib.bgzip_tabix_vcf(None, output_vcf_germline_loh_gz),                    ]),                ], name="varscan2_readcount_fpfilter." + tumor_pair.name + "." + sequence['name']))        return jobs    def merge_varscan2(self):        """        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.        """        jobs = []        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            varscan_directory = os.path.join(pair_directory, "rawVarscan2")            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")                          for sequence in self.sequence_dictionary_variant()]            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")            jobs.append(concat_jobs([                pipe_jobs([                    bcftools.concat(all_inputs, None),                    tools.fix_varscan_output(None, None),                    Job([None], [None],                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                    Job([None], [None],                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                    htslib.bgzip_tabix_vcf(None, all_output),                ]),                pipe_jobs([                    vt.decompose_and_normalize_mnps(all_output, None),                    htslib.bgzip_tabix_vcf(None, all_output_vt),                ]),                pipe_jobs([                    bcftools.view(all_output_vt, None,                                  config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),                    htslib.bgzip_tabix_vcf(None, somtic_output_vt),                ]),                pipe_jobs([                    bcftools.view(all_output_vt, None,                                  config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),                    htslib.bgzip_tabix_vcf(None, germline_output_vt),                ]),            ], name="merge_varscan2." + tumor_pair.name))        return jobs    def paired_mutect2(self):        """        GATK MuTect2 caller for SNVs and Indels.        """        jobs = []        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')        if nb_jobs > 50:            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            mutect_directory = os.path.join(pair_directory, "rawMuTect2")            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")            mkdir_job = Job(command="mkdir -p " + mutect_directory, removable_files=[mutect_directory])            if nb_jobs == 1:                jobs.append(concat_jobs([                    # Create output directory since it is not done by default by GATK tools                    mkdir_job,                    gatk.mutect2(input_normal, input_tumor,                                 os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"))                ], name="gatk_mutect2." + tumor_pair.name))            else:                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(                    self.sequence_dictionary_variant(), nb_jobs - 1)                # Create one separate job for each of the first sequences                for idx, sequences in enumerate(unique_sequences_per_job):                    outprefix = tumor_pair.name + "." + str(idx) + ".mutect2"                    jobs.append(concat_jobs([                        # Create output directory since it is not done by default by GATK tools                        mkdir_job,                        gatk.mutect2(input_normal, input_tumor, os.path.join(mutect_directory, outprefix + ".vcf.gz"),                                     intervals=sequences)                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))                # Create one last job to process the last remaining sequences and 'others' sequences                jobs.append(concat_jobs([                    # Create output directory since it is not done by default by GATK tools                    mkdir_job,                    gatk.mutect2(input_normal, input_tumor,                                 os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"),                                 exclude_intervals=unique_sequences_per_job_others)                ], name="gatk_mutect2." + tumor_pair.name + ".others"))        return jobs    def merge_mutect2(self):        """        Merge SNVs and indels for mutect2        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names        Generate a somatic vcf containing only PASS variants                """        jobs = []        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            mutect_directory = os.path.join(pair_directory, "rawMuTect2")            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")            if nb_jobs == 1:                input_vcf = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")                jobs.append(concat_jobs([                    Job([input_vcf], [output_gz], command="ln -s -f " + input_vcf + " " + output_gz),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output_gz, None),                        Job([output_vt_gz], [None],                            command="zcat " + output_vt_gz + " | sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),                        htslib.bgzip_tabix_vcf(None, output_somatic_vt),                    ]),                ], name="symlink_mutect_vcf." + tumor_pair.name))            elif nb_jobs > 1:                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(                    self.sequence_dictionary_variant(), nb_jobs - 1)                # Create one separate job for each of the first sequences                inputs = []                for idx, sequences in enumerate(unique_sequences_per_job):                    inputs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))                inputs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))                jobs.append(concat_jobs([                    pipe_jobs([                        bcftools.concat(inputs, None, config.param('merge_filter_mutect2', 'bcftools_options')),                        Job([None], [None],                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),                        htslib.bgzip_tabix_vcf(None, output_gz),                    ]),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output_gz, None),                        htslib.bgzip_tabix_vcf(None, output_vt_gz),                    ]),                    pipe_jobs([                        bcftools.view(output_vt_gz, None, config.param('merge_filter_mutect2', 'filter_options')),                        htslib.bgzip_tabix_vcf(None, output_somatic_vt),                    ]),                ], name="merge_filter_mutect2." + tumor_pair.name))        return jobs    def samtools_paired(self):        """        Samtools caller for SNVs and Indels using verison 0.1.19.        """        jobs = []        nb_jobs = config.param('samtools_paired', 'nb_jobs', type='posint')        if nb_jobs > 50:            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            samtools_directory = os.path.join(pair_directory, "rawSamtools")            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")            paired_sample = [input_normal, input_tumor]            if nb_jobs == 1:                jobs.append(concat_jobs([                    Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory]),                    pipe_jobs([                        samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'),                                         ini_section="samtools_paired"),                        samtools.bcftools_call_pair("-", os.path.join(samtools_directory, tumor_pair.name + ".bcf"),                                                    config.param('samtools_paired', 'bcftools_view_options'),                                                    pair_calling=True),                    ]),                ], name="samtools_paired." + tumor_pair.name))            else:                for region in self.generate_approximate_windows(                        nb_jobs):  # for idx,sequences in enumerate(unique_sequences_per_job):                    jobs.append(concat_jobs([                        Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory]),                        pipe_jobs([                            samtools.mpileup(paired_sample, None,                                             config.param('samtools_paired', 'mpileup_other_options'), region,                                             ini_section="samtools_paired"),                            samtools.bcftools_call_pair("-", os.path.join(samtools_directory,                                                                          tumor_pair.name + "." + region + ".bcf"),                                                        config.param('samtools_paired', 'bcftools_view_options'),                                                        pair_calling=True),                        ]),                    ], name="samtools_paired." + tumor_pair.name + "." + region))        return jobs    def merge_filter_paired_samtools(self):        """        bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.        The output of bcftools is fed to varfilter, which does an additional filtering of the variants        and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls        for all samples in the experiment.        Additional somatic filters are performed to reduce the number of FPs:         1. vcflibs vcfsamplediff tags each variant with <tag>={germline,somatic,loh} to specify the type         of variant given the genotype difference between the two samples.        2. bcftools filter is used to retain only variants with CLR>=15 and have STATUS=somatic from         vcfsamplediff        3. bcftools filter is used to retain only variants that have STATUS=germline or STATUS=loh from        vcfsamplediff        """        jobs = []        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            samtools_directory = os.path.join(pair_directory, "rawSamtools")            output = os.path.join(samtools_directory, tumor_pair.name + ".samtools.bcf")            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".samtools.vcf.gz")            output_vcf_vt = os.path.join(pair_directory, tumor_pair.name + ".samtools.vt.vcf.gz")            output_somatics = os.path.join(pair_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")            if nb_jobs == 1:                inputs = os.path.join(samtools_directory, tumor_pair.name + ".bcf")                jobs.append(concat_jobs([                    samtools.bcftools_cat_pair(inputs, output),                    pipe_jobs([                        samtools.bcftools_view_pair(output, None),                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                        htslib.bgzip_tabix_vcf(None, output_vcf),                    ]),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output_vcf, None),                        htslib.bgzip_tabix_vcf(None, output_vcf_vt),                    ]),                    pipe_jobs([                        bcftools.filter(output_vcf_vt, None,                                        config.param('merge_filter_paired_samtools', 'somatic_filter_options')),                        vcflib.vcffilter(None, None,                                         config.param('merge_filter_paired_samtools', 'somatic_vcffilter_options')),                        htslib.bgzip_tabix_vcf(None, output_somatics),                    ]),                    pipe_jobs([                        bcftools.filter(output_vcf_vt, None,                                        config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_germline_loh),                    ]),                ], name="merge_filter_paired_samtools." + tumor_pair.name))            else:                inputs = [os.path.join(samtools_directory, tumor_pair.name + "." + region + ".bcf") for region in                          self.generate_approximate_windows(                              nb_jobs)]  # for idx,sequences in enumerate(unique_sequences_per_job):                jobs.append(concat_jobs([                    samtools.bcftools_cat_pair(inputs, output),                    pipe_jobs([                        samtools.bcftools_view_pair(output, None),                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                        htslib.bgzip_tabix_vcf(None, output_vcf),                    ]),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output_vcf, None),                        htslib.bgzip_tabix_vcf(None, output_vcf_vt),                    ]),                    pipe_jobs([                        bcftools.filter(output_vcf_vt, None,                                        config.param('merge_filter_paired_samtools', 'somatic_filter_options')),                        # vcflib.vcffilter(None, None, config.param('merge_filter_paired_samtools', 'somatic_vcffilter_options')),                        htslib.bgzip_tabix_vcf(None, output_somatics),                    ]),                    pipe_jobs([                        bcftools.filter(output_vcf_vt, None,                                        config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_germline_loh),                    ]),                ], name="merge_filter_paired_samtools." + tumor_pair.name))        report_file = os.path.join("report", "DnaSeq.merge_filter_bcf.md")        jobs.append(            Job(                [output_vcf],                [report_file],                command="""\mkdir -p report && \\cp \\  {report_template_dir}/{basename_report_file} \\  {report_template_dir}/HumanVCFformatDescriptor.tsv \\  report/ && \\sed 's/\t/|/g' report/HumanVCFformatDescriptor.tsv | sed '2i-----|-----' >> {report_file}""".format(                    report_template_dir=self.report_template_dir,                    basename_report_file=os.path.basename(report_file),                    report_file=report_file                ),                report_files=[report_file],                name="merge_filter_bcf_report")        )        return jobs    def vardict_paired(self):        """        vardict caller for SNVs and Indels.        Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is AUPAC nomenclature         """        ##TO DO - the BED system needs to be revisted !!         jobs = []        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')        if nb_jobs > 50:            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')        bed_file_list = []        if use_bed:            bed = self.samples[0].readsets[0].beds[0]            bed_intervals, interval_size = bed_file.parse_bed_file(bed)            last_bed_file = 'vardict.tmp.' + str(nb_jobs - 1) + '.bed'            if not os.path.exists(last_bed_file):                bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="./vardict.tmp")            else:                for idx in range(nb_jobs):                    bed_file_list.append(os.path.join("vardict.tmp." + str(idx) + ".bed"))        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            vardict_directory = os.path.join(pair_directory, "rawVardict")            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")            mkdir_job = Job(command="mkdir -p " + vardict_directory, removable_files=[vardict_directory])            if use_bed:                idx = 0                for bf in bed_file_list:                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")                    jobs.append(concat_jobs([                        mkdir_job,                        pipe_jobs([                            vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, bf),                            vardict.testsomatic(None, None),                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),                            htslib.bgzip_tabix_vcf(None, output),                        ]),                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))                    idx += 1            else:                beds = []                for idx in range(nb_jobs):                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))                if nb_jobs == 1:                    bedjob = vardict.dict2beds(genome_dictionary, beds)                    output = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")                    jobs.append(concat_jobs([                        mkdir_job,                        bedjob,                        pipe_jobs([                            vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, beds.pop()),                            vardict.testsomatic(None, None),                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),                            htslib.bgzip_tabix_vcf(None, output),                        ]),                    ], name="vardict_paired." + tumor_pair.name + ".0"))                else:                    bedjob = vardict.dict2beds(genome_dictionary, beds)                    jobs.append(concat_jobs([mkdir_job, bedjob], name="vardict.genome.beds." + tumor_pair.name))                    for idx in range(nb_jobs):                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")                        jobs.append(concat_jobs([                            mkdir_job,                            pipe_jobs([                                vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, beds[idx]),                                vardict.testsomatic(None, None),                                vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),                                htslib.bgzip_tabix_vcf(None, output),                            ]),                        ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))        return jobs    def merge_filter_paired_vardict(self):        """        The fully merged vcf is filtered using following steps:        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic        2. Somatics identified in step 1 must have PASS filter        """        jobs = []        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')        for tumor_pair in self.tumor_pairs.itervalues():            pair_directory = os.path.join("pairedVariants", tumor_pair.name)            vardict_directory = os.path.join(pair_directory, "rawVardict")            output_tmp = os.path.join(pair_directory, tumor_pair.name + ".vardict.tmp.vcf.gz")            output = os.path.join(pair_directory, tumor_pair.name + ".vardict.vcf.gz")            output_vt = os.path.join(pair_directory, tumor_pair.name + ".vardict.vt.vcf.gz")            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")            if nb_jobs == 1:                inputs = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")                jobs.append(concat_jobs([                    Job([inputs], [output_tmp], command="ln -s -f " + inputs + " " + output_tmp),                    pipe_jobs([                        Job([output_tmp], [None],                            command="zcat {output} | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                        htslib.bgzip_tabix_vcf(None, output)                    ]),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output, None),                        htslib.bgzip_tabix_vcf(None, output_vt),                    ]),                    pipe_jobs([                        bcftools.view(output_vt, None,                                      config.param('merge_filter_paired_vardict', 'somatic_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_somatic),                    ]),                    pipe_jobs([                        bcftools.view(output_vt, None,                                      config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_germline_loh),                    ]),                ], name="symlink_vardict_vcf." + tumor_pair.name))            else:                inputVCFs = []                for idx in range(nb_jobs):                    inputVCFs.append(                        os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))                jobs.append(concat_jobs([                    pipe_jobs([                        bcftools.concat(inputVCFs, None),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),                        Job([None], [None],                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),                        htslib.bgzip_tabix_vcf(None, output),                    ]),                    pipe_jobs([                        vt.decompose_and_normalize_mnps(output, None),                        htslib.bgzip_tabix_vcf(None, output_vt),                    ]),                    pipe_jobs([                        bcftools.view(output_vt, None,                                      config.param('merge_filter_paired_vardict', 'somatic_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_somatic),                    ]),                    pipe_jobs([                        bcftools.view(output_vt, None,                                      config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),                        htslib.bgzip_tabix_vcf(None, output_germline_loh),                    ]),                ], name="merge_filter_paired_vardict." + tumor_pair.name))        return jobs    def ensemble_somatic(self):        """        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls        Filter ensemble calls to retain only calls overlapping 2 or more callers        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)            input_directory = os.path.join("pairedVariants", tumor_pair.name)            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")            inputs_somatic = [input_mutect2, input_vardict, input_varscan2, input_samtools]            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")            #output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[                os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic-work"), output_ensemble])            jobs.append(concat_jobs([                # Create output directory since it is not done by default by GATK tools                mkdir_job,                bcbio_variation_recall.ensemble(inputs_somatic, output_ensemble,                                                config.param('bcbio_ensemble_somatic', 'options')),            ], name="bcbio_ensemble_somatic." + tumor_pair.name))        return jobs    def ensemble_germline_loh(self):        """        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls        Filter ensemble calls to retain only calls overlapping 2 or more callers        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)            input_directory = os.path.join("pairedVariants", tumor_pair.name)            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")            inputs_germline = [input_vardict, input_varscan2, input_samtools]            output_ensemble = os.path.join(paired_ensemble_directory,                                           tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[                os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh-work"),                output_ensemble])            jobs.append(concat_jobs([                # Create output directory since it is not done by default by GATK tools                mkdir_job,                bcbio_variation_recall.ensemble(inputs_germline, output_ensemble,                                                config.param('bcbio_ensemble_germline_loh', 'options')),            ], name="bcbio_ensemble_germline_loh." + tumor_pair.name))        return jobs    def gatk_variant_annotator_somatic(self):        """        Add vcf annotations to ensemble vcf: Standard and Somatic annotations        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name,                                                  tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")            for sequence in self.sequence_dictionary_variant():                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name,                                                       tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence[                                                           'name'] + ".vcf.gz")                mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_somatic_variants])                jobs.append(concat_jobs([                    mkdir_job,                    gatk.variant_annotator(input_normal, input_tumor, input_somatic_variants, output_somatic_variants,                                           intervals=[sequence['name']]),                ], name="gatk_variant_annotator.somatic." + sequence['name'] + "." + tumor_pair.name))        return jobs    def gatk_variant_annotator_germline(self):        """        Add vcf annotations to ensemble vcf: most importantly the AD field        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            input_normal = os.path.join("alignment", tumor_pair.normal.name,                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,                                        tumor_pair.tumor.name + ".sorted.dup.recal.bam")            input_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name,                                                       tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")            for sequence in self.sequence_dictionary_variant():                output_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name,                                                            tumor_pair.name + ".ensemble.germline_loh.vt.annot." +                                                            sequence['name'] + ".vcf.gz")                mkdir_job = Job(command="mkdir -p " + ensemble_directory,                                removable_files=[output_germline_loh_variants])                jobs.append(concat_jobs([                    mkdir_job,                    gatk.variant_annotator(input_normal, input_tumor, input_germline_loh_variants,                                           output_germline_loh_variants, intervals=[sequence['name']]),                ], name="gatk_variant_annotator.germline." + sequence['name'] + "." + tumor_pair.name))        return jobs    def merge_gatk_variant_annotator_somatic(self):        """        Merge annotated somatic vcfs        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            inputs_somatic = [os.path.join(ensemble_directory, tumor_pair.name,                                           tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence[                                               'name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]            ouputs_somatic = os.path.join(ensemble_directory, tumor_pair.name,                                          tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")            jobs.append(concat_jobs([                pipe_jobs([                    bcftools.concat(inputs_somatic, None),                    htslib.bgzip_tabix_vcf(None, ouputs_somatic),                ]),            ], name="merge_gatk_variant_annotator.somatic." + tumor_pair.name))        return jobs    def merge_gatk_variant_annotator_germline(self):        """        Merge annotated germline and LOH vcfs        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            germline_inputs = [os.path.join(ensemble_directory, tumor_pair.name,                                            tumor_pair.name + ".ensemble.germline_loh.vt.annot." + sequence[                                                'name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]            germline_output = os.path.join(ensemble_directory, tumor_pair.name,                                           tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")            jobs.append(concat_jobs([                pipe_jobs([                    bcftools.concat(germline_inputs, None),                    htslib.bgzip_tabix_vcf(None, germline_output),                ]),            ], name="merge_gatk_variant_annotator.germline." + tumor_pair.name))        return jobs    def compute_cancer_effects_somatic(self):        """        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).        Modified arguments to consider paired cancer data.        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        if not os.path.exists(ensemble_directory):            os.makedirs(ensemble_directory)        for tumor_pair in self.tumor_pairs.itervalues():            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)            if not os.path.exists(paired_directory):                os.makedirs(paired_directory)            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")            output_somatic_gz = os.path.join(paired_directory,                                             tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')            cancer_pair = open(cancer_pair_filename, 'w')            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_somatic])            jobs.append(concat_jobs([                mkdir_job,                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename,                                       options=config.param('compute_cancer_effects_somatic', 'options')),                htslib.bgzip_tabix_vcf(output_somatic, output_somatic_gz),            ], name="compute_cancer_effects_somatic." + tumor_pair.name))        return jobs    def compute_cancer_effects_germline(self):        """        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).        Modified arguments to consider paired cancer data.        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        for tumor_pair in self.tumor_pairs.itervalues():            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")            output_germline = os.path.join(paired_directory,                                           tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf")            output_germline_gz = os.path.join(paired_directory,                                              tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz")            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')            cancer_pair = open(cancer_pair_filename, 'w')            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_germline])            jobs.append(concat_jobs([                mkdir_job,                snpeff.compute_effects(input_germline, output_germline,                                       options=config.param('compute_cancer_effects_germline', 'options')),                htslib.bgzip_tabix_vcf(output_germline, output_germline_gz),            ], name="compute_cancer_effects_germline." + tumor_pair.name))        return jobs    def combine_tumor_pairs_somatic(self):        """        Combine numerous ensemble vcfs into one vcf for gemini annotations        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input_merged_vcfs = [            os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for            tumor_pair in self.tumor_pairs.itervalues()]        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])        if len(input_merged_vcfs) == 1:            jobs.append(concat_jobs([                mkdir_job,                Job([input_merged_vcfs[0]], [output],                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)            ], name="gatk_combine_variants.somatic.allPairs"))        else:            jobs.append(concat_jobs([                mkdir_job,                gatk.combine_variants(input_merged_vcfs, output)            ], name="gatk_combine_variants.somatic.allPairs"))        return jobs    def combine_tumor_pairs_germline(self):        """        Combine numerous ensemble vcfs into one vcf for gemini annotations        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name,                                          tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz") for tumor_pair in                             self.tumor_pairs.itervalues()]        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])        if len(input_merged_vcfs) == 1:            jobs.append(concat_jobs([                mkdir_job,                Job([input_merged_vcfs[0]], [output],                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)            ], name="gatk_combine_variants.germline_loh.allPairs"))        else:            jobs.append(concat_jobs([                mkdir_job,                gatk.combine_variants(input_merged_vcfs, output)            ], name="gatk_combine_variants.germline_loh.allPairs"))        return jobs    def decompose_and_normalize_mnps_somatic(self):        """        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])        jobs.append(concat_jobs([            mkdir_job,            vt.decompose_and_normalize_mnps(input, output)        ], name="decompose_and_normalize_mnps.somatic.allPairs"))        return jobs    def decompose_and_normalize_mnps_germline(self):        """        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vcf.gz")        output_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")        job = vt.decompose_and_normalize_mnps(input_vcf, output_vcf)        job.name = "decompose_and_normalize_mnps.germline.allPairs"        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_vcf])        jobs.append(concat_jobs([            mkdir_job,            vt.decompose_and_normalize_mnps(input_vcf, output_vcf)        ], name="decompose_and_normalize_mnps.somatic.allPairs"))        return jobs    def all_pairs_compute_effects_somatic(self):        """        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).        Modified arguments to consider paired cancer data.        Applied to all tumor pairs.        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')        cancer_pair = open(cancer_pair_filename, 'w')        for tumor_pair in self.tumor_pairs.itervalues():            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])        jobs.append(concat_jobs([            mkdir_job,            snpeff.compute_effects(input, output, cancer_sample_file=cancer_pair_filename,                                   options=config.param('compute_cancer_effects_somatic', 'options')),            htslib.bgzip_tabix_vcf(output, output_gz),        ], name="compute_effects.somatic.allPairs"))        return jobs    def all_pairs_compute_effects_germline(self):        """        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).        Modified arguments to consider paired cancer data.        Applied to all tumor pairs.        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf")        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf.gz")        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])        jobs.append(concat_jobs([            mkdir_job,            snpeff.compute_effects(input, output, options=config.param('compute_cancer_effects_germline', 'options')),            htslib.bgzip_tabix_vcf(output, output_gz),        ], name="compute_effects.germline.allPair"))        return jobs    def gemini_annotations_somatic(self):        """        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        temp_dir = os.path.join(os.getcwd(), ensemble_directory)        gemini_prefix = os.path.join(ensemble_directory, "allPairs")        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])        jobs.append(concat_jobs([            Job(command="mkdir -p " + ensemble_directory),            gemini.gemini_annotations(gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",                                      gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)        ], name="gemini_annotations.somatic.allPairs"))        return jobs    def gemini_annotations_germline(self):        """        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html        """        jobs = []        ensemble_directory = os.path.join("pairedVariants", "ensemble")        temp_dir = os.path.join(os.getcwd(), ensemble_directory)        gemini_prefix = os.path.join(ensemble_directory, "allPairs")        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])        jobs.append(concat_jobs([            Job(command="mkdir -p " + ensemble_directory),            gemini.gemini_annotations(gemini_prefix + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz",                                      gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)        ], name="gemini_annotations.germline.allPairs"))        return jobs    @property    def steps(self):        return [            self.picard_sam_to_fastq,            self.trimmomatic,            self.merge_trimmomatic_stats,            self.skewer_trimming,            self.bwa_mem_picard_sort_sam,            self.sambamba_merge_sam_files,            self.gatk_indel_realigner,            self.sambamba_merge_realigned,            # self.fix_mate_by_coordinate,            self.sambamba_mark_duplicates,            self.recalibration,            self.conpair_concordance_contamination,            self.rawmpileup_panel,            self.paired_varscan2_panel,            self.merge_varscan2_panel,            self.preprocess_vcf_panel,            self.snp_effect_panel,            self.gemini_annotations_panel,            self.metrics_dna_damage_estimation,            self.metrics_dna_picard_metrics,            self.metrics_dna_sample_qualimap,            self.metrics_dna_fastqc,            self.run_pair_multiqc,            self.picard_calculate_hs_metrics,            self.gatk_callable_loci,            self.extract_common_snp_freq,            self.baf_plot,            self.rawmpileup,            self.paired_varscan2,            # self.varscan2_fpfilter,            self.merge_varscan2,            self.paired_mutect2,            self.merge_mutect2,            self.samtools_paired,            self.merge_filter_paired_samtools,            self.vardict_paired,            self.merge_filter_paired_vardict,            self.ensemble_somatic,            self.gatk_variant_annotator_somatic,            self.merge_gatk_variant_annotator_somatic,            self.compute_cancer_effects_somatic,            self.combine_tumor_pairs_somatic,            self.decompose_and_normalize_mnps_somatic,            self.all_pairs_compute_effects_somatic,            self.gemini_annotations_somatic,            self.ensemble_germline_loh,            self.gatk_variant_annotator_germline,            self.merge_gatk_variant_annotator_germline,            self.compute_cancer_effects_germline,            self.combine_tumor_pairs_germline,            self.decompose_and_normalize_mnps_germline,            self.all_pairs_compute_effects_germline,            self.gemini_annotations_germline        ]if __name__ == '__main__':    TumorPair()
>>>>>>> 5449fa52 (updates to config)
