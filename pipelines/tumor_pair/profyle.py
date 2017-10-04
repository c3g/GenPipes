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
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.sample_tumor_pairs import *
from bfx.sequence_dictionary import *

from bfx import picard
from bfx import samtools
from bfx import vcflib
from bfx import htslib
from bfx import varscan
from bfx import bamreadcount
from bfx import bcftools
from bfx import gq_seq_utils
from bfx import gatk
from bfx import tools
from bfx import bed_file
from bfx import vardict
from bfx import bcbio_variation
from bfx import vt
from bfx import snpeff
from bfx import gemini
from pipelines.dnaseq import dnaseq
from pipelines.tumor_pair import tumor_pair

log = logging.getLogger(__name__)

class profyle_dna(tumor_pair.TumorPair):
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
        self._protocol=protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)
        super(TumorPair, self).__init__(protocol)

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(self.args.pairs.name, self.samples)
        return self._tumor_pairs
    @property
    def fastpass_bed(self):
        if not hasattr(self, "_fastpass_bed"):
            self._fastpass_bed = config.param('DEFAULT', 'fastpas_bed', type='filepath')
        return self._fastpass_bed

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', type='filepath'), variant=True)
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
            normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
            tumor_alignment_directory = os.path.join("alignment", tumor_pair.tumor.name)
            tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")
            
            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.join(tumor_realign_directory, "all.intervals")
                bam_postfix = ".all.realigned.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.all.realigned.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".realigned.qsorted.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_bam = os.path.join( tumor_pair.tumor.name + ".sorted.all.realigned.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".realigned.qsorted.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, optional=bam_postfix),
                    # Move sample realign 
                    Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(tumor_realign_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned." + str(idx)  + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned." + str(idx)  + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),
                        Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),
                        gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor, intervals=intervals),
                        gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, intervals=intervals, optional=bam_postfix),
                        Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                        Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(tumor_realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor, exclude_intervals=unique_sequences_per_job_others),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam, output_tum_dep=tumor_output_bam, target_intervals=realign_intervals, exclude_intervals=unique_sequences_per_job_others, optional=bam_postfix),
                    Job([input_normal], [normal_output_bam], command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam], command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))

        return jobs

    def merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")
              
            if nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                inputBAMs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned." + str(idx) + ".bam"))
                inputBAMs.append(os.path.join(realign_directory, sample.name + ".sorted.realigned.others.bam"))

                job = picard.merge_sam_files(inputBAMs, merged_realigned_bam)
                job.name = "merge_realigned." + sample.name
                jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.gatk_indel_realigner.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".realigned.qsorted.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_realigned_report")
        )

        return jobs

    def rawmpileup(self,fast_bed=None,fast_name=None):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                normal_output = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + "." + fast_name + ".mpileup" if fast_name else tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                tumor_output = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + "." + fast_name + ".mpileup" if fast_name else tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")], normal_output, config.param('rawmpileup', 'mpileup_other_options'), sequence['name'],fast_bed),
                ], name="rawmpileup." + tumor_pair.name + "." + sequence['name'] + "." + fast_name + ".normal" if fast_name else "rawmpileup." + tumor_pair.name + "." + sequence['name'] + ".normal"))
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")], tumor_output, config.param('rawmpileup', 'mpileup_other_options'), sequence['name'],fast_bed),
                ], name="rawmpileup." + tumor_pair.name + "." + sequence['name'] + "." + fast_name + ".tumor" if fast_name else "rawmpileup." + tumor_pair.name + "." + sequence['name'] + ".tumor" ))

        return jobs
    
    def rawmpileup_fastpass(self):
        jobs=self.rawmpileup(fast_bed=self.fastpass_bed, fast_name="fast_pass")
        return jobs
    
    def rawmpileup_slowpass(self):
        jobs=self.rawmpileup()
        return jobs

    def rawmpileup_cat(self, fast_name=None):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
            mpileup_suffix =  fast_name + ".mpileup" if fast_name else "mpileup"
            job_suffix = "_" + fast_name if fast_name else "" 

            mpileup_normal_file_prefix = os.path.join(varscan_directory, tumor_pair.normal.name + ".")
            mpileup_normal_inputs = [mpileup_normal_file_prefix + sequence['name'] + "." +  mpileup_suffix for sequence in self.sequence_dictionary_variant()]

            mpileup_tumor_file_prefix = os.path.join(varscan_directory, tumor_pair.tumor.name + ".")
            mpileup_tumor_inputs = [mpileup_tumor_file_prefix + sequence['name'] + "." +  mpileup_suffix for sequence in self.sequence_dictionary_variant()]

            normal_output = mpileup_normal_file_prefix + mpileup_suffix 
            tumor_output = mpileup_tumor_file_prefix + mpileup_suffix

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),
                Job(mpileup_normal_inputs, [normal_output], command = "cat \\\n  " + " \\\n  ".join(mpileup_normal_inputs) + " \\\n  > " + normal_output),
            ], name = "rawmpileup_cat." + tumor_pair.name + job_suffix + "_normal"))
                Job(mpileup_tumor_inputs, [tumor_output], command = "cat \\\n  " + " \\\n  ".join(mpileup_tumor_inputs) + " \\\n  > " + tumor_output)
            ], name = "rawmpileup_cat." + tumor_pair.name + job_suffix + "_tumor"))
          
        return jobs

    def rawmpileup_cat_fastpass(self):
        jobs=self.rawmpileup_cat(fast_name="fast_pass")
        return jobs
    
    def rawmpileup_cat_lowpass(self):
        jobs=self.rawmpileup_cat()
        return jobs

    def paired_varscan2(self, fast_name=None):
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
            mpileup_suffix =  fast_name + ".mpileup" if fast_name else "mpileup"
            output_suffix = "." + fast_name if fast_name else ""
            job_suffix = "_" + fast_name if fast_name else ""

            for sequence in self.sequence_dictionary_variant():
                input_normal = os.path.join(varscan_directory, tumor_pair.normal.name + "." + sequence['name'] + "." + mpileup_suffix)
                input_tumor = os.path.join(varscan_directory, tumor_pair.tumor.name + "." + sequence['name'] + "." + mpileup_suffix)

                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + output_suffix)
                output_snp =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + output_suffix + ".snp.vcf")
                output_indel =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + output_suffix + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + output_suffix + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + output_suffix + ".varscan2.vcf.gz")
           
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory]),
                    varscan.somatic(input_normal, input_tumor, output, config.param('varscan2_somatic', 'other_options'), output_snp_dep=output_snp, output_indel_dep=output_indel),
                    htslib.bgzip_tabix_vcf(output_snp, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")),
                    htslib.bgzip_tabix_vcf(output_indel, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")),
                    pipe_jobs([
                        bcftools.concat([output_snp + ".gz",output_indel + ".gz")], None),
                        Job([None], [output_vcf], command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > " + output_vcf),
                    ]),
                    htslib.bgzip_tabix_vcf(output_vcf, output_vcf_gz),
                ], name = "varscan2_somatic." + tumor_pair.name + "." + sequence['name'] + job_suffix) )

        return jobs
   
    def paired_varscan2_fastpass(self):
        jobs=self.paired_varscan2(fast_name="fast_pass")
        return jobs
    
    def paired_varscan2_lowpass(self):
        jobs=self.paired_varscan2()
        return jobs
        

    def varscan2_fpfilter(self, fast_name=None):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 filtering thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        Somatic and germline calls are filter using bam-readcount info from tumor and normal bam, respectively. 
        Variants denoted as PASS and somatic (SS=1) or germline (SS=2) and loh (SS=3) are retained for further analysis.
        """

        jobs = []
        file_suffix = "." + fast_name if fast_name else ""
        job_suffix = "_" + fast_name if fast_name else ""
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
        
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")    

            for sequence in self.sequence_dictionary_variant():
                input_vcf =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.vcf")
                output_bed = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.bed")
                output_normal_readcount = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".normal.readcount")
                output_tumor_readcount =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".tumor.readcount")

                output_fpfilter_somatic =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.fpfilter.somatic.vcf")
                output_fpfilter_somatic_gz =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.fpfilter.somatic.vcf.gz")
                output_fpfilter_germline =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.fpfilter.germline_loh.vcf")
                output_fpfilter_germline_gz =  os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.fpfilter.germline_loh.vcf.gz") 
                output_vcf_somatic_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.somatic.vcf.gz")            
                output_vcf_germline_loh_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.germline_loh.vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline]),
                    tools.vcf2bed(input_vcf, output_bed),
                ], name = "varscan2_vcf2bed." + tumor_pair.name + "." + sequence['name'] + job_suffix ) )
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory,
                    removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline]),
                    bamreadcount.readcount(inputNormal, output_bed, output_normal_readcount),
                ], name = "varscan2_readcount." + tumor_pair.name + "." + sequence['name']+ job_suffix + "_normal" ) )
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory,
                    removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline]),
                    bamreadcount.readcount(inputTumor, output_bed, output_tumor_readcount),
                ], name = "varscan2_readcount." + tumor_pair.name + "." + sequence['name']+ job_suffix + "_tumor" ) )
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory,
                    removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline]),
                    varscan.fpfilter_somatic(input_vcf, output_tumor_readcount, output_fpfilter_somatic),
                    htslib.bgzip_tabix_vcf(output_fpfilter_somatic, output_fpfilter_somatic_gz),
                    pipe_jobs([
                        bcftools.view(output_fpfilter_somatic_gz, None, config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_vcf_somatic_gz),
                    ]),
                ], name = "varscan2_fpfilter." + tumor_pair.name + "." + sequence['name']+ job_suffix + "_tumor" ) )
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory,
                    removable_files=[varscan_directory,output_fpfilter_somatic,output_fpfilter_germline]),
                    varscan.fpfilter_somatic(input_vcf, output_normal_readcount, output_fpfilter_germline),
                    htslib.bgzip_tabix_vcf(output_fpfilter_germline, output_fpfilter_germline_gz),
                    pipe_jobs([
                        bcftools.view(output_fpfilter_germline_gz, None, config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_vcf_germline_loh_gz),
                    ]),             
                ], name = "varscan2_fpfilter." + tumor_pair.name + "." + sequence['name'] + job_suffix + "_normal" ) )

        return jobs
    
    def varscan2_fpfilter_fastpass(self):
        jobs=self.varscan2_fpfilter(fast_name="fast_pass")
        return jobs
    
    def varscan2_fpfilter_slowpass(self):
        jobs=self.varscan2_fpfilter()
        return jobs

    def merge_varscan2(self, fast_name=None):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        file_suffix = "." + fast_name if fast_name else ""
        job_suffix = "_" + fast_name if fast_name else ""
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.vcf.gz") for sequence in self.sequence_dictionary_variant()]    
        
            somatic_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.somatic.vcf.gz") for sequence in self.sequence_dictionary_variant()]

            germline_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + file_suffix + ".varscan2.germline_loh.vcf.gz") for sequence in self.sequence_dictionary_variant()]

            all_output = os.path.join(pair_directory, tumor_pair.name + file_suffix + ".varscan2.vcf.gz")
            somtic_output = os.path.join(pair_directory, tumor_pair.name + file_suffix + ".varscan2.somatic.vcf.gz")
            germline_output = os.path.join(pair_directory, tumor_pair.name + file_suffix + ".varscan2.germline_loh.vcf.gz")

            jobs.append(concat_jobs([
                pipe_jobs([
                    bcftools.concat(all_inputs, None),
                    htslib.bgzip_tabix_vcf(None, all_output),
                ]),
            ], name = "merge_varscan2." + tumor_pair.name + job_suffix + "_all"))
            jobs.append(concat_jobs([
                pipe_jobs([
                    bcftools.concat(somatic_inputs, None),
                    htslib.bgzip_tabix_vcf(None, somtic_output),
                ]),
            ], name = "merge_varscan2." + tumor_pair.name + job_suffix + "_somatic"))
            jobs.append(concat_jobs([
                 pipe_jobs([
                    bcftools.concat(germline_inputs, None),
                    htslib.bgzip_tabix_vcf(None, germline_output),                 
                ]),
            ], name = "merge_varscan2." + tumor_pair.name + job_suffix + "_germline"))

        return jobs
    
    def merge_varscan2_fastpass(self):
        jobs=self.merge_varscan2(fast_name="fast_pass")
        return jobs
    
    def merge_varscan2_slowpass(self):
        jobs=self.merge_varscan2()
        return jobs


    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + mutect_directory, removable_files=[mutect_directory])
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"))
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    outPrefix =  tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        mkdir_job,
                        gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, outPrefix + ".vcf.gz"), intervals=sequences)               
                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(inputNormal, inputTumor, os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"), exclude_intervals=unique_sequences_per_job_others)
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
            output = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_somatic =  os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vcf.gz")

            if nb_jobs == 1:
                input = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output], command="ln -s -f " + input + " " + output), 
                    pipe_jobs([
                        Job([output_gz], [None], command="zcat " + output + " | sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),
                        bcftools.view(None, None, config.param('gatk_cat_mutect', 'filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                ], name="symlink_mutect_vcf." + tumor_pair.name))
    
            elif nb_jobs > 1:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputVCFs = []
                for idx,sequences in enumerate(unique_sequences_per_job):
                    inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputVCFs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                jobs.append(concat_jobs([
                    gatk.cat_variants(inputVCFs, output),
                    pipe_jobs([
                        Job([output_gz], [None], command="zcat " + output + " | sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),
                        bcftools.view(None, None, config.param('gatk_cat_mutect', 'filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatic),
                     ]),
                 ], name = "gatk_cat_mutect." + tumor_pair.name))

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
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
 
            paired_sample = [inputNormal,inputTumor]
           
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory]),
                    pipe_jobs([
                        samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options')),
                        samtools.bcftools_view("-", os.path.join(samtools_directory,  tumor_pair.name + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                    ]),
                ], name="samtools_paired." + tumor_pair.name))

            else:
                for region in self.generate_approximate_windows(nb_jobs): #for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory]),
                        pipe_jobs([
                            samtools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'), region),
                            samtools.bcftools_view("-", os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf"), config.param('samtools_paired', 'bcftools_view_options'), pair_calling=True),
                        ]),
                ], name="samtools_paired." + tumor_pair.name + "." + region))

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
        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory,  tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory,  tumor_pair.name + ".samtools.vcf.gz")
            output_somatics = os.path.join(pair_directory,  tumor_pair.name + ".samtools.somatic.vcf.gz")
            output_germline_loh = os.path.join(pair_directory,  tumor_pair.name + ".samtools.germline_loh.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory,  tumor_pair.name + ".bcf") 
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf, None, config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf, None, config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory,  tumor_pair.name + "." + region + ".bcf") for region in self.generate_approximate_windows(nb_jobs)] #for idx,sequences in enumerate(unique_sequences_per_job):
                jobs.append(concat_jobs([
                    samtools.bcftools_cat(inputs, output),
                    pipe_jobs([
                        samtools.bcftools_view(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        htslib.bgzip_tabix_vcf(None, output_vcf),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf, None, config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_somatics),   
                    ]),
                     pipe_jobs([
                        bcftools.filter(output_vcf, None, config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_samtools." + tumor_pair.name))                    
        
        report_file = os.path.join("report", "DnaSeq.merge_filter_bcf.md")
        jobs.append(
            Job(
                [output_vcf],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/HumanVCFformatDescriptor.tsv \\
  report/ && \\
sed 's/\t/|/g' report/HumanVCFformatDescriptor.tsv | sed '2i-----|-----' >> {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_filter_bcf_report")
        )

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
        if use_bed :
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(bed)
            bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="./vardict.tmp")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + vardict_directory, removable_files=[vardict_directory])

            if use_bed :
                idx = 0
                for bf in bed_file_list:
                    output=os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, bf),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                            Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
                            htslib.bgzip_tabix_vcf(None, output),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
                if nb_jobs == 1:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    output =os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        bedJob,
                        pipe_jobs([
                            vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds.pop()),
                            vardict.testsomatic(None,None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None ),
                            Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
                            htslib.bgzip_tabix_vcf(None, output),
                        ]),
                     ], name="vardict_paired." + tumor_pair.name + ".0"))
                else:
                    bedJob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([mkdir_job,bedJob], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                        jobs.append(concat_jobs([
                            mkdir_job,
                            pipe_jobs([
                                vardict.paired_java(inputNormal, inputTumor, tumor_pair.name, None, beds[idx]),
                                vardict.testsomatic(None, None),
                                vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                                Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, 'N', $4) } {print}'"),
                                htslib.bgzip_tabix_vcf(None, output), 
                            ]),
                        ], name="vardict_paired." + tumor_pair.name + "." + str(idx) ))                        
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
            output = os.path.join(pair_directory,  tumor_pair.name + ".vardict.vcf.gz")
            output_somatic = os.path.join(pair_directory,  tumor_pair.name + ".vardict.somatic.vcf.gz")
            output_germline_loh = os.path.join(pair_directory,  tumor_pair.name + ".vardict.germline_loh.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory,  tumor_pair.name + ".0.vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input], [output], command="ln -s -f " + input + " " + output),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatics),
                    ]),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], removable_files=[output], name="symlink_vardict_vcf." + tumor_pair.name ))
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(os.path.join(vardict_directory,  tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))
                
                jobs.append(concat_jobs([
                    #gatk.cat_variants(inputVCFs, output),
                    pipe_jobs([
                        bcftools.concat(inputVCFs, None),
                        htslib.bgzip_tabix_vcf(None, output),
                    ]),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_somatic),
                    ]),
                    pipe_jobs([
                       bcftools.view(output, None, config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                       htslib.bgzip_tabix_vcf(None, output_germline_loh),
                    ]),
                ], name = "merge_filter_paired_vardict." + tumor_pair.name ))

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

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.somatic.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.somatic.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")        
            inputSomaticVCFs = [input_mutect2, input_vardict, input_samtools, input_varscan2]
                        
        
            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vcf")
            output_gz = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vcf.gz")
            output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic-work"), output_ensemble] )   
       
            jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    bcbio_variation.tumor_pair_ensemble(inputSomaticVCFs, output_ensemble, config.param('bcbio_ensemble_somatic', 'config_yaml')),
                    htslib.bgzip_tabix_vcf(output_ensemble, output_gz),
                    pipe_jobs([
                        Job([output_gz], [None], command="zgrep -Pv \"set=samtools\\t\" " + output_gz + " | grep -Pv \"set=vardict\\t\" " + " | grep -Pv \"set=varscan2\\t\" "),
                        htslib.bgzip_tabix_vcf(None, output_flt),
                    ]),
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

            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline_loh.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.germline_loh.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline_loh.vcf.gz")
            inputGermlineVCFs = [input_vardict, input_samtools, input_varscan2]

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.vcf")
            output_gz = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.vcf.gz")
            output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.flt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh-work"), output_ensemble])

            jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    bcbio_variation.tumor_pair_ensemble(inputGermlineVCFs, output_ensemble, config.param('bcbio_ensemble_germline_loh', 'config_yaml')),
                    htslib.bgzip_tabix_vcf(output_ensemble, output_gz),
                    pipe_jobs([
                        Job([output_gz], [None], command="zgrep -Pv \"set=samtools\\t\" " + output_gz + " | grep -Pv \"set=vardict\\t\" " + " | grep -Pv \"set=varscan2\\t\" "),
                        htslib.bgzip_tabix_vcf(None, output_flt),
                    ]),
                ], name="bcbio_ensemble_germline_loh." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")
            output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.flt.annot.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_somatic_variants])

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.variant_annotator( inputNormal, inputTumor, input_somatic_variants, output_somatic_variants ),
            ], name = "gatk_variant_annotator.somatic." + tumor_pair.name))            
        
        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputNormal = os.path.join("alignment",tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment",tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.flt.vcf.gz")
            output_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline_loh.flt.annot.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_germline_loh_variants])

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.variant_annotator( inputNormal, inputTumor, input_germline_loh_variants, output_germline_loh_variants ),
            ], name = "gatk_variant_annotator.germline." + tumor_pair.name))

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

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.flt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.flt.annot.snpeff.vcf")         
            output_somatic_gz = os.path.join(paired_directory, tumor_pair.name + ".ensemble.annot.snpeff.vcf.gz")              

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write( tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n" )
            
            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_somatic])

            jobs.append(concat_jobs([
                mkdir_job,           
                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename, options=config.param('compute_cancer_effects_somatic', 'options')),
                htslib.bgzip_tabix_vcf(output_somatic, output_somatic_gz),
            ],name = "compute_cancer_effects_somatic." + tumor_pair.name))

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

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.flt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.flt.annot.snpeff.vcf")
            output_germline_gz = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write( tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n" )

            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_germline])

            jobs.append(concat_jobs([
                mkdir_job,
                snpeff.compute_effects(input_germline, output_germline, options=config.param('compute_cancer_effects_germline', 'options')),
                htslib.bgzip_tabix_vcf(output_germline, output_germline_gz),
            ],name = "compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []
        
        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name , tumor_pair.name + ".ensemble.somatic.flt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output], command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ],name="gatk_combine_variants.somatic.allPairs"))
 
        else:
            
            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.somatic.allPairs"))

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name , tumor_pair.name + ".ensemble.germline_loh.flt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output], command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ],name="gatk_combine_variants.germline_loh.allPairs"))

        else:

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.germline_loh.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vt.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])
    
        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input, output)
        ],name = "decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vt.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input, output)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])

        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input, output)
        ],name = "decompose_and_normalize_mnps.somatic.allPairs"))

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
        input =  os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vt.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vt.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vt.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        for tumor_pair in self.tumor_pairs.itervalues():
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")
        
        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])
    
        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, cancer_sample_file=cancer_pair_filename,options=config.param('compute_cancer_effects_somatic', 'options')),
            htslib.bgzip_tabix_vcf(output, output_gz),
        ],name = "compute_effects.somatic.allPairs"))

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
        input =  os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vt.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vt.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vt.snpeff.vcf.gz")      

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output])

        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, options=config.param('compute_cancer_effects_germline', 'options')),
            htslib.bgzip_tabix_vcf(output, output_gz),
        ],name = "compute_effects.germline.allPair"))

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []
        
        ensemble_directory = os.path.join("pairedVariants", "ensemble")        
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")
        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory),
            gemini.gemini_annotations( gemini_prefix + ".ensemble.somatic.annot.vt.snpeff.vcf.gz", gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
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
        gemini_module=config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2],gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory),
            gemini.gemini_annotations( gemini_prefix + ".ensemble.germline_loh.annot.vt.snpeff.vcf.gz", gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.germline.allPairs"))

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.bwa_mem_picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.fix_mate_by_coordinate,
            self.picard_mark_duplicates,
            self.recalibration,
            self.metrics,
            self.picard_calculate_hs_metrics,
            self.gatk_callable_loci,
            self.extract_common_snp_freq,
            self.baf_plot,
            self.rawmpileup_fastpass,
            self.rawmpileup_cat_fastpass,
            self.paired_varscan2_fastpass,
            self.varscan2_fpfilter_fastpass,
            self.merge_varscan2_fastpass,
            self.paired_mutect2_fastpass,
            self.merge_mutect2_fastpass,
            self.samtools_paired_fastpass,
            self.merge_filter_paired_samtools_fastpass,
            self.vardict_paired_fastpass,
            self.merge_filter_paired_vardict_fastpass,
            self.ensemble_somatic_fastpass,
            self.gatk_variant_annotator_somatic_fastpass,
            self.compute_cancer_effects_somatic_fastpass,
            self.combine_tumor_pairs_somatic_fastpass,
            self.decompose_and_normalize_mnps_somatic_fastpass,
            self.all_pairs_compute_effects_somatic_fastpass,
            self.gemini_annotations_somatic_fastpass,
            self.ensemble_germline_loh_fastpass,
            self.gatk_variant_annotator_germline_fastpass,
            self.compute_cancer_effects_germline_fastpass,
            self.combine_tumor_pairs_germline_fastpass,
            self.decompose_and_normalize_mnps_germline_fastpass,
            self.all_pairs_compute_effects_germline_fastpass,
            self.gemini_annotations_germline_fastpass,
            self.rawmpileup_slowpass
            self.rawmpileup_cat_slowpass,
            self.paired_varscan2_slowpass,
            self.varscan2_fpfilter_slowpass,
            self.merge_varscan2_slowpass,
            self.paired_mutect2_slowpass,
            self.merge_mutect2_slowpass,
            self.samtools_paired_slowpass,
            self.merge_filter_paired_samtools_slowpass,
            self.vardict_paired_slowpass,
            self.merge_filter_paired_vardict_slowpass,
            self.ensemble_somatic_slowpass,
            self.gatk_variant_annotator_somatic_slowpass,
            self.compute_cancer_effects_somatic_slowpass,
            self.combine_tumor_pairs_somatic_slowpass,
            self.decompose_and_normalize_mnps_somatic_slowpass,
            self.all_pairs_compute_effects_somatic_slowpass,
            self.gemini_annotations_somatic_slowpass,
            self.ensemble_germline_loh_slowpass,
            self.gatk_variant_annotator_germline_slowpass,
            self.compute_cancer_effects_germline_slowpass,
            self.combine_tumor_pairs_germline_slowpass,
            self.decompose_and_normalize_mnps_germline_slowpass,
            self.all_pairs_compute_effects_germline_slowpass,
            self.gemini_annotations_germline_slowpass
        ]

if __name__ == '__main__': 
    profyle_dna()
