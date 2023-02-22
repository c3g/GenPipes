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
import os
import re
import sys
import itertools

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs
import utils.utils
from core.readset import parse_illumina_readset_file

from bfx import bvatools
from bfx import bismark
from bfx import picard2 as picard
from bfx import bedtools
from bfx import samtools
from bfx import sambamba
from bfx import gatk
from bfx import igvtools
from bfx import bissnp
from bfx import tools
from bfx import ucsc
from bfx import fgbio
from bfx import metrics
from bfx import bash_cmd as bash
from bfx import dragen

from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class MethylSeq(dnaseq.DnaSeqRaw):
    """
    Methyl-Seq Pipeline
    ================

The GenPIpes Methyl-Seq pipeline now has three protocols.
 1. bismark
 2. hybrid
 3. dragen

The "bismark" protocol uses Bismark to align reads to the reference genome.
Picard is used to mark and remove duplicates and generate metric files. The
"hybrid" protocl uses [Illumina Dragen Bio-IT processor](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)
and [dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm)
software to align the reads to the reference genome. All the other steps are common
with bismark protocol. The "dragen" protocol uses Dragen to align reads to the
reference genome, call methylation, mark and remove duplicates.

Although dragen provides higher rate of mapping percentage with in a very short time
duration (approximately three hours compared to 30 hours from bismark), it only
accessible through McGill Genome Center cluster Abacus and The jobs cannot be
submitted to any of the HPCs from the
[Digital Research Aliance](https://status.computecanada.ca/).
Importantly, the user needs to have permission to submit jobs to Abacus. Therefore,
other users may continue to use only bismark protocol since it works in all the clusters.

However, if you would like to setup and use dragen in own cluster please refer our
[GenPipes Documentation](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_wgs_methylseq.html)

    """

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_illumina_readset_file(self.args.readsets.name)
                for readset in self._readsets:
                    if readset._run == "":
                        _raise(SanitycheckError("Error: no run was provided for readset \"" + readset.name +
                            "\"... Run has to be provided for all the readsets in order to use this pipeline."))
                    if readset._lane == "":
                        _raise(SanitycheckError("Error: no lane provided for readset \"" + readset.name +
                            "\"... Lane has to be provided for all the readsets in order to use this pipeline."))
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=argparse.FileType('r'))
        self.argparser.add_argument("-t", "--type", help="MethylSeq analysis type", choices=["bismark", "hybrid", "dragen"],
                                    default="bismark")
        super(MethylSeq, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'ihec_metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'ihec_metrics'), self.output_dir),
            'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
            'methylation_call_directory': os.path.relpath(os.path.join(self.output_dir, 'methylation_call'), self.output_dir),
            'methylkit_directory': os.path.relpath(os.path.join(self.output_dir, 'methylkit'), self.output_dir),
            'corrected_umi_directory': os.path.relpath(os.path.join(self.output_dir, 'corrected_umi'), self.output_dir),
            'tracks_directory': os.path.relpath(os.path.join(self.output_dir, 'tracks'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
        }
        return dirs

    def bismark_align(self):
        """
        Align reads with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf).
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
            no_readgroup_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted_noRG.bam")
            output_bam = re.sub("_noRG.bam", ".bam", no_readgroup_bam)
            index_bam = re.sub("_noRG.bam", ".bam.bai", no_readgroup_bam)
            report_suffix = "_bismark_bt2_report.txt"

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                # Defining the bismark output files (bismark sets the names of its output files from the basename of fastq1)
                # Note : these files will then be renamed (using a "mv" command) to fit with the mugqic pipelines nomenclature (cf. no_readgroup_bam)
                bismark_out_bam = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_pe.bam", os.path.basename(fastq1)))
                bismark_out_report = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_PE_report.txt", os.path.basename(fastq1)))
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                # Defining the bismark output files (bismark sets the names of its output files from the basename of fastq1)
                # Note : these files will then be renamed (using a "mv" command) to fit with the mugqic pipelines nomenclature (cf. no_readgroup_bam)
                bismark_out_bam = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2.bam", os.path.basename(fastq1)))
                bismark_out_report = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_SE_report.txt", os.path.basename(fastq1)))
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + os.path.dirname(output_bam)),
                    bismark.align(
                        fastq1,
                        fastq2,
                        os.path.dirname(no_readgroup_bam),
                        [no_readgroup_bam, re.sub(".bam", report_suffix, no_readgroup_bam)],
                    ),
                    Job(command="mv " + bismark_out_bam + " " + no_readgroup_bam),
                    Job(command="mv " + bismark_out_report + " " + re.sub(".bam", report_suffix, no_readgroup_bam)),
                ], name="bismark_align." + readset.name, samples=[readset.sample])
            )
            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + alignment_directory),
                    picard.add_read_groups(
                        no_readgroup_bam,
                        output_bam,
                        readset.name,
                        readset.library if readset.library else readset.sample.name,
                        readset.run + "_" + readset.lane,
                        readset.sample.name
                    )
                ], name="picard_add_read_groups." + readset.name, samples=[readset.sample])
            )
            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + alignment_directory),
                    sambamba.index(
                        output_bam,
                        index_bam
                    )
                ], name="sambamba_index." + readset.name, samples=[readset.sample])
            )
            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + alignment_directory),
                    sambamba.flagstat(
                        output_bam,
                        re.sub(".bam", "_flagstat.txt", output_bam),
                        config.param('sambamba_flagstat', 'flagstat_options')
                    )
                ], name="sambamba_flagstat." + readset.name, samples=[readset.sample])
            )

        report_file = os.path.join(self.output_dirs["report_directory"], "MethylSeq.bismark_align.md")
        jobs.append(
            Job(
                [os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name, readset.name + ".sorted.bam") for readset in self.readsets],
                [report_file],
                [['bismark_align', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable scientific_name="{scientific_name}" \\
  --variable assembly="{assembly}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    scientific_name=config.param('bismark_align', 'scientific_name'),
                    assembly=config.param('bismark_align', 'assembly'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="bismark_align_report")
        )

        return jobs


    def add_bam_umi(self):
        """
        Add read UMI tag to individual bam files using [fgbio](https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html).
        """

        jobs = []
        for readset in self.readsets:
            if readset.umi:
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
                input_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
                output_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.UMI.bam")
                output_bai = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.UMI.bai")
                input_umi = readset.umi
                input_umi_corrected = os.path.join(self.output_dirs["corrected_umi_directory"], readset.name, readset.name + ".corrected.fastq.gz")

                #correct umi fastq name (removing space)

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.join(self.output_dirs["corrected_umi_directory"], readset.name)),
                            fgbio.correct_readname(
                                input_umi,
                                input_umi_corrected
                            ),
                    ], name="fgbio_correct_readname." + readset.name, samples=[readset.sample])
                )

                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + os.path.dirname(output_bam)),
                        fgbio.addumi(
                            input_bam,
                            input_umi_corrected,
                            output_bam,
                            output_bai
                        ),
                    ], name="fgbio_addumi." + readset.name, samples=[readset.sample])
                )

        #report_file = os.path.join(self.output_dirs["report_directory"], "MethylSeq.bismark_align.md")
        #jobs.append(
            #Job(
                #[os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name, readset.name + ".sorted.bam") for readset in self.readsets],
                #[report_file],
                #[['bismark_align', 'module_pandoc']],
                #command="""\
#mkdir -p report && \\
#pandoc --to=markdown \\
  #--template {report_template_dir}/{basename_report_file} \\
  #--variable scientific_name="{scientific_name}" \\
  #--variable assembly="{assembly}" \\
  #{report_template_dir}/{basename_report_file} \\
  #> {report_file}""".format(
                    #scientific_name=config.param('bismark_align', 'scientific_name'),
                    #assembly=config.param('bismark_align', 'assembly'),
                    #report_template_dir=self.report_template_dir,
                    #basename_report_file=os.path.basename(report_file),
                    #report_file=report_file
                #),
                #report_files=[report_file],
                #name="bismark_align_report")
        #)

        return jobs


    def picard_remove_duplicates(self):
        """
        Remove duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be removed as a duplicate in the BAM file. Removing duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        This step is used in bismark and hybrid protocols.
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".")
            input = alignment_file_prefix + "sorted.bam"
            bam_output = alignment_file_prefix + "sorted.dedup.bam"
            metrics_file = alignment_file_prefix + "sorted.dedup.metrics"

            job = picard.mark_duplicates(
                [input],
                bam_output,
                metrics_file,
                remove_duplicates="true",
                ini_section='mark_duplicates'
            )
            job.name = "mark_duplicates." + sample.name
            job.samples = [sample]
            jobs.append(job)

            job = samtools.flagstat(
                bam_output,
                re.sub(".bam", "_flagstat.txt", bam_output),
            )
            job.name = "samtools_flagstat_dedup." + sample.name
            job.samples = [sample]
            jobs.append(job)

        report_file = os.path.join(self.output_dirs["report_directory"], "MethylSeq.picard_remove_duplicates.md")
        jobs.append(
              Job(
                [os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.bam") for sample in self.samples],
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
                name="picard_remove_duplicates_report")
        )

        return jobs


    def metrics(self):

        """
        Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
        Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
        Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
        covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
        whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
        bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
        for easy visualization of coverage in the IGV browser.

        """

        # check the library status
        library, bam = {}, {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample]="SINGLE_END"
            if readset.run_type == "PAIRED_END" :
                library[readset.sample]="PAIRED_END"
            if not readset.sample in bam:
                bam[readset.sample]=""
            if readset.bam:
                bam[readset.sample]=readset.bam

        jobs = []
        created_interval_lists = []
        for sample in self.samples:
            file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.")
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            candidate_input_files = [[file_prefix + "bam"]]
            if bam[sample]:
                candidate_input_files.append([bam[sample]])
            [input] = self.select_input_files(candidate_input_files)

            job = picard.collect_multiple_metrics(
                input,
                re.sub("bam", "all.metrics", input),
                library_type=library[sample]
            )
            job.name = "picard_collect_multiple_metrics." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Compute genome coverage with GATK
            job = gatk.depth_of_coverage(
                input,
                re.sub("bam", "all.coverage", input),
                coverage_bed
            )
            job.name = "gatk_depth_of_coverage.genome." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Compute genome or target coverage with BVATools
            job = bvatools.depth_of_coverage(
                input,
                re.sub("bam", "coverage.tsv", input),
                coverage_bed,
                other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            )
            job.name = "bvatools_depth_of_coverage." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Get reads# raw and after duplication
            dedup_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.bam")
            raw_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.bam")

            count_dedup_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".dedup.count")
            count_raw_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".raw.count")


            job = samtools.mapped_count(
                raw_bam,
                count_raw_output,
                None
            )
            job.name = "count.raw." + sample.name
            job.removable_files=[count_dedup_output,count_raw_output]
            job.samples = [sample]
            jobs.append(job)

            job = samtools.mapped_count(
                dedup_bam,
                count_dedup_output,
                None
            )
            job.name = "count.dedup." + sample.name
            job.removable_files=[count_dedup_output,count_raw_output]
            job.samples = [sample]
            jobs.append(job)

            if coverage_bed:
                # Get on-target reads (if on-target context is detected) raw and after duplication

                count_dedup_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".onTarget.dedup.count")
                count_raw_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".onTarget.raw.count")
                job = samtools.mapped_count(
                    raw_bam,
                    count_raw_output,
                    coverage_bed
                )
                job.name = "ontarget_count.raw." + sample.name
                job.removable_files=[count_dedup_output,count_raw_output]
                job.samples = [sample]
                jobs.append(job)

                job = samtools.mapped_count(
                    dedup_bam,
                    count_dedup_output,
                    coverage_bed
                )
                job.name = "ontarget_count.dedup." + sample.name
                job.removable_files=[count_dedup_output,count_raw_output]
                job.samples = [sample]
                jobs.append(job)

                # Compute on target percent of hybridisation based capture
                interval_list = re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(
                        coverage_bed,
                        interval_list
                     )
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)
                file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.")
                job = picard.calculate_hs_metrics(
                    file_prefix + "bam",
                    file_prefix + "onTarget.tsv",
                    interval_list
                )
                job.name = "picard_calculate_hs_metrics." + sample.name
                job.samples = [sample]
                jobs.append(job)

            # Calculate the number of reads with higher mapping quality than the threshold passed in the ini file
            job = concat_jobs([
                samtools.view(
                    input,
                    re.sub(".bam", ".filtered_reads.counts.txt", input),
                    "-c " + config.param('mapping_quality_filter', 'quality_threshold')
                )
            ])
            job.name = "mapping_quality_filter." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Calculate GC bias
            gc_content_file = re.sub(".bam", ".gc_cov.1M.txt", input)
            job = bedtools.coverage(
                input,
                gc_content_file
            )
            # For captured analysis
            if coverage_bed:
                gc_content_on_target_file = re.sub(".bam", ".gc_cov.1M.on_target.txt", input)
                gc_content_target_job = bedtools.intersect(
                    gc_content_file,
                    gc_content_on_target_file,
                    coverage_bed
                )
                job = concat_jobs([
                    job,
                    gc_content_target_job
                ])
                gc_content_file = gc_content_on_target_file
            job = concat_jobs([
                job,
                metrics.gc_bias(
                    gc_content_file,
                    re.sub(".bam", ".GCBias_all.txt", input)
                )
            ])
            job.name = "GC_bias." + sample.name
            job.samples = [sample]
            jobs.append(job)

            job = igvtools.compute_tdf(input, input + ".tdf")
            job.name = "igvtools_compute_tdf." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def methylation_call(self):
        """
        The script reads in a bisulfite read alignment file produced by the Bismark bisulfite mapper
        and extracts the methylation information for individual cytosines.
        The methylation extractor outputs result files for cytosines in CpG, CHG and CHH context.
        It also outputs bedGraph, a coverage file from positional methylation data and cytosine methylation report.
        """

        # Check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample]="SINGLE_END"
            if readset.run_type == "PAIRED_END" :
                library[readset.sample]="PAIRED_END"

        jobs = []

        methylseq_protocol = self.args.type

        methylation_protocol = config.param('dragen_align', 'methylation_protocol', param_type='string')
        mapping_implementation = config.param('dragen_align', 'mapping_implementation', param_type='string')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            input_file = os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")

            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            outputs = [
                os.path.join(methyl_directory, "CpG_context_" + re.sub( ".bam", ".txt.gz", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", ".bedGraph.gz", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", ".CpG_report.txt.gz", os.path.basename(input_file)))
            ]

            if (methylseq_protocol == "hybrid" and methylation_protocol == "directional" and mapping_implementation=="single-pass"):
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + methyl_directory),
                        samtools.view(
                            input_file,
                            re.sub("sorted","sorted.filtered", input_file),
                            "-d XM -b"
                            ),
                        picard.sort_sam(
                            re.sub("sorted", "sorted.filtered", input_file),
                            re.sub("sorted", "readset_sorted", input_file),
                            "queryname"
                            )
                    ], name="picard_sort_sam." + sample.name, samples=[sample])
                )
                outputs = [re.sub("sorted", "readset_sorted", output) for output in outputs]
                bismark_job = bismark.methyl_call(
                    re.sub("sorted", "readset_sorted", input_file),
                    outputs,
                    library[sample]
                )
                bismark_job.name = "bismark_methyl_call." + sample.name
                bismark_job.samples = [sample]
                jobs.append( bismark_job )

            else:
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + methyl_directory),
                        picard.sort_sam(
                            input_file,
                            re.sub("sorted", "readset_sorted", input_file),
                            "queryname"
                            )
                    ], name="picard_sort_sam." + sample.name, samples=[sample])
                )
                outputs = [re.sub("sorted", "readset_sorted", output) for output in outputs]
                bismark_job = bismark.methyl_call(
                    re.sub("sorted", "readset_sorted", input_file),
                    outputs,
                    library[sample]
                )
                bismark_job.name = "bismark_methyl_call." + sample.name
                bismark_job.samples = [sample]
                jobs.append( bismark_job )
    
        return jobs

    def wiggle_tracks(self):
        """
        Generate wiggle tracks suitable for multiple browsers, to show coverage and methylation data.
        When using GRCh37 build of Human genome, to be compatible with the UCSC Genome Browser we only keep chromosomes 1-22, X, Y and MT,
        and we also rename them by prefixing "chr" to the chromosome anme (e.g. "1" becomes "chr1"), and changing the mitocondrial chromosome from "MT" to "chrM", and keeping the GRCh37 coordinates.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            # Generation of a bedGraph and a bigWig track to show the genome coverage
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")
            bed_graph_prefix = os.path.join(self.output_dirs["tracks_directory"], sample.name, sample.name)
            big_wig_prefix = os.path.join(self.output_dirs["tracks_directory"], "bigWig", sample.name)
            bed_graph_output = bed_graph_prefix + ".bedGraph"
            big_wig_output = big_wig_prefix + ".coverage.bw"

            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + os.path.join(self.output_dirs["tracks_directory"], sample.name) + " " + os.path.join(self.output_dirs["tracks_directory"], "bigWig"), removable_files=["tracks"]),
                    bedtools.graph(
                        input_bam,
                        bed_graph_output,
                        ""
                    )
                ], name="bed_graph." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output)), samples=[sample])
            )
            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + os.path.join(self.output_dirs["tracks_directory"], "bigWig")),
                    ucsc.bedGraphToBigWig(
                        bed_graph_output,
                        big_wig_output,
                        False
                    )
                ], name="wiggle." + re.sub(".bw", "", os.path.basename(big_wig_output)), samples=[sample])
            )

            # Generation of a bigWig from the methylation bedGraph
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_bed_graph = os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.bedGraph.gz")
            output_wiggle = os.path.join(self.output_dirs["tracks_directory"], "bigWig", sample.name + ".methylation.bw")

            jobs.append(
                concat_jobs([
                    Job(command="mkdir -p " + methyl_directory),
                    ucsc.bedGraphToBigWig(
                        input_bed_graph,
                        output_wiggle
                    )
                ], name = "bismark_bigWig." + sample.name, samples=[sample])
            )

        return jobs

    def methylation_profile(self):
        """
        Generation of a CpG methylation profile by combining both forward and reverse strand Cs.
        Also generating of all the methylatoin metrics : CpG stats, pUC19 CpG stats, lambda conversion rate, median CpG coverage, GC bias.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            cpg_input_file = os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.CpG_report.txt.gz")
            cpg_profile = re.sub(".CpG_report.txt.gz", ".CpG_profile.strand.combined.csv", cpg_input_file)

            # Generate CpG methylation profile
            job = tools.bismark_combine(
                cpg_input_file,
                cpg_profile
            )
            job.name = "methylation_profile." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Generate stats for lambda, pUC19 and regular CpGs
            cg_stats_output = re.sub(".CpG_report.txt.gz", ".profile.cgstats.txt", cpg_input_file)
            lambda_stats_output = re.sub(".CpG_report.txt.gz", ".profile.lambda.conversion.rate.tsv", cpg_input_file)
            puc19_stats_output = re.sub(".CpG_report.txt.gz", ".profile.pUC19.txt", cpg_input_file)
            job = tools.cpg_stats(
                cpg_profile,
                cg_stats_output,
                lambda_stats_output,
                puc19_stats_output
            )
            job.name = "CpG_stats." + sample.name
            job.samples = [sample]
            jobs.append(job)

            target_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if target_bed:
                # Create targeted combined file
                target_cpg_profile = re.sub("combined", "combined.on_target", cpg_profile)
                job = bedtools.intersect(
                    cpg_profile,
                    target_cpg_profile,
                    target_bed,
                    include_header=True
                )
                job.name = "extract_target_CpG_profile." + sample.name
                job.samples = [sample]
                jobs.append(job)
                cpg_profile = target_cpg_profile

                target_cpg_profile_count=os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.CpG_profile.strand.combined.on_target.count")
                job = metrics.target_cpg_profile(
                    target_cpg_profile,
                    target_cpg_profile_count,
                    sample.name
                )
                job.name = "metrics_target_CpG_profile." + sample.name
                job.samples = [sample]
                jobs.append(job)



            # Caluculate median & mean CpG coverage
            median_CpG_coverage = re.sub(".CpG_report.txt.gz", ".median_CpG_coverage.txt", cpg_input_file)
            job = tools.cpg_cov_stats(
                cpg_profile,
                median_CpG_coverage
            )
            job.name = "median_CpG_coverage." + sample.name
            job.samples = [sample]
            if target_bed:
                job.removable_files = [target_cpg_profile]
            jobs.append(job)

        return jobs

    def ihec_sample_metrics_report(self):
        """
        Retrieve the computed metrics which fit the IHEC standards and build a tsv report table for IHEC.
        Note: The dragen protocol does not generate a metric for estimated library size. You will have to run Picard separately for this metric.
        """

        jobs = []

        target_bed = bvatools.resolve_readset_coverage_bed(self.samples[0].readsets[0])
        metrics_all_file = os.path.join(self.output_dirs["metrics_directory"], "IHEC.sampleMetrics.stats")
        report_metrics_file = os.path.join(self.output_dirs["report_directory"], "IHEC.sampleMetricsTable.tsv")

        if target_bed:
            report_file = os.path.join(self.output_dirs["report_directory"], "MethylSeq.ihec_sample_metrics_targeted_report.md")
        else:
            report_file = os.path.join(self.output_dirs["report_directory"], "MethylSeq.ihec_sample_metrics_report.md")

        # Create the list of input files to handle job dependencies
        sample_list = []
        counter=0
        metrics_output_list=[]
        metrics_output_list.append(None)
        for sample in self.samples:
            inputs = []
            sample_list.append(sample.name)
            metrics_file =  os.path.join(self.output_dirs["ihec_metrics_directory"], sample.name + ".read_stats.txt")

            metrics_output_list.append(metrics_file)
            #add dependency for previous job in order to fill the allSample output correctly
            inputs.append(metrics_output_list[counter])

            # Trim log files
            for readset in sample.readsets:
                inputs.append((os.path.join(self.output_dirs["trim_directory"], sample.name, readset.name + ".trim.log")))

            # Aligned pre-deduplicated bam files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.bam"))

            # Deduplicated bam files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.bam"))

            # Coverage summary files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.all.coverage.sample_summary"))

            # Filtered reads count files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.filtered_reads.counts.txt"))

            # GC bias files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.GCBias_all.txt"))

            if self.protocol == "bismark":
                # Bismark alignment files
                for readset in sample.readsets:
                    inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, readset.name, readset.name + ".sorted_noRG_bismark_bt2_report.txt"))

            # CpG coverage files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".readset_sorted.dedup.median_CpG_coverage.txt"))

            # pUC19 methylation files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".readset_sorted.dedup.profile.pUC19.txt"))

            # Lambda conversion rate files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".readset_sorted.dedup.profile.lambda.conversion.rate.tsv"))

            # CG stat files
            [cgstats_file] = self.select_input_files([
                [os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".sorted.dedup.profile.cgstats.txt")],
                [os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".readset_sorted.dedup.profile.cgstats.txt")]
            ])
            inputs.append(cgstats_file)

            #Estimated library sizes
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".sorted.dedup.metrics"))

            # read count (raw, dedup)
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".dedup.count"))
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".raw.count"))

            # read count (raw, dedup)  and CpG_profyle in targeted context
            if target_bed :
                inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".onTarget.dedup.count"))
                inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".onTarget.raw.count"))
                inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".readset_sorted.dedup.CpG_profile.strand.combined.on_target.count"))


            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(self.output_dirs["ihec_metrics_directory"]),
                        bash.mkdir(self.output_dirs["metrics_directory"]),
                        tools.methylseq_ihec_metrics_report(
                            sample.name,
                            inputs,
                            metrics_file,
                            metrics_all_file,
                            target_bed,
                            counter
                        ),
                    ],
                    name="ihec_sample_metrics." + sample.name,
                    samples=[sample]
                )
            )
            counter+=1


        jobs.append(
            concat_jobs([
                Job(command="mkdir -p metrics"),
                Job(
                    [metrics_all_file],
                    [report_file],
                    [['ihec_sample_metrics_report', 'module_pandoc']],
                    command="""\
mkdir -p report && \\
cp {metrics_all_file} {report_metrics_file} && \\
metrics_table_md=`sed 's/\t/|/g' {metrics_file}`
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable sequence_alignment_table="$metrics_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                        report_template_dir=self.report_template_dir,
                        metrics_all_file=metrics_all_file,
                        metrics_file=metrics_file,
                        basename_report_file=os.path.basename(report_file),
                        report_metrics_file=report_metrics_file,
                        report_file=report_file
                    ),
                    report_files=[report_file]
                )
            ], name="ihec_sample_metrics_report")
        )

        return jobs

    def bis_snp(self):
        """
        SNP calling with [BisSNP](https://people.csail.mit.edu/dnaase/bissnp2011/).
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            input_file = os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")

            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            cpg_output_file = os.path.join(variant_directory, sample.name + ".cpg.vcf")
            snp_output_file = os.path.join(variant_directory, sample.name + ".snp.vcf")

            jobs.append(
                concat_jobs([
                    bash.mkdir(variant_directory),
                    bissnp.bisulfite_genotyper(
                        input_file,
                        cpg_output_file,
                        snp_output_file
                    )
                ], name="bissnp." + sample.name, samples=[sample])
            )

        return jobs

    def filter_snp_cpg(self):
        """
        SNP CpGs filtering.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_file = os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.CpG_profile.strand.combined.csv")

            output_file = re.sub("CpG_profile.strand.combined.csv", "filtered.bed", input_file)

            jobs.append(
                concat_jobs([
                    bash.mkdir(methyl_directory),
                    tools.filter_snp_cpg(
                        input_file,
                        output_file
                    )
                ], name="filter_snp_cpg." + sample.name, samples=[sample])
            )
        return jobs

    def prepare_methylkit(self):
        """
        Prepare input file for [methylKit](https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) differential analysis.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_file = os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.filtered.bed")

            output_directory = os.path.join(self.output_dirs["methylkit_directory"], "inputs")
            output_file = os.path.join(output_directory, re.sub("filtered.bed", "map.input", os.path.basename(input_file)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(output_directory),
                    tools.prepare_methylkit(
                        input_file,
                        output_file
                    )
                ], name="prepare_methylkit." + sample.name, samples=[sample])
            )
        return jobs

    def methylkit_differential_analysis(self):
        """
        Run methylKit to get DMCs & DMRs for different design comparisons.
        """

        jobs = []
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)

        input_directory = os.path.join(self.output_dirs["methylkit_directory"], "inputs")
        input_files = []
        for contrast in self.contrasts:
            input_files.extend([[os.path.join(input_directory, sample.name + ".readset_sorted.dedup.map.input") for sample in group] for group in [contrast.controls, contrast.treatments]])

        input_files = list(itertools.chain.from_iterable(input_files))

        output_directory = os.path.join(self.output_dirs["methylkit_directory"], "results")
        output_files = [os.path.join(output_directory, "Rdata_files", contrast.name, "perbase.testingresults.txt.gz") for contrast in self.contrasts]

        methylkit_job = tools.methylkit_differential_analysis(
            design_file,
            input_files,
            output_files,
            output_directory
        )

        return [concat_jobs([
            bash.mkdir(output_directory),
            methylkit_job
        ], name="methylkit_differential_analysis")]

    def dragen_align(self):
        """
        Align reads with [Dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm) both hybrid and dragen protocols use this step to align reads.
        The Dragen parameters can be changed using other_options of the ini configuration.
        """
       # duplicate_marking = config.param('dragen_align', 'duplicate_marking', param_type='string').lower()

        jobs = []
        methylseq_protocol = self.args.type

        methylation_protocol = config.param('dragen_align', 'methylation_protocol', param_type='string')
        mapping_implementation = config.param('dragen_align', 'mapping_implementation', param_type='string')
        # if the protocol is hybrid and methylation_protocol is "directional-complement" and
        # mapping_implementation is "sigle-pass" pipeline will not be generating genpipes file. Dragen protocol
        # should be used in this case.
        # In order to run the hybrid protocol for directional libraries and generate dragen outfiles that are compatible with bismark methylation calling,
        # the following should be set in the config file for dragen_align:
        # methylation_protocol=directional
        # OR, for multi-pass mapping_implementation:
        # methylation_protocol=directional, sort=false, mapping_implementation=multi-pass, duplicate_marking=false, remove_duplicates=false.
        # note that the "multi-pass" mapping implementation is deprecated and no longer recommended by Illumina. 
        if not (methylseq_protocol == "hybrid" and methylation_protocol == "directional-complement" and mapping_implementation=="single-pass"):

            for readset in self.readsets:
                trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
                # All the input files should copy to dragen work folder because IO operations with default locations
                # are really slow or have permission issues. After processing files all the files should move back to
                # the default GenPipes folder and remove files from dragen work folder
                dragen_inputfolder = os.path.join(config.param('dragen_align', 'work_folder'), "reads", readset.name)
                dragen_workfolder = os.path.join(config.param('dragen_align', 'work_folder'), "alignment", readset.name)
                dragen_tmp_bam = os.path.join(dragen_workfolder, readset.name + ".bam")
                # dragen output file name
                dragen_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
                #output_bam = dragen_bam
                index_bam = dragen_bam + ".bai"

                # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
                if readset.run_type == "PAIRED_END":
                    candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                    if readset.fastq1 and readset.fastq2:
                         candidate_input_files.append([readset.fastq1, readset.fastq2])
                    if readset.bam:
                        candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                      re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                    [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                    # first fastqs need to be copy to the dragen work_folder
                    dragen_tmp_fastq1 = os.path.join(dragen_inputfolder, readset.name + ".pair1.fastq.gz")
                    dragen_tmp_fastq2 = os.path.join(dragen_inputfolder, readset.name + ".pair2.fastq.gz")
                    cp_dragen_fastq_job = concat_jobs([
                        bash.mkdir(dragen_inputfolder),
                        bash.cp(os.path.abspath(fastq1), dragen_tmp_fastq1),
                        bash.cp(os.path.abspath(fastq2), dragen_tmp_fastq2)
                    ], name="dragen_copy_fastq." + readset.name, samples=[readset.sample])
                    rm_dragen_fastq_job = concat_jobs([
                        bash.rm(dragen_tmp_fastq1, recursive=True, force=True),
                        bash.rm(dragen_tmp_fastq2, recursive=True, force=True)
                    ], name="dragen_remove_fastq." + readset.name, samples=[readset.sample])
                elif readset.run_type == "SINGLE_END":
                    candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                    if readset.fastq1:
                        candidate_input_files.append([readset.fastq1])
                    if readset.bam:
                        candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                    [fastq1] = self.select_input_files(candidate_input_files)
                    fastq2 = None
                    # first fastqs need to be copy to the dragen work_folder
                    dragen_tmp_fastq1 = os.path.join(dragen_workfolder, readset.name + ".single.fastq.gz")
                    dragen_tmp_fastq2 = fastq2
                    cp_dragen_fastq_job = concat_jobs([
                        bash.mkdir(dragen_inputfolder ),
                        bash.cp(os.path.abspath(fastq1), dragen_tmp_fastq1)
                    ], name="dragen_copy_fastq." + readset.name, samples=[readset.sample])
                    rm_dragen_fastq_job = concat_jobs([
                        bash.rm(dragen_tmp_fastq1, recursive=True, force=True)
                    ], name="dragen_remove_fastq." + readset.name, samples=[readset.sample])
                else:
                    _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                                            "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
                dragen_align = concat_jobs([
                    #to track the dragen output files and done file input dependency and output dependency of the original
                    #location should be defined
                    cp_dragen_fastq_job,
                    bash.mkdir(dragen_workfolder),
                    bash.mkdir(os.path.join(config.param('dragen_align', 'work_folder'), "job_output", "dragen_align")),
                    bash.mkdir(os.path.abspath(alignment_directory)),
                    dragen.align_methylation(
                        dragen_tmp_fastq1,
                        dragen_tmp_fastq2,
                        dragen_workfolder,
                        readset.name,
                        readset.sample.name,
                        readset.library if readset.library else readset.sample.name,
                        readset.name + "_" + readset.run + "_" + readset.lane, protocol=methylseq_protocol

                    ),
                    bash.cp(dragen_workfolder, os.path.abspath(alignment_directory) + "/", recursive=True),
                    bash.rm(dragen_workfolder, recursive=True, force=True),
                    # bash.rm(dragen_workfolder, recursive=True, force=True, input_dependency=[fastq1,fastq2], output_dependency=[dragen_bam]),
                    rm_dragen_fastq_job
                ], name="dragen_align." + readset.name, samples=[readset.sample], input_dependency=[fastq1, fastq2], output_dependency=[dragen_bam])
                jobs.append(
                    dragen_align
                )
                jobs.append(
                    concat_jobs([
                        bash.mkdir(alignment_directory),
                        samtools.flagstat(
                            dragen_bam,
                            re.sub(".bam", "_flagstat.txt", dragen_bam),
                        )
                    ], name="samtools_flagstat." + readset.name, samples=[readset.sample], output_dependency = [re.sub(".bam", "_flagstat.txt", dragen_bam)])

                )
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + alignment_directory),
                        sambamba.index(
                            dragen_bam,
                            index_bam
                        )
                    ], name="sambamba_index." + readset.name, samples=[readset.sample]))

            return jobs

        else:
            _raise(SanitycheckError("Please use \"dragen\" protocol when using directional-complement on single-pass mode. Skipping generating the genpipes file..."))

    def dragen_methylation_call(self):
        """
        Call methylation with Dragen using the 2nd run of Dragen alignment.
        """
       # duplicate_marking = config.param('dragen_align', 'duplicate_marking', param_type='string').lower()

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            methylation_call_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            dragen_inputfolder = os.path.join(config.param('dragen_align', 'work_folder'), "reads", sample.name)
            dragen_workfolder = os.path.join(config.param('dragen_align', 'work_folder'), "dragen_methylation_call", sample.name)
            dragen_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            output_report = os.path.join(methylation_call_directory, sample.name + ".CX_report.txt")

            dragen_tmp_bam = os.path.join(dragen_inputfolder, sample.name + ".sorted.bam")
            cp_dragen_bam_job = concat_jobs([
                bash.mkdir(dragen_inputfolder),
                bash.cp(os.path.abspath(dragen_bam), dragen_tmp_bam)
            ], name="dragen_copy_bam." + sample.name, samples=[sample])
            rm_dragen_bam_job = concat_jobs([
                bash.rm(dragen_tmp_bam, recursive=True, force=True)
            ], name="dragen_remove_bam." + sample.name, samples=[sample])


           # candidate_input_files2 = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            dragen_methylationn_call = concat_jobs([
                cp_dragen_bam_job,
                bash.mkdir(dragen_workfolder),
                bash.mkdir(os.path.join(config.param('dragen_align', 'work_folder'), "job_output", "dragen_methylation_call")),
                bash.mkdir(os.path.abspath(methylation_call_directory)),
                dragen.call_methylation(
                    dragen_tmp_bam,
                    dragen_workfolder,
                    sample.name, output=output_report),
                bash.cp(dragen_workfolder, os.path.abspath(self.output_dirs["methylation_call_directory"]) + "/", recursive=True),
                bash.rm(dragen_workfolder, recursive=True, force=True ),
                rm_dragen_bam_job
            ], name="dragen_methylation_call." + sample.name, samples=[sample],  input_dependency=[dragen_bam],
                         output_dependency=[output_report])
            jobs.append(
                dragen_methylationn_call
            )
        return jobs

    def sort_dragen_sam(self):
        """
        Coordinate sorting the bam file resulted from dragen and create an index.
        """
        jobs = []

        duplicate_marking = config.param('dragen_align', 'duplicate_marking', param_type='boolean')

        if duplicate_marking == True:

            for sample in self.samples:
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

                input_file = os.path.join(alignment_directory, sample.name + ".sorted.bam")
                output_file = os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")
                output_file_index = os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam.bai")
                # empty metric file will be created for "ihec_sample_metrics_report" steps. Otherwise it will be failed.
                # Therefore added an empty file with ESTIMATED_LIBRARY_SIZE = NA
                # so the dragen protocol will not estimate the library size
                empty_metric_file= os.path.join(alignment_directory, sample.name + ".sorted.dedup.metrics")

                jobs.append(

                    concat_jobs([
                        bash.mkdir(alignment_directory),
                        picard.sort_sam(
                            input_file,
                            output_file
                        ), Job(output_files=[empty_metric_file], command="""\
                        printf "ESTIMATED_LIBRARY_SIZE\\nNA" > {output}""".
                            format(
                            output=empty_metric_file
                        ) )
                    ], name="picard_sort_sam." + sample.name, samples=[sample]

                    )

                )
                jobs.append(
                    concat_jobs([
                        picard.build_bam_index(
                        output_file,
                        output_file_index
                    )
                ], name = "build_bam_index." + sample.name, samples = [sample]
                ))
        else:
            log.info("skipping symlinks creation for duplicate marked bams....")

        return jobs

    def split_dragen_methylation_report(self):

        """
        Dragen methylation report contains all three methylation context.
To create combined CSV CpGs should be extracted from the dragen methylation report.
        """

        jobs = []


        for sample in self.samples:
            methylation_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            input_file = os.path.join(methylation_directory, sample.name + ".CX_report.txt")
            output_file = os.path.join(methylation_directory, sample.name + ".readset_sorted.dedup.CpG_report.txt.gz")

            jobs.append(
                concat_jobs([

                    dragen.split_dragen_methylation_report(
                        input_file,
                        output_file,
                        meth_contex="CG"
                    )
                ],
                    name="split_dragen_methylation_report." + sample.name,
                    samples=[sample]
                )
            )

        return jobs

    def dragen_bedgraph(self):

        """
        Creates bedgraph file from combined strand CpG file
        """
        jobs = []

        for sample in self.samples:
            methylation_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            input_file = os.path.join(methylation_directory, sample.name + ".readset_sorted.dedup.CpG_profile.strand.combined.csv")
            output_file = os.path.join(methylation_directory, sample.name + ".readset_sorted.dedup.bedGraph.gz")

            jobs.append(
                concat_jobs([

                    dragen.dragen_bedgraph(
                        input_file,
                        output_file
                    )
                ],
                    name="dragen_bedgraph." + sample.name,
                    samples=[sample]
                )
            )

        return jobs


    @property
    def steps(self):
        return [
            [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.bismark_align,
            self.add_bam_umi,               # step 5
            self.sambamba_merge_sam_files,
            self.picard_remove_duplicates,
            self.metrics,
            self.methylation_call,
            self.wiggle_tracks,             # step 10
            self.methylation_profile,
            self.ihec_sample_metrics_report,
            self.bis_snp,
            self.filter_snp_cpg,
            self.prepare_methylkit,         # step 15
            self.methylkit_differential_analysis,
            self.cram_output
        ], [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.dragen_align,
            self.add_bam_umi,               # step 5
            self.sambamba_merge_sam_files,
            self.picard_remove_duplicates,
            self.metrics,
            self.methylation_call,
            self.wiggle_tracks,  # step 10
            self.methylation_profile,
            self.ihec_sample_metrics_report,
            self.bis_snp,
            self.filter_snp_cpg,
            self.prepare_methylkit,  # step 15
            self.methylkit_differential_analysis,
            self.cram_output
        ],
            [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.dragen_align,
            self.add_bam_umi,               # step 5
            self.sambamba_merge_sam_files,
            self.sort_dragen_sam,
            self.metrics,
            self.dragen_methylation_call,
            self.split_dragen_methylation_report, # step 10
            self.methylation_profile,
            self.dragen_bedgraph,
            self.wiggle_tracks,
            self.ihec_sample_metrics_report,
            self.bis_snp, # step 15
            self.filter_snp_cpg,
            self.prepare_methylkit,
            self.methylkit_differential_analysis,
            self.cram_output
            ]
        ]

if __name__ == '__main__':
    argv = sys.argv

    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MethylSeq(protocol=['bismark', 'hybrid', 'dragen'])
