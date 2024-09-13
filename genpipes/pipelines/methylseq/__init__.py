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
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import itertools
import logging
import os
import re

# GenPipes Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs
from .. import dnaseq

from ...bfx import (
    bash_cmd as bash,
    bedtools,
    bismark,
    bissnp,
    bvatools,
    dragen,
    fgbio,
    gatk,
    gembs,
    htslib,
    igvtools,
    metrics,
    multiqc,
    picard2 as picard,
    sambamba,
    samtools,
    tools,
    ucsc
    )

log = logging.getLogger(__name__)

class MethylSeq(dnaseq.DnaSeqRaw):
    """
Methyl-Seq Pipeline
================

The GenPIpes Methyl-Seq pipeline now has four protocols.
 1. bismark
 2. gembs
 3. hybrid
 4. dragen

The "bismark" protocol uses Bismark to align reads to the reference genome. Picard is used to mark and remove duplicates and generate metric files. 

The "gembs" procotol uses GemBS for mapping and methylation and variant calling (http://statgen.cnag.cat/GEMBS/UserGuide/_build/html/index.html).

The "hybrid" protocl uses [Illumina Dragen Bio-IT processor](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html) and [dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm) software to align the reads to the reference genome. All the other steps are common with bismark protocol. The "dragen" protocol uses Dragen to align reads to the reference genome, call methylation, mark and remove duplicates.

Although dragen provides higher rate of mapping percentage with in a very short time duration (approximately three hours compared to 30 hours from bismark), it only accessible through McGill Genome Center cluster Abacus and The jobs cannot be submitted to any of the HPCs from the [Digital Research Aliance](https://status.computecanada.ca/). Importantly, the user needs to have permission to submit jobs to Abacus. Therefore, other users may continue to use only bismark protocol since it works in all the clusters.

However, if you would like to setup and use dragen in own cluster please refer to our [GenPipes Documentation](https://genpipes.readthedocs.io/en/latest/user_guide/pipelines/gp_wgs_methylseq.html)

A pipeline for processing and analyzing bisulfite sequencing data. The pipeline uses Bismark to align reads and extract methylation information, and Picard to remove duplicates, add read groups and index the BAM files. The pipeline also computes metrics and generates coverage tracks per sample. The pipeline currently supports the following protocols: bismark, hybrid and dragen.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.

Attributes:
    readsets (list): A list of readsets to process.
    output_dirs (dict): A dictionary of output directories.
    multiqc_inputs (list): A list of files to include in the MultiQC report.
Methods:
    bismark_align: Align reads with Bismark.
    add_bam_umi: Add read UMI tag to individual bam files using fgbio.
    picard_remove_duplicates: Remove duplicates using Picard.
    metrics: Compute metrics and generate coverage tracks per sample.
    methylation_call: Extract methylation information for individual cytosines.
    wiggle_tracks: Generate wiggle tracks suitable for multiple browsers.
    methylation_profile: Generate a CpG methylation profile by combining both forward and reverse strand Cs.
    ihec_sample_metrics_report: Retrieve the computed metrics which fit the IHEC standards and build a tsv report table for IHEC.
Parameters:
    protocol (str): Type of pipeline (default bismark).
    """

    def __init__(self, *args, protocol='bismark', **kwargs):
        self._protocol = protocol
        # Add pipeline specific arguments
        super(MethylSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-t", "--type", help="Type of pipeline (default chipseq)",
                                    choices=["bismark", "gembs", "hybrid", "dragen"], default="bismark", dest='protocol')
        return cls._argparser
    @property
    def readsets(self):
        """
        A list of readsets to process.
        Returns:
            list: A list of readsets.
        """
        # extra check on readsets
        for readset in super().readsets:
            if readset._run == "":
                _raise(SanitycheckError(f"Error: no run was provided for readset {readset.name}... Run has to be provided for all the readsets in order to use this pipeline."))
            if readset._lane == "":
                _raise(SanitycheckError(f"Error: no lane provided for readset {readset.name}... Lane has to be provided for all the readsets in order to use this pipeline."))

        return self._readsets

    @property
    def output_dirs(self):
        """
        Output directory paths.
        Returns:
            dict: Output directory paths.
        """
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

    @property
    def multiqc_inputs(self):
        """
        List of MultiQC input files.
        Returns:
            list: List of MultiQC input files.
        """
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    def bismark_align(self):
        """
        Align reads with [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf).
        Returns:
            list: A list of bismark jobs.
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, f"{readset.name}.trim.")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
            no_readgroup_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted_noRG.bam")
            output_bam = re.sub("_noRG.bam", ".bam", no_readgroup_bam)
            index_bam = re.sub("_noRG.bam", ".bam.bai", no_readgroup_bam)

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam), re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                # Defining the bismark output files (bismark sets the names of its output files from the basename of fastq1)
                # Note : these files will then be renamed (using a "mv" command) to fit with the GenPipes nomenclature (cf. no_readgroup_bam)
                bismark_out_bam = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_pe.bam", os.path.basename(fastq1)))
                bismark_out_report = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_PE_report.txt", os.path.basename(fastq1)))
                report_suffix = "_bismark_bt2_PE_report.txt"
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                # Defining the bismark output files (bismark sets the names of its output files from the basename of fastq1)
                # Note : these files will then be renamed (using a "mv" command) to fit with the GenPipes nomenclature (cf. no_readgroup_bam)
                bismark_out_bam = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2.bam", os.path.basename(fastq1)))
                bismark_out_report = os.path.join(alignment_directory, readset.name, re.sub(r'(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$', "_bismark_bt2_SE_report.txt", os.path.basename(fastq1)))
                report_suffix = "_bismark_bt2_SE_report.txt"
            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(output_bam)),
                    bash.mkdir(link_directory),
                    bismark.align(
                        fastq1,
                        fastq2,
                        os.path.dirname(no_readgroup_bam),
                        [no_readgroup_bam, re.sub(".bam", report_suffix, no_readgroup_bam)],
                    ),
                    Job(command=f"mv {bismark_out_bam} {no_readgroup_bam}"),
                    Job(command=f"mv {bismark_out_report} {re.sub('.bam', report_suffix, no_readgroup_bam)}"),
                    bash.ln(
                        os.path.relpath(re.sub(".bam", report_suffix, no_readgroup_bam), link_directory),
                        os.path.join(link_directory, readset.name + report_suffix),
                        input=re.sub(".bam", report_suffix, no_readgroup_bam)
                        )
                ],
                name=f"bismark_align.{readset.name}",
                samples=[readset.sample])
            )
            self.multiqc_inputs.append(re.sub(".bam", report_suffix, no_readgroup_bam))

            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    picard.add_read_groups(
                        no_readgroup_bam,
                        output_bam,
                        readset.name,
                        readset.library if readset.library else readset.sample.name,
                        f"{readset.run}_{readset.lane}",
                        readset.sample.name
                    )
                ],
                name=f"picard_add_read_groups.{readset.name}",
                samples=[readset.sample])
            )
            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    sambamba.index(
                        output_bam,
                        index_bam
                    )
                ],
                name=f"sambamba_index.{readset.name}",
                samples=[readset.sample])
            )
            jobs.append(
                concat_jobs([
                    bash.mkdir(alignment_directory),
                    sambamba.flagstat(
                        output_bam,
                        re.sub(".bam", "_flagstat.txt", output_bam),
                        global_conf.global_get('sambamba_flagstat', 'flagstat_options')
                    )
                ],
                name=f"sambamba_flagstat.{readset.name}",
                samples=[readset.sample])
            )

        return jobs

    def gembs_prepare(self):
        """
        Prepare metadata and config files for mapping with gemBS.
        """

        jobs = []

        metadata_file = os.path.join(self.output_dir, "metadata.csv")
        gembs_config_file = os.path.join(self.output_dir, "gembs.config")

        metadata_list = []
        trim_files = []
        
        for readset in self.readsets:
                
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
                
            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam), re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                trim_files.extend([fastq1, fastq2])
                
                metadata = ','.join([readset.sample.name,readset.name,readset.library,readset.sample.name,os.path.basename(fastq1),os.path.basename(fastq2)])
                
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                trim_files.append(fastq1)
    
                metadata = ','.join([readset.sample.name,readset.name,readset.library,readset.sample.name,fastq1])
    
            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))
                
            metadata_list.append(metadata)

        jobs.append(
                concat_jobs(
                    [
                        gembs.make_metadata(
                            metadata_list,
                            metadata_file
                            ),
                        gembs.make_config(
                            self.output_dir,
                            gembs_config_file
                            ),
                        gembs.prepare(
                            metadata_file,
                            gembs_config_file,
                            self.output_dir
                            )
                        ],
                    name="gembs_prepare",
                    input_dependency=[self.readsets_file.name] + trim_files
                    )
                )

        return jobs

    def gembs_map(self):
        """
        Map reads to reference genome with GemBS's gem-mapper.
        """

        jobs = []

        gembs_config = os.path.join(self.output_dir, ".gemBS/gemBS.mp")

        for sample in self.samples:
            alignment_dir = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            config_dir = os.path.join(alignment_dir, ".gemBS")
            trim_directory = os.path.join(self.output_dirs["trim_directory"], sample.name)
            trim_files = []

            for readset in sample.readsets:
                trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
                if readset.run_type == "PAIRED_END":
                    trim_fastqs = [trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]
                    trim_files.extend(trim_fastqs)
                elif readset.run_type == "SINGLE_END":
                    trim_fastq = trim_file_prefix + "single.fastq.gz"
                    trim_files.append(trim_fastq)
                else:
                    _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))


            jobs.append(
                    concat_jobs(
                        [
                            bash.rm(alignment_dir),
                            bash.mkdir(alignment_dir),
                            bash.mkdir(config_dir),
                            bash.cp(
                                gembs_config,
                                config_dir
                                ),
                            gembs.map(
                                sample.name,
                                alignment_dir
                                ),
                            bash.ln(
                                sample.name + ".bam",
                                os.path.join(alignment_dir, sample.name + ".sorted.bam"),
                                input = os.path.join(alignment_dir, sample.name + ".bam"))
                        ],
                        name = "gembs_map." + sample.name,
                        samples = [sample],
                        input_dependency=[gembs_config] + trim_files
                        )
                    )
        return jobs

    def add_bam_umi(self):
        """
        Add read UMI tag to individual bam files using [fgbio](https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html).
        Returns:
            list: A list of fgbio jobs.
        """

        jobs = []
        for readset in self.readsets:
            if readset.umi:
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
                input_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")
                output_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.UMI.bam")
                output_bai = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.UMI.bai")
                input_umi = readset.umi
                input_umi_corrected = os.path.join(self.output_dirs["corrected_umi_directory"], readset.name, f"{readset.name}.corrected.fastq.gz")

                #correct umi fastq name (removing space)

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.join(self.output_dirs["corrected_umi_directory"], readset.name)),
                            fgbio.correct_readname(
                                input_umi,
                                input_umi_corrected
                            ),
                    ],
                    name=f"fgbio_correct_readname.{readset.name}",
                    samples=[readset.sample]
                    )
                )

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_bam)),
                        fgbio.addumi(
                            input_bam,
                            input_umi_corrected,
                            output_bam,
                            output_bai
                        ),
                    ],
                    name=f"fgbio_addumi.{readset.name}",
                    samples=[readset.sample]
                    )
                )

        return jobs


    def picard_remove_duplicates(self):
        """
        Remove duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be removed as a duplicate in the BAM file. Removing duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        This step is used in bismark and hybrid protocols.
        Returns:
            list: A list of picard jobs.
        """

        jobs = []
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            alignment_file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.")
            input_file = f"{alignment_file_prefix}sorted.bam"
            bam_output = f"{alignment_file_prefix}sorted.dedup.bam"
            metrics_file = f"{alignment_file_prefix}sorted.dedup.metrics"
            flagstat_file = re.sub(".bam", "_flagstat.txt", bam_output)

            job = concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        picard.mark_duplicates(
                            [input_file],
                            bam_output,
                            metrics_file,
                            remove_duplicates="true",
                            ini_section='mark_duplicates'
                            ),
                        bash.ln(
                            os.path.relpath(metrics_file, link_directory),
                            os.path.join(link_directory, f"{sample.name}.sorted.dedup.metrics"),
                            input=metrics_file
                            )
                    ]
                )
            job.name = "mark_duplicates." + sample.name
            job.samples = [sample]
            self.multiqc_inputs.append(metrics_file)
            jobs.append(job)

            job = concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        samtools.flagstat(
                            bam_output,
                            flagstat_file
                            ),
                        bash.ln(
                            os.path.relpath(flagstat_file, link_directory),
                            os.path.join(link_directory, f"{sample.name}.sorted.dedup_flagstat.txt"),
                            input=flagstat_file
                            )
                    ]
                )
            job.name = f"samtools_flagstat_dedup.{sample.name}"
            job.samples = [sample]
            self.multiqc_inputs.append(flagstat_file)
            jobs.append(job)

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
        Returns:
            list: A list of metrics jobs.
        """

        # Check the library status
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

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        for sample in self.samples:
            file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.")
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            candidate_input_files = [[file_prefix + "bam"]]
            if bam[sample]:
                candidate_input_files.append([bam[sample]])
            [input_file] = self.select_input_files(candidate_input_files)

            # picard collect multiple metrics
            job = concat_jobs(
                    [
                        picard.collect_multiple_metrics(
                            input_file,
                            re.sub("bam", "all.metrics", input_file),
                            library_type=library[sample]
                        ),
                        bash.mkdir(link_directory)
                    ]
                )
            for outfile in job.report_files:
                self.multiqc_inputs.append(outfile)
                job = concat_jobs(
                        [
                            job,
                            bash.ln(
                                os.path.relpath(outfile, link_directory),
                                os.path.join(link_directory, os.path.basename(outfile)),
                                input=outfile
                            )
                        ]
                    )

            job.name = f"picard_collect_multiple_metrics.{sample.name}"
            job.samples = [sample]
            jobs.append(job)

            # Compute genome coverage with GATK
            job = gatk.depth_of_coverage(
                input_file,
                re.sub("bam", "all.coverage", input_file),
                coverage_bed
            )
            job.name = f"gatk_depth_of_coverage.genome.{sample.name}"
            job.samples = [sample]
            jobs.append(job)

            # Compute genome or target coverage with BVATools
            job = bvatools.depth_of_coverage(
                input_file,
                re.sub("bam", "coverage.tsv", input_file),
                coverage_bed,
                other_options=global_conf.global_get('bvatools_depth_of_coverage', 'other_options', required=False)
            )
            job.name = "bvatools_depth_of_coverage." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Get reads# raw and after duplication
            dedup_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.bam")
            raw_bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.bam")

            count_dedup_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.dedup.count")
            count_raw_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.raw.count")


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

                count_dedup_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.onTarget.dedup.count")
                count_raw_output = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.onTarget.raw.count")
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
                interval_list = re.sub(r"\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(
                        coverage_bed,
                        interval_list
                     )
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)
                file_prefix = os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.")
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
                    input_file,
                    re.sub(".bam", ".filtered_reads.counts.txt", input_file),
                    "-c " + global_conf.global_get('mapping_quality_filter', 'quality_threshold')
                )
            ])
            job.name = "mapping_quality_filter." + sample.name
            job.samples = [sample]
            jobs.append(job)

            # Calculate GC bias
            gc_content_file = re.sub(".bam", ".gc_cov.1M.txt", input_file)
            job = bedtools.coverage(
                input_file,
                gc_content_file
            )
            # For captured analysis
            if coverage_bed:
                gc_content_on_target_file = re.sub(".bam", ".gc_cov.1M.on_target.txt", input_file)
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
                    re.sub(".bam", ".GCBias_all.txt", input_file)
                )
            ])
            job.name = "GC_bias." + sample.name
            job.samples = [sample]
            jobs.append(job)

            job = igvtools.compute_tdf(input_file, input_file + ".tdf")
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
        Returns:
            list: A list of bismark methylation call jobs.
        """

        # Check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample]="SINGLE_END"
            if readset.run_type == "PAIRED_END" :
                library[readset.sample]="PAIRED_END"

        jobs = []

        methylseq_protocol = self.protocol

        methylation_protocol = global_conf.global_get('dragen_align', 'methylation_protocol', param_type='string', required=False)
        mapping_implementation = global_conf.global_get('dragen_align', 'mapping_implementation', param_type='string', required=False)
        link_directory = os.path.join(self.output_dirs["metrics_directory"],"multiqc_inputs")

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            input_file = os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")

            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            outputs = [
                os.path.join(methyl_directory, "CpG_context_" + re.sub( ".bam", ".txt.gz", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", ".bedGraph.gz", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", ".CpG_report.txt.gz", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", "_splitting_report.txt", os.path.basename(input_file))),
                os.path.join(methyl_directory, re.sub(".bam", ".M-bias.txt", os.path.basename(input_file)))
            ]

            if (methylseq_protocol == "hybrid" and methylation_protocol == "directional" and mapping_implementation=="single-pass"):
                jobs.append(
                    concat_jobs([
                        bash.mkdir(methyl_directory),
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

            else:
                jobs.append(
                    concat_jobs([
                        bash.mkdir(methyl_directory),
                        picard.sort_sam(
                            input_file,
                            re.sub("sorted", "readset_sorted", input_file),
                            "queryname"
                            )
                    ], name="picard_sort_sam." + sample.name, samples=[sample])
                )

            outputs = [re.sub("sorted", "readset_sorted", output) for output in outputs]
            bismark_job = concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        bismark.methyl_call(
                            re.sub("sorted", "readset_sorted", input_file),
                            outputs,
                            library[sample]
                            ),
                        bash.ln(
                            os.path.relpath(os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup_splitting_report.txt"), link_directory),
                            os.path.join(link_directory, f"{sample.name}.dedup_splitting_report.txt"),
                            input=os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup_splitting_report.txt")
                            ),
                        bash.ln(
                            os.path.relpath(os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.M-bias.txt"), link_directory),
                            os.path.join(link_directory, f"{sample.name}.dedup.M-bias.txt"),
                            input=os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.M-bias.txt")
                            )
                    ]
                )

            bismark_job.name = "bismark_methyl_call." + sample.name
            bismark_job.samples = [sample]
            self.multiqc_inputs.append(os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.M-bias.txt"))
            self.multiqc_inputs.append(os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup_splitting_report.txt"))
            jobs.append(bismark_job)

        return jobs

    def wiggle_tracks(self):
        """
        Generate wiggle tracks suitable for multiple browsers, to show coverage and methylation data.
        When using GRCh37 build of Human genome, to be compatible with the UCSC Genome Browser we only keep chromosomes 1-22, X, Y and MT,
        and we also rename them by prefixing "chr" to the chromosome anme (e.g. "1" becomes "chr1"), and changing the mitocondrial chromosome from "MT" to "chrM", and keeping the GRCh37 coordinates.
        Returns:
            list: A list of wiggle tracks jobs.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            # Generation of a bedGraph and a bigWig track to show the genome coverage
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")
            bed_graph_prefix = os.path.join(self.output_dirs["tracks_directory"], sample.name, sample.name)
            big_wig_prefix = os.path.join(self.output_dirs["tracks_directory"], "bigWig", sample.name)
            bed_graph_output = bed_graph_prefix + ".bedGraph"
            big_wig_output = big_wig_prefix + ".coverage.bw"

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.join(self.output_dirs["tracks_directory"], sample.name)),
                    bash.mkdir(os.path.join(self.output_dirs["tracks_directory"], "bigWig")),
                    bedtools.graph(
                        input_bam,
                        bed_graph_output,
                        ""
                    )
                ], name="bed_graph." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output)), samples=[sample])
            )
            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.join(self.output_dirs["tracks_directory"], "bigWig")),
                    ucsc.bedGraphToBigWig(
                        bed_graph_output,
                        big_wig_output,
                        False
                    )
                ], name="wiggle." + re.sub(".bw", "", os.path.basename(big_wig_output)), samples=[sample])
            )

            # Generation of a bigWig from the methylation bedGraph
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_bed_graph = os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.bedGraph.gz")
            output_wiggle = os.path.join(self.output_dirs["tracks_directory"], "bigWig", f"{sample.name}.methylation.bw")

            jobs.append(
                concat_jobs([
                    bash.mkdir(methyl_directory),
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
        Returns:
            list: A list of methylation profile jobs.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            cpg_input_file = os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.CpG_report.txt.gz")
            cpg_profile = re.sub(".CpG_report.txt.gz", ".CpG_profile.strand.combined.csv", cpg_input_file)

            # Generate CpG methylation profile
            job = tools.bismark_combine(
                cpg_input_file,
                cpg_profile
            )
            job.name = f"methylation_profile.{sample.name}"
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

                target_cpg_profile_count=os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.CpG_profile.strand.combined.on_target.count")
                job = metrics.target_cpg_profile(
                    target_cpg_profile,
                    target_cpg_profile_count,
                    sample.name
                )
                job.name = "metrics_target_CpG_profile." + sample.name
                job.samples = [sample]
                jobs.append(job)



            # Caluculate median & mean CpG coverage
            median_cpg_coverage = re.sub(".CpG_report.txt.gz", ".median_CpG_coverage.txt", cpg_input_file)
            job = tools.cpg_cov_stats(
                cpg_profile,
                median_cpg_coverage
            )
            job.name = f"median_CpG_coverage.{sample.name}"
            job.samples = [sample]
            if target_bed:
                job.removable_files = [target_cpg_profile]
            jobs.append(job)

        return jobs

    def ihec_sample_metrics_report(self):
        """
        Retrieve the computed metrics which fit the IHEC standards and build a tsv report table for IHEC.
        Note: The dragen protocol does not generate a metric for estimated library size. You will have to run Picard separately for this metric.
        Returns:
            list: A list of IHEC sample metrics report jobs.
        """

        jobs = []

        target_bed = bvatools.resolve_readset_coverage_bed(self.samples[0].readsets[0])
        metrics_all_file = os.path.join(self.output_dirs["metrics_directory"], "IHEC.sampleMetrics.stats")
        report_metrics_file = os.path.join(self.output_dirs["report_directory"], "IHEC.sampleMetricsTable.tsv")
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        ihec_multiqc_file = os.path.join(link_directory, "IHEC.sampleMetrics_mqc.tsv")

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
            metrics_file =  os.path.join(self.output_dirs["ihec_metrics_directory"], f"{sample.name}.read_stats.txt")

            metrics_output_list.append(metrics_file)
            # Add dependency for previous job in order to fill the allSample output correctly
            inputs.append(metrics_output_list[counter])

            # Trim log files
            for readset in sample.readsets:
                inputs.append((os.path.join(self.output_dirs["trim_directory"], sample.name, f"{readset.name}.trim.log")))

            # Aligned pre-deduplicated bam files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.bam"))

            # Deduplicated bam files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.bam"))

            # Coverage summary files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.all.coverage.sample_summary"))

            # Filtered reads count files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.filtered_reads.counts.txt"))

            # GC bias files
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.GCBias_all.txt"))

            if self.protocol == "bismark":
                # Bismark alignment files
                for readset in sample.readsets:
                    [bismark_report] = self.select_input_files([
                        [os.path.join(self.output_dirs["alignment_directory"], sample.name, readset.name, f"{readset.name}.sorted_noRG_bismark_bt2_PE_report.txt")],
                        [os.path.join(self.output_dirs["alignment_directory"], sample.name, readset.name, f"{readset.name}.sorted_noRG_bismark_bt2_SE_report.txt")]
                        ])
                    inputs.append(bismark_report)

            # CpG coverage files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.readset_sorted.dedup.median_CpG_coverage.txt"))

            # pUC19 methylation files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.readset_sorted.dedup.profile.pUC19.txt"))

            # Lambda conversion rate files
            inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.readset_sorted.dedup.profile.lambda.conversion.rate.tsv"))

            # CG stat files
            [cgstats_file] = self.select_input_files([
                [os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.sorted.dedup.profile.cgstats.txt")],
                [os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.readset_sorted.dedup.profile.cgstats.txt")]
            ])
            inputs.append(cgstats_file)

            #Estimated library sizes
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.dedup.metrics"))

            # read count (raw, dedup)
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.dedup.count"))
            inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.raw.count"))

            # read count (raw, dedup)  and CpG_profyle in targeted context
            if target_bed :
                inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.onTarget.dedup.count"))
                inputs.append(os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.onTarget.raw.count"))
                inputs.append(os.path.join(self.output_dirs["methylation_call_directory"], sample.name, f"{sample.name}.readset_sorted.dedup.CpG_profile.strand.combined.on_target.count"))


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

        report_multiqc_format_job = Job(
                [metrics_all_file],
                [ihec_multiqc_file],
                command="""\
echo -e "# plot_type: 'table'
# section_name: 'IHEC'
# description: 'Sequencing, Alignment and Methylation Metrics per Sample'
# headers:
#   raw_reads:
#       title: 'Raw Reads'
#       description: 'total number of reads obtained from the sequencer'
#       format: '{{:,.0f}}'
#       placement: 900
#   trimmed_reads:
#       title: 'Trimmed Reads'
#       description: 'number of remaining reads after the trimming step'
#       format: '{{:,.0f}}'
#       placement: 910
#   perc_survival_rate:
#       title: '% Survival Rate'
#       description: 'trimmed_reads / raw_reads * 100'
#       placement: 920
#   aligned_reads:
#       title: 'Aligned Reads'
#       description: 'number of aligned reads to the reference'
#       format: '{{:,.0f}}'
#       placement: 930
#   perc_mapping_efficiency:
#       title: '% Mapping Efficiency'
#       description: 'aligned_reads / trimmed_reads * 100'
#       placement: 940
#   duplicated_reads:
#       title: 'Duplicated Reads'
#       description: 'number of duplicated read entries providing the same mapping coordinates (due to PCR duplicates)'
#       format: '{{:,.0f}}'
#       placement: 950
#   perc_duplication_rate:
#       title: '% Duplication Rate'
#       description: 'duplicated_reads / aligned_reads * 100'
#       placement: 960
#   deduplicated_aligned_reads:
#       title: 'Deduplicated Aligned Reads'
#       description: 'aligned_reads - duplicated_reads'
#       format: '{{:,.0f}}'
#       placement: 970
#   perc_useful_aligned_rate:
#       title: '% Useful Aligned Rate'
#       description: 'deduplicated_aligned_reads / raw_reads * 100'
#       placement: 980
#   perc_proportion_unique_filtered_reads_MAPQ>10:
#       title: '% Unique Filtered Reads MAPQ>10'
#       description: 'deduplicated_aligned_reads with mapping quality > 10 / trimmed_reads * 100'
#       placement: 990
#   GC_bias:
#       title: 'GC bias'
#       description: 'the Pearson correlation between coverage values and GC content for 1000 bins of 100 base pair across genome'
#       placement: 1000
#   perc_pUC19_methylation_rate:
#       title: '% pUC19 methylation rate'
#       description: '100 - C->T conversion rate on pUC19 * 100'
#       placement: 1010
#   perc_lambda_conversion_rate:
#       title: '% Lambda Conversion Rate'
#       description: 'C->T conversion rate on the lambda phage * 100 '
#       placement: 1020
#   perc_human_conversion:
#       title: '% Human Conversion Rate'
#       description: 'estimation of non-CpG methylation * 100'
#       placement: 1030
#   estimated_average_genome_coverage:
#       title: 'Est. Avg. Genome Coverage'
#       description: 'aligned_reads / genome size'
#       placement: 1040
#   median_CpG_coverage:
#       title: 'Median CpG Coverage'
#       description: 'median read coverage for on-target CpGs'
#       placement: 1050
#   num_CpG_1X:
#       title: 'CpGs at 1X'
#       description: 'total number of CpGs with a coverage >= 1x'
#       format: '{{:,.0f}}'
#       placement: 1060
#   num_CpG_10X:
#       title: 'CpGs at 10X'
#       description: 'total number of CpGs with a coverage >= 10x'
#       format: '{{:,.0f}}'
#       placement: 1070
#   num_CpG_30X:
#       title: 'CpGs at 30X'
#       description: 'total number of CpGs with a coverage >= 30x'
#       placement: 1080
#   estimate_library_size:
#       title: 'Est. Library Size'
#       description: 'estimate library size'
#       format: '{{:,.0f}}'
#       placement: 1090" > {ihec_multiqc_file}

cat {metrics_all_file} | sed 's/%_/perc_/g' | sed 's/#_/num_/g' >> {ihec_multiqc_file}""".format(
                ihec_multiqc_file=ihec_multiqc_file,
                metrics_all_file=metrics_all_file
                )
            )

        jobs.append(
            concat_jobs([
                bash.mkdir(self.output_dirs["metrics_directory"]),
                bash.mkdir(link_directory),
                report_multiqc_format_job
                ],
            name="ihec_sample_metrics_report"
            )
        )

        self.multiqc_inputs.append(ihec_multiqc_file)

        return jobs

    def bis_snp(self):
        """
        SNP calling with [BisSNP](https://people.csail.mit.edu/dnaase/bissnp2011/).
        Returns:
            list: A list of BisSNP jobs.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            input_file = os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")

            variant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name)
            cpg_output_file = os.path.join(variant_directory, f"{sample.name}.cpg.vcf")
            snp_output_file = os.path.join(variant_directory, f"{sample.name}.snp.vcf")

            jobs.append(
                concat_jobs([
                    bash.mkdir(variant_directory),
                    bissnp.bisulfite_genotyper(
                        input_file,
                        cpg_output_file,
                        snp_output_file
                        ),
                    htslib.bgzip(
                        cpg_output_file,
                        f"{cpg_output_file}.gz"
                        ),
                    htslib.bgzip(
                        snp_output_file,
                        f"{snp_output_file}.gz"
                        )
                ],
                name=f"bissnp.{sample.name}",
                samples=[sample], 
                output_dependency=[f"{snp_output_file}.gz", f"{cpg_output_file}.gz"],
                removable_files=[cpg_output_file, snp_output_file]
                )
            )

        return jobs

    def gembs_call(self):
        """
        Methylation calling with bs_call as part of GemBS pipeline
        """
        jobs = []
        
        for sample in self.samples:
            bam = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".bam")
            output_dir = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            output_prefix = os.path.join(output_dir, sample.name)
            config_dir = os.path.join(output_dir, ".gemBS")
            gembs_config = os.path.join(self.output_dir, ".gemBS/gemBS.mp")
            
            jobs.append(
                    concat_jobs(
                        [
                            bash.rm(output_dir),
                            bash.mkdir(output_dir),
                            bash.mkdir(config_dir),
                            bash.cp(
                                gembs_config,
                                config_dir
                                ),
                            gembs.call(
                                sample.name,
                                bam,
                                output_prefix
                                )
                            ],
                        name = "gembs_call." + sample.name,
                        samples = [sample],
                        input_dependency=[bam,gembs_config]
                        )
                    )
            
            variants_dir = os.path.join(self.output_dirs["variants_directory"], sample.name)
            config_dir = os.path.join(variants_dir, ".gemBS")
            input = output_prefix + ".bcf"

            jobs.append(
                    concat_jobs(
                        [
                            bash.rm(variants_dir),
                            bash.mkdir(variants_dir),
                            bash.mkdir(config_dir),
                            bash.cp(
                                gembs_config,
                                config_dir
                                ),
                            gembs.extract(
                                input,
                                sample.name,
                                variants_dir
                                )
                            ],
                        name = "gembs_extract." + sample.name,
                        samples = [sample],
                        input_dependency = [input]
                        )
                    )
        return jobs

    def gembs_bcf_to_vcf(self):
        """
        Create vcf of SNPs with bedtools intersect, by intersecting gemBS .bcf with SNP DB.
        """
        jobs = []

        # create vcf of snps instead of using gemBS snp output format

        for sample in self.samples:
            bcf_dir = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            bcf_prefix = os.path.join(bcf_dir, sample.name)
            variants_dir = os.path.join(self.output_dirs["variants_directory"], sample.name)
            snp_output = os.path.join(variants_dir, sample.name + "_snps.vcf")
            input = bcf_prefix + ".bcf"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            variants_dir
                            ),
                        tools.gembs_bcf_to_vcf(
                            input,
                            snp_output
                            ),
                        bash.gzip(
                            snp_output,
                            None,
                            '-f '
                            )
                    ],
                    name = "gembs_bcf_to_vcf." + sample.name,
                    samples = [sample],
                    input_dependency=[input],
                    output_dependency=[snp_output + ".gz"]
                        )
                    )

        return jobs

    def gembs_report(self):
        """
        GemBS report
        """
        jobs = []

        report_dir = self.output_dirs["report_directory"]
        project = global_conf.global_get("gembs_report", "project_name")
        report = os.path.join(report_dir, project + "_QC_Report.html")
        inputs = []

        for sample in self.samples:
            map_json = os.path.join(self.output_dirs["alignment_directory"], sample.name, sample.name + ".json")
            call_json = os.path.join(self.output_dirs["methylation_call_directory"], sample.name, sample.name + ".json")
            inputs.extend([map_json, call_json])

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(report_dir),
                        gembs.report(
                            inputs,
                            report
                            )
                    ],
                    name = "gembs_report"
                )
            )
        
        return jobs

    def gembs_format_cpg_report(self):
        """
        Reformat gemBS output to match bismark and dragen output, so following steps can be followed. 
        """
        
        jobs = []

        for sample in self.samples:
            cpg_report_input = os.path.join(self.output_dirs["variants_directory"], sample.name, sample.name + "_cpg.bed.gz")
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            cpg_output = os.path.join(methyl_directory, sample.name + ".readset_sorted.dedup.CpG_report.txt.gz")

            job = tools.gembs_format_cpg_report(
                        cpg_report_input,
                        cpg_output
                        )
            job.name = "gembs_format_cpg_report." + sample.name
            job.samples = [sample]

            jobs.append(job)

        return jobs

    def filter_snp_cpg(self):
        """
        SNP CpGs filtering.
        Returns:
            list: A list of SNP CpGs filtering jobs.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_file = os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.CpG_profile.strand.combined.csv")

            output_file = re.sub("CpG_profile.strand.combined.csv", "filtered.bed", input_file)

            jobs.append(
                concat_jobs([
                    bash.mkdir(methyl_directory),
                    tools.filter_snp_cpg(
                        input_file,
                        output_file
                    )
                ],
                name=f"filter_snp_cpg.{sample.name}",
                samples=[sample]
                )
            )
        return jobs

    def prepare_methylkit(self):
        """
        Prepare input file for [methylKit](https://www.bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html) differential analysis.
        Returns:
            list: A list of methylKit preparation jobs.
        """

        jobs = []
        for sample in self.samples:
            methyl_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            input_file = os.path.join(methyl_directory, f"{sample.name}.readset_sorted.dedup.filtered.bed")

            output_directory = os.path.join(self.output_dirs["methylkit_directory"], "inputs")
            output_file = os.path.join(output_directory, re.sub("filtered.bed", "map.input", os.path.basename(input_file)))

            jobs.append(
                concat_jobs([
                    bash.mkdir(output_directory),
                    tools.prepare_methylkit(
                        input_file,
                        output_file
                    )
                ],
                name=f"prepare_methylkit.{sample.name}",
                samples=[sample]
                )
            )
        return jobs

    def methylkit_differential_analysis(self):
        """
        Run methylKit to get DMCs & DMRs for different design comparisons.
        Returns:
            list: A list of methylKit differential analysis jobs.
        """

        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)

        input_directory = os.path.join(self.output_dirs["methylkit_directory"], "inputs")
        input_files = []
        for contrast in self.contrasts:
            input_files.extend([[os.path.join(input_directory, f"{sample.name}.readset_sorted.dedup.map.input") for sample in group] for group in [contrast.controls, contrast.treatments]])

        input_files = list(itertools.chain.from_iterable(input_files))

        output_directory = os.path.join(self.output_dirs["methylkit_directory"], "results")
        output_files = [os.path.join(output_directory, "Rdata_files", contrast.name, "perbase.testingresults.txt.gz") for contrast in self.contrasts]

        methylkit_job = tools.methylkit_differential_analysis(
            design_file,
            input_files,
            output_files,
            output_directory
        )

        return [
            concat_jobs([
                bash.mkdir(output_directory),
                methylkit_job
                ],
                name="methylkit_differential_analysis"
                )]

    def dragen_align(self):
        """
        Align reads with [Dragen](https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/FrontPages/DRAGEN.htm) both hybrid and dragen protocols use this step to align reads.
        The Dragen parameters can be changed using other_options of the ini configuration.
        Returns:
            list: A list of Dragen alignment jobs.
        """

        jobs = []
        methylseq_protocol = self.protocol

        methylation_protocol = global_conf.global_get('dragen_align', 'methylation_protocol', param_type='string')
        mapping_implementation = global_conf.global_get('dragen_align', 'mapping_implementation', param_type='string')
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
                trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, f"{readset.name}.trim")
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
                # All the input files should copy to dragen work folder because IO operations with default locations
                # are really slow or have permission issues. After processing files all the files should move back to
                # the default GenPipes folder and remove files from dragen work folder
                dragen_inputfolder = os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "reads", readset.name)
                dragen_workfolder = os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "alignment", readset.name)
                # dragen output file name
                dragen_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")
                dragen_report_files = [
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.mapping_metrics.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.wgs_coverage_metrics.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.wgs_fine_hist.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.wgs_contig_mean_cov.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.fragment_length_hist.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.trimmer_metrics.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.time_metrics.csv"),
                        os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.fastqc_metrics.csv")
                        ]
                index_bam = f"{dragen_bam}.bai"

                # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
                if readset.run_type == "PAIRED_END":
                    candidate_input_files = [[f"{trim_file_prefix}.pair1.fastq.gz", f"{trim_file_prefix}.pair2.fastq.gz"]]
                    if readset.fastq1 and readset.fastq2:
                        candidate_input_files.append([readset.fastq1, readset.fastq2])
                    if readset.bam:
                        candidate_input_files.append(
                            [
                                re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam),
                                re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)
                            ]
                        )
                    [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                    # first fastqs need to be copy to the dragen work_folder
                    dragen_tmp_fastq1 = os.path.join(dragen_inputfolder, f"{readset.name}.pair1.fastq.gz")
                    dragen_tmp_fastq2 = os.path.join(dragen_inputfolder, f"{readset.name}.pair2.fastq.gz")
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
                    candidate_input_files = [[f"{trim_file_prefix}.single.fastq.gz"]]
                    if readset.fastq1:
                        candidate_input_files.append([readset.fastq1])
                    if readset.bam:
                        candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                    [fastq1] = self.select_input_files(candidate_input_files)
                    fastq2 = None
                    # first fastqs need to be copy to the dragen work_folder
                    dragen_tmp_fastq1 = os.path.join(dragen_workfolder, f"{readset.name}.single.fastq.gz")
                    dragen_tmp_fastq2 = fastq2
                    cp_dragen_fastq_job = concat_jobs([
                        bash.mkdir(dragen_inputfolder ),
                        bash.cp(os.path.abspath(fastq1), dragen_tmp_fastq1)
                    ],
                    name=f"dragen_copy_fastq.{readset.name}",
                    samples=[readset.sample]
                    )
                    rm_dragen_fastq_job = concat_jobs([
                        bash.rm(dragen_tmp_fastq1, recursive=True, force=True)
                    ],
                    name=f"dragen_remove_fastq.{readset.name}",
                    samples=[readset.sample]
                    )
                else:
                    _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))
                dragen_align = concat_jobs([
                    # To track the dragen output files and done file input dependency and output dependency of the original location should be defined
                    cp_dragen_fastq_job,
                    bash.mkdir(dragen_workfolder),
                    bash.mkdir(os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "job_output", "dragen_align")),
                    bash.mkdir(os.path.abspath(alignment_directory)),
                    dragen.align_methylation(
                        dragen_tmp_fastq1,
                        dragen_tmp_fastq2,
                        dragen_workfolder,
                        readset.name,
                        readset.sample.name,
                        readset.library if readset.library else readset.sample.name,
                        f"{readset.name}_{readset.run}_{readset.lane}",
                        protocol=methylseq_protocol

                    ),
                    bash.cp(dragen_workfolder, f"{os.path.abspath(alignment_directory)}/", recursive=True),
                    bash.rm(dragen_workfolder, recursive=True, force=True),
                    rm_dragen_fastq_job
                ],
                name=f"dragen_align.{readset.name}",
                samples=[readset.sample],
                input_dependency=[fastq1, fastq2],
                output_dependency=[dragen_bam] + dragen_report_files)
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
                    ],
                    name=f"samtools_flagstat.{readset.name}",
                    samples=[readset.sample],
                    output_dependency = [re.sub(".bam", "_flagstat.txt", dragen_bam)]
                    )

                )
                # Create symlinks to dragen_align metrics outputs to be used for multiqc
                link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
                symlink_job = bash.mkdir(link_directory)
                for outfile in dragen_align.report_files:
                    self.multiqc_inputs.append(os.path.join(alignment_directory, readset.name, outfile))
                    symlink_job = concat_jobs(
                            [
                                symlink_job,
                                bash.ln(
                                    os.path.relpath(os.path.join(alignment_directory, readset.name, outfile), link_directory),
                                    os.path.join(link_directory, outfile),
                                    input=os.path.join(alignment_directory, readset.name, outfile)
                                    )
                                ]
                            )
                symlink_job.name = f"symlink_dragen_metrics.{readset.name}"
                symlink_job.samples = [readset.sample]
                jobs.append(symlink_job)

                jobs.append(
                    concat_jobs([
                        bash.mkdir(alignment_directory),
                        sambamba.index(
                            dragen_bam,
                            index_bam
                        )
                    ],
                    name=f"sambamba_index.{readset.name}",
                    samples=[readset.sample])
                    )

            return jobs

        else:
            _raise(SanitycheckError("Please use \"dragen\" protocol when using directional-complement on single-pass mode. Skipping generating the genpipes file..."))

    def dragen_methylation_call(self):
        """
        Call methylation with Dragen using the 2nd run of Dragen alignment.
        Returns:
            list: A list of Dragen methylation call jobs.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            methylation_call_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)
            dragen_inputfolder = os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "reads", sample.name)
            dragen_workfolder = os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "dragen_methylation_call", sample.name)
            dragen_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            output_report = os.path.join(methylation_call_directory, f"{sample.name}.CX_report.txt")

            dragen_tmp_bam = os.path.join(dragen_inputfolder, f"{sample.name}.sorted.bam")
            cp_dragen_bam_job = concat_jobs([
                bash.mkdir(dragen_inputfolder),
                bash.cp(os.path.abspath(dragen_bam), dragen_tmp_bam)
                ],
                name=f"dragen_copy_bam.{sample.name}",
                samples=[sample]
                )
            rm_dragen_bam_job = concat_jobs([
                bash.rm(dragen_tmp_bam, recursive=True, force=True)
                ],
                name=f"dragen_remove_bam.{sample.name}",
                samples=[sample]
                )

            dragen_methylationn_call = concat_jobs([
                cp_dragen_bam_job,
                bash.mkdir(dragen_workfolder),
                bash.mkdir(os.path.join(global_conf.global_get('dragen_align', 'work_folder'), "job_output", "dragen_methylation_call")),
                bash.mkdir(os.path.abspath(methylation_call_directory)),
                dragen.call_methylation(
                    dragen_tmp_bam,
                    dragen_workfolder,
                    sample.name, output=output_report),
                bash.cp(dragen_workfolder, os.path.abspath(self.output_dirs["methylation_call_directory"]) + "/", recursive=True),
                bash.rm(dragen_workfolder, recursive=True, force=True ),
                rm_dragen_bam_job
                ],
                name=f"dragen_methylation_call.{sample.name}",
                samples=[sample],
                input_dependency=[dragen_bam],
                output_dependency=[output_report]
                )
            jobs.append(
                dragen_methylationn_call
            )
        return jobs

    def sort_dragen_sam(self):
        """
        Coordinate sorting the bam file resulted from dragen and create an index.
        Returns:
            list: A list of Dragen sam sorting jobs.
        """
        jobs = []

        duplicate_marking = global_conf.global_get('dragen_align', 'duplicate_marking', param_type='boolean')

        if duplicate_marking is True:

            for sample in self.samples:
                alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

                input_file = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
                output_file = os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")
                output_file_index = os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam.bai")
                # empty metric file will be created for "ihec_sample_metrics_report" steps. Otherwise it will be failed.
                # Therefore added an empty file with ESTIMATED_LIBRARY_SIZE = NA
                # so the dragen protocol will not estimate the library size
                empty_metric_file= os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.metrics")

                jobs.append(
                    concat_jobs([
                        bash.mkdir(alignment_directory),
                        picard.sort_sam(
                            input_file,
                            output_file
                            ),
                        Job(
                            output_files=[empty_metric_file],
                            command=f"""printf "ESTIMATED_LIBRARY_SIZE\\nNA" > {empty_metric_file}"""
                            )
                        ],
                        name=f"picard_sort_sam.{sample.name}",
                        samples=[sample]
                        )
                    )
                jobs.append(
                    concat_jobs([
                        picard.build_bam_index(
                            output_file,
                            output_file_index
                            )
                        ],
                        name=f"build_bam_index.{sample.name}",
                        samples=[sample]
                        )
                    )
        else:
            log.info("skipping symlinks creation for duplicate marked bams....")

        return jobs

    def split_dragen_methylation_report(self):

        """
        Dragen methylation report contains all three methylation context.
        To create combined CSV CpGs should be extracted from the dragen methylation report.
        Returns:
            list: A list of split dragen methylation report
        """

        jobs = []

        for sample in self.samples:
            methylation_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            input_file = os.path.join(methylation_directory, f"{sample.name}.CX_report.txt")
            output_file = os.path.join(methylation_directory, f"{sample.name}.readset_sorted.dedup.CpG_report.txt.gz")

            jobs.append(
                concat_jobs([
                    dragen.split_dragen_methylation_report(
                        input_file,
                        output_file,
                        meth_contex="CG"
                    )
                ],
                    name=f"split_dragen_methylation_report.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def dragen_bedgraph(self):

        """
        Creates bedgraph file from combined strand CpG file
        Returns:
            list: A list of dragen bedgraph jobs
        """
        jobs = []

        for sample in self.samples:
            methylation_directory = os.path.join(self.output_dirs["methylation_call_directory"], sample.name)

            input_file = os.path.join(methylation_directory, f"{sample.name}.readset_sorted.dedup.CpG_profile.strand.combined.csv")
            output_file = os.path.join(methylation_directory, f"{sample.name}.readset_sorted.dedup.bedGraph.gz")

            jobs.append(
                concat_jobs([

                    dragen.dragen_bedgraph(
                        input_file,
                        output_file
                    )
                ],
                    name=f"dragen_bedgraph.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report.
        MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
        perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).
        Returns:
            list: A list of MultiQC jobs.
        """
        jobs = []

        input_links = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
        output = os.path.join(self.output_dirs['report_directory'], f"MethylSeq.{self.protocol}.multiqc")

        job = concat_jobs(
            [
                bash.mkdir(os.path.join(self.output_dirs['report_directory'])),
                multiqc.run(
                    [input_links],
                    output
                )
            ]
        )
        job.name = "multiqc"
        job.input_files = self.multiqc_inputs
        jobs.append(job)

        return jobs


    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        """
        Returns the protocol for the pipeline.
        Returns:
            dict: A dictionary of protocols for the pipeline.
        """
        return {'bismark':
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.bismark_align,
                self.add_bam_umi,              # step 5
                self.sambamba_merge_sam_files,
                self.picard_remove_duplicates,
                self.metrics,
                self.methylation_call,
                self.wiggle_tracks,           # step 10
                self.methylation_profile,
                self.ihec_sample_metrics_report,
                self.bis_snp,
                self.filter_snp_cpg,
                self.prepare_methylkit,           # step 15
                self.methylkit_differential_analysis,
                self.multiqc,
                self.cram_output
            ], 'gembs':
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.gembs_prepare,
                self.gembs_map,   
            #    self.add_bam_umi, # important to include? Slightly annoying because merge step not needed for gemBS
                self.picard_remove_duplicates,
                self.metrics,
                self.gembs_call,
                self.gembs_bcf_to_vcf,
                self.gembs_format_cpg_report,
                self.methylation_profile,
                self.dragen_bedgraph,
                self.wiggle_tracks,
                self.ihec_sample_metrics_report,
                self.gembs_report,
                self.filter_snp_cpg,
                self.prepare_methylkit,
                self.methylkit_differential_analysis,
                self.multiqc,
                self.cram_output
            ], 'hybrid': 
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.dragen_align,
                self.add_bam_umi,              # step 5
                self.sambamba_merge_sam_files,
                self.picard_remove_duplicates,
                self.metrics,
                self.methylation_call,
                self.wiggle_tracks,             # step 10
                self.methylation_profile,
                self.ihec_sample_metrics_report,
                self.bis_snp,
                self.filter_snp_cpg,
                self.prepare_methylkit,            # step 15
                self.methylkit_differential_analysis,
                self.multiqc,
                self.cram_output
            ], 'dragen':
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
                self.multiqc,
                self.cram_output
            ]
        }

def main(parsed_args):
    """
    The function that will call this pipeline!
    Parameters:
        parsed_args (argparse.Namespace): The parsed arguments from the command line.
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    no_json = parsed_args.no_json
    json_pt = parsed_args.json_pt
    force = parsed_args.force
    force_mem_per_cpu = parsed_args.force_mem_per_cpu
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file
    protocol = parsed_args.protocol

    pipeline = MethylSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container, protocol=protocol)

    pipeline.submit_jobs()
