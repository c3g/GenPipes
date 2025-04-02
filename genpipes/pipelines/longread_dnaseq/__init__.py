################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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
import os
import logging
import re

# GenPipes Modules
from ...core.config import global_conf, SanitycheckError, _raise
from ...core.job import Job, concat_jobs, pipe_jobs
from .. import common

from ...bfx import (
    annotsv,
    bash_cmd as bash,
    bcftools,
    bvatools,
    cpsr,
    deepvariant,
    gatk4,
    hificnv,
    hiphase,
    htslib,
    job2json_project_tracking,
    minimap2,
    mosdepth,
    multiqc,
    nanoplot,
    pbmm2,
    pycoqc,
    sambamba,
    samtools,
    sawfish,
    svim,
    tools,
    trgt
    )

log = logging.getLogger(__name__)

class LongReadDnaSeq(common.LongRead):
    """
LongRead DNA-Seq Pipeline
==============

The LongRead Pipeline is used to analyse long reads produced by the Oxford Nanopore Technologies (ONT) 
and PacBio Revio sequencers. The protocols used are nanopore and revio, respectively.

Currently, the nanopore protocol of the pipeline uses minimap2 to align reads to the reference genome.
Additionally, it produces a QC report that includes an interactive dashboard with data from the basecalling
summary file as well as the alignment. A step aligning random reads to the NCBI nt database and reporting 
the species of the highest hits is also done as QC.

Once the QC and alignments have been produced, Picard is used to merge readsets coming from the same
sample. Finally, SVIM is used to detect Structural Variants (SV) including deletions, insertions and
translocations. For a full summary of the types of SVs detected, please consult the following [site](
https://github.com/eldariont/svim#background-on-structural-variants-and-long-reads).

The SV calls produced by SVIM are saved as VCFs for each sample, which can then be used in downstream
analyses. No filtering is performed on the SV calls.

This pipeline currently does not perform base calling and requires both FASTQ and a sequencing_summary
file produced by a ONT supported basecaller (we recommend Guppy). Additionally, the testing and
development of the pipeline were focused on genomics applications, and functionality has not been tested
for transcriptomics or epigenomics datasets.

For more information on using ONT data for structural variant detection, as well as an alternative
approach, please consult [this GitHub repository](https://github.com/nanoporetech/pipeline-structural-variation).

The Revio protocol uses pbmm2 to align reads to the reference genome, followed by variant calling with DeepVariant
and structural variant calling with HiFiCNV, TRGT, and Sawfish. Variants are annotated with AnnotSV and phased
with HiPhase. A CPSR report can be produced from the phased variants. Metrics on the raw and mapped reads are
collected with NanoPlot and mosdepth, respectively. 

Both protocols require as input a readset file, which provides sample metadata and paths to input data (FASTQ, FAST5 or BAM).
For information on the structure and contents of the LongRead readset file, please consult [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/readset_file.html).
    """

    def __init__(self, *args, protocol='nanopore', **kwargs):
        self.protocol = protocol
        super(LongReadDnaSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument(
            "-t",
            "--type",
            help="Type of pipeline (default nanopore)",
            dest='protocol',
            choices=["nanopore", "revio"],
            default="nanopore"
            )
        return cls._argparser

    @property
    def output_dirs(self):
        dirs = {
            'blastqc_directory': os.path.relpath(os.path.join(self.output_dir, 'blastQC'), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'pycoqc_directory': os.path.relpath(os.path.join(self.output_dir, 'pycoQC'), self.output_dir),
            'svim_directory': os.path.relpath(os.path.join(self.output_dir, 'svim'), self.output_dir),
            'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
            'SVariants_directory': os.path.relpath(os.path.join(self.output_dir, 'SVariants'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
            'annotsv_directory': os.path.relpath(os.path.join(self.output_dir, 'annotSV'), self.output_dir),
            'hiphase_directory': os.path.relpath(os.path.join(self.output_dir, 'hiphase'), self.output_dir),
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
    
    def guppy(self):
        """
        Use the Guppy basecaller to perform basecalling on all raw fast5 files.
        Uses the 'flip-flop' basecalling model by default.
        IN DEVELOPMENT
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

            blast_directory = os.path.join(self.output_dirs["blastqc_directory"], readset.name)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            job = tools.sh_blastQC_ONT(
                blast_directory,
                reads_fastq_dir,
                readset.name
            )
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs
    
    def metrics_nanoplot(self):
        """
        Collect QC metrics on unaligned bam or fastq files with nanoplot.
        """
        jobs =[]

        for readset in self.readsets:
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'], readset.sample.name)
            nanoplot_directory = os.path.join(metrics_directory, "nanoplot")
            nanoplot_prefix = f"{readset.name}."

            if readset.fastq_files:
                input_fastq = os.path.join(nanoplot_directory, os.path.basename(readset.fastq_files))
                link_job = bash.ln(
                    os.path.abspath(readset.fastq_files),
                    input_fastq,
                    readset.fastq_files
                    )
                input_bam = None
            elif readset.bam:
                input_bam = os.path.join(nanoplot_directory, os.path.basename(readset.bam))
                link_job = bash.ln(
                    os.path.abspath(readset.bam),
                    input_bam,
                    readset.bam
                    )
                input_fastq = None
            else:
                _raise(SanitycheckError(f"Error: Neither BAM nor FASTQ file available for readset {readset.name} !"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(nanoplot_directory),
                        link_job,
                        nanoplot.qc(
                            nanoplot_directory,
                            nanoplot_prefix,
                            input_bam,
                            input_fastq
                        )
                    ],
                    name=f"nanoplot.{readset.name}",
                    samples=[readset.sample],
                    readsets=[readset]
                )
            )

            self.multiqc_inputs.append(
                os.path.join(nanoplot_directory, f"{nanoplot_prefix}NanoStats.txt")
                )

        return jobs

    def minimap2_align(self):
        """
        Uses minimap2 to align the Fastq reads that passed the minimum QC threshold to
        the provided reference genome. By default, it aligns to GRCh38.
        """
        jobs = []

        for readset in self.readsets:

            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name)
            out_bam = os.path.join(alignment_directory, readset.name + ".sorted.bam")
            out_bai = re.sub(r"\.bam$", ".bam.bai", out_bam)

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
            job = concat_jobs(
                [
                    pipe_jobs(
                        [
                            bash.mkdir(os.path.dirname(out_bam)),
                            minimap2.minimap2_ont(
                                reads_fastq_dir,
                                read_group,
                                ini_section= "minimap2_align"
                            ),
                            sambamba.view(
                                "/dev/stdin",
                                None,
                                options="-S -f bam"
                            ),
                            sambamba.sort(
                                "/dev/stdin",
                                out_bam,
                                tmp_dir=global_conf.global_get('minimap2_align', 'tmp_dir', required=True),
                            )
                        ]
                    ),
                    samtools.quickcheck(
                        out_bam
                    ),
                    sambamba.index(
                        out_bam,
                        out_bai,
                    )
                ],
                name="minimap2_align." + readset.name,
                samples=[readset.sample],
                input_dependency=[reads_fastq_dir]
            )
            jobs.append(job)

        return jobs

    def pbmm2_align(self):
        """
        Uses pbmm2 to align fastq files or the raw hifi bam to the reference.
        """
        jobs = []

        for readset in self.readsets:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name)
            output_bam = os.path.join(alignment_directory, readset.name + ".sorted.bam")
            
            if readset.bam:
                read_group = None
                input_file = readset.bam
            elif readset.fastq_files:
                read_group = "'@RG" + \
                "\\tID:" + readset.name + \
                "\\tSM:" + readset.sample.name + \
                "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                "\\tPU:1\\tPL:pacbio_revio'"
                input_file = readset.fastq_files
            else:
                _raise(SanitycheckError(f"Error: Neither BAM nor FASTQ file available for readset {readset.name} !"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        pbmm2.align(
                            input_file,
                            read_group,
                            readset.sample.name,
                            output_bam,
                            sort=True
                        )
                    ],
                    name="pbmm2_align." + readset.name,
                    samples=[readset.sample],
                    readsets=[readset]
                )
            )

        return jobs
    
    def pycoqc(self):
        """
        Use pycoQC to produce an interactive quality report based on the summary file and
        alignment outputs.
        """
        jobs = []

        for readset in self.readsets:

            pycoqc_directory = os.path.join(self.output_dirs["pycoqc_directory"], readset.name)

            if readset.summary_file:
                in_summary = readset.summary_file
            else:
                _raise(SanitycheckError("Error: summary file not available for readset \"" + readset.name + "\"!"))

            align_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name)
            in_bam = os.path.join(align_directory, readset.name + ".sorted.bam")

            jobs.append(
                concat_jobs([
                    bash.mkdir(pycoqc_directory),
                    pycoqc.pycoqc(
                        readset_name=readset.name,
                        input_summary=in_summary,
                        output_directory=pycoqc_directory,
                        input_barcode=None,
                        input_bam=in_bam
                        )
                ],
                    name="pycoqc." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample.
        Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:
        Aligned and sorted BAM output files from previous minimap2_align or pbmm2_align step
        """
        jobs = []

        for sample in self.samples:

            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)

            # Find input readset BAMs first from previous minimap2_align or pbmm2_align job,
            readset_bams = self.select_input_files([
                [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets]
            ])

            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            mkdir_job = bash.mkdir(os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:

                readset_bam = readset_bams[0]

                readset_index = re.sub(r"\.bam$", ".bam.bai", readset_bam)
                sample_index = re.sub(r"\.bam$", ".bam.bai", sample_bam)

                job = concat_jobs(
                    [
                        mkdir_job,
                        bash.ln(
                            os.path.relpath(readset_bam, os.path.dirname(sample_bam)),
                            sample_bam,
                            input_file=readset_bam
                        ),
                        bash.ln(
                            os.path.relpath(readset_index, os.path.dirname(sample_index)),
                            sample_index,
                            input_file=readset_index
                        )
                    ],
                    name="symlink_readset_sample_bam." + sample.name,
                    samples=[sample],
                )
                job.samples = [sample]

            elif len(sample.readsets) > 1:

                job = concat_jobs(
                    [
                        mkdir_job,
                        gatk4.merge_sam_files(readset_bams, sample_bam)
                    ],
                    samples=[sample],
                    name="picard_merge_sam_files." + sample.name
                )

            jobs.append(job)

        return jobs

    def metrics_mosdepth(self):
        """
        Calculate depth stats with [Mosdepth](https://github.com/brentp/mosdepth)
        Returns:
            list: A list of mosdepth jobs.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)

            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")],
                ]
            )
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'], sample.name)
            mosdepth_directory = os.path.join(metrics_directory, "mosdepth")
            output_prefix = os.path.join(mosdepth_directory, sample.name)
            region = None
            output_dist = f"{output_prefix}.mosdepth.global.dist.txt"
            output_summary = f"{output_prefix}.mosdepth.summary.txt"

            job_name = f"mosdepth.{sample.name}"

            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        mosdepth.parse_dedup_coverage_metrics_pt(f"{output_prefix}.mosdepth.summary.txt"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.mosdepth.summary.txt",
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="dedup_coverage=$dedup_coverage"
                        )
                    ]
                )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(mosdepth_directory),
                        mosdepth.run(
                            input_file,
                            output_prefix,
                            True,
                            region
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[output_dist, output_summary],
                    removable_files=[]
                )
            )
            self.multiqc_inputs.extend(
                [
                    os.path.join(mosdepth_directory, os.path.basename(output_dist)),
                    os.path.join(mosdepth_directory, os.path.basename(output_summary))
                ]
            )

        return jobs
    
    def set_deepvariant_regions(self):
        """
        Create an interval list with ScatterIntervalsByNs from GATK: [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360041416072-ScatterIntervalsByNs-Picard).
        Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. 
        By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering 
        the variants with the resulting interval list will be the same as calling with one large region.
        Returns:
            list: A list of set interval list jobs.
        """
        jobs = []

        reference = global_conf.global_get('deepvariant', 'genome_fasta', param_type='filepath')
        scatter_jobs = global_conf.global_get('deepvariant', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            interval_directory = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant", "regions")
            output = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.interval_list'))
            interval_list_acgt_noalt = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            if coverage_bed:
                region = coverage_bed
                coverage_bed_noalt = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.noALT.bed'))

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(interval_directory),
                            pipe_jobs(
                                [
                                    bash.grep(
                                        region,
                                        None,
                                        '-Ev "_GL|_K"'
                                        ),
                                    bash.grep(
                                        None,
                                        coverage_bed_noalt,
                                        '-v "EBV"'
                                    )
                                ]
                            ),
                        ],
                        name=f"gatk_scatterIntervalsByNs.{sample.name}",
                        samples=[sample],
                        readsets=[*list(sample.readsets)]
                    )
                )
            elif scatter_jobs == 1:
                log.info("Number of jobs set to 1, skipping region creation for variant calling...")
  
            else:

                job = concat_jobs(
                        [
                            bash.mkdir(interval_directory),
                            gatk4.scatterIntervalsByNs(
                                reference,
                                output
                            ),
                            pipe_jobs(
                                [
                                    bash.grep(
                                        output,
                                        None,
                                        '-Ev "_GL|_K"'
                                    ),
                                    bash.grep(
                                        None,
                                        interval_list_acgt_noalt,
                                        '-v "EBV"'
                                    )
                                ]
                            ),
                            gatk4.splitInterval(
                                interval_list_acgt_noalt,
                                interval_directory,
                                scatter_jobs
                            )
                        ]
                    )
                
                for idx in range(scatter_jobs):
                    interval_file = os.path.join(interval_directory, str(idx).zfill(4) + "-scattered.interval_list")
                    bed_file = os.path.join(interval_directory, str(idx).zfill(4) + "-region.bed")
                    job = concat_jobs(
                        [
                            job,
                            gatk4.interval_list2bed(
                                interval_file,
                                bed_file,
                                ini_section='gatk_splitInterval'
                            )
                        ]
                    )

                job.name=f"gatk_scatterIntervalsByNs.{sample.name}"
                job.samples=[sample]
                job.readsets=[*list(sample.readsets)]

                jobs.append(job)

        return jobs


    def deepvariant(self):
        """
        Germline variant calling with DeepVariant.
        """
        jobs = []

        nb_jobs = global_conf.global_get('deepvariant', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            deepvariant_dir = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant")
            region_directory = os.path.join(deepvariant_dir, "regions")
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            if coverage_bed:
                region = coverage_bed
            elif nb_jobs == 1:
                region = global_conf.global_get('deepvariant', 'region') if global_conf.global_get('deepvariant', 'region') else None
            
            if nb_jobs == 1 or coverage_bed:
                
                output_vcf = os.path.join(deepvariant_dir, f"{sample.name}.deepvariant.vcf.gz")
                tmp_dir = os.path.join(deepvariant_dir, "tmp")
                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(tmp_dir),
                            deepvariant.run(
                                input_bam,
                                output_vcf,
                                tmp_dir,
                                region
                            )
                        ],
                        name=f"deepvariant.{sample.name}",
                        samples=[sample],
                        readsets=[*list(sample.readsets)]
                    )
                )
            else:
                regions = [os.path.join(region_directory, f"{idx:04d}-region.bed") for idx in range(nb_jobs)]

                for idx, region in enumerate(regions):

                    output_vcf = os.path.join(deepvariant_dir, f"{sample.name}.deepvariant.{str(idx)}.vcf.gz")
                    tmp_dir = os.path.join(deepvariant_dir, "tmp", str(idx))

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(tmp_dir),
                                deepvariant.run(
                                    input_bam,
                                    output_vcf,
                                    tmp_dir,
                                    region
                                )
                            ],
                            name=f"deepvariant.{sample.name}.{str(idx)}",
                            input_dependency=[input_bam, region],
                            samples=[sample],
                            readsets=[*list(sample.readsets)]
                        )
                    )

        return jobs
    
    def merge_filter_deepvariant(self):
        """
        Merge deepvariant outputs, if applicable, and filter vcf.
        """
        jobs = []

        nb_jobs = global_conf.global_get('deepvariant', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            
            deepvariant_dir = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant")
            deepvariant_prefix = os.path.join(deepvariant_dir, f"{sample.name}.deepvariant")
            deepvariant_vcf = os.path.join(deepvariant_dir, f"{sample.name}.deepvariant.vcf.gz")
            deepvariant_filtered = os.path.join(deepvariant_dir, f"{sample.name}.deepvariant.flt.vcf.gz")

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            job = concat_jobs(
                [
                    bcftools.view(
                        deepvariant_vcf,
                        deepvariant_filtered,
                        "-f PASS -Oz"
                    ),
                    htslib.tabix(
                        deepvariant_filtered,
                        "-f -pvcf"
                    )
                ]
            )

            if nb_jobs == 1 or coverage_bed:
                job.name = f"merge_filter_deepvariant.{sample.name}"
                job.samples = [sample]
                job.readsets = [*list(sample.readsets)]
                jobs.append(job)

            else:
                vcfs_to_merge = [f"{deepvariant_prefix}.{str(idx)}.vcf.gz" for idx in range(nb_jobs)]
                jobs.append(
                    concat_jobs(
                        [
                            bcftools.concat(
                                vcfs_to_merge,
                                deepvariant_vcf,
                                "-oZ"
                            ),
                            htslib.tabix(
                                deepvariant_vcf,
                                "-f -pvcf"
                            ),
                            job
                        ],
                        name = f"merge_filter_deepvariant.{sample.name}",
                        samples = [sample],
                        readsets = [*list(sample.readsets)],
                        removable_files=vcfs_to_merge
                    )
                )

        return jobs
    
    def svim(self):
        """
        Use SVIM to perform SV calling on each sample.
        """
        jobs = []

        for sample in self.samples:

            align_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            in_bam = os.path.join(align_directory, sample.name + ".sorted.bam")

            svim_directory = os.path.join(self.output_dirs["svim_directory"], sample.name)

            job = svim.svim_ont(in_bam, svim_directory)
            job.name = "svim." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs
    
    def sawfish(self):
        """
        Call structural variants from mapped HiFi sequencing reads with Sawfish.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            in_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            sawfish_directory = os.path.join(self.output_dirs["SVariants_directory"], sample.name, "sawfish")
            discover_directory = os.path.join(sawfish_directory, "discover")
            call_directory = os.path.join(sawfish_directory, "call")

            sawfish_output = os.path.join(call_directory, "genotyped.sv.vcf.gz")
            sawfish_output_filtered = os.path.join(sawfish_directory, f"{sample.name}.sawfish.flt.vcf.gz")
            
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(sawfish_directory),
                        sawfish.discover(
                            in_bam,
                            discover_directory
                        ),
                        sawfish.joint_call(
                            discover_directory,
                            call_directory
                        ),
                        bcftools.view(
                            sawfish_output,
                            sawfish_output_filtered,
                            "-f PASS -Oz"
                        ),
                        htslib.tabix(
                            sawfish_output_filtered,
                            "-pvcf"
                        )
                    ],
                    name=f"sawfish.{sample.name}",
                    samples=[sample],
                    output_dependency=[sawfish_output_filtered],
                    readsets=[*list(sample.readsets)]
                )
            )

        return jobs
    
    def hificnv(self):
        """
        Call copy number variation and visualise results with HiFiCNV
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            deepvariant_directory = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant")
            in_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            in_maf = os.path.join(deepvariant_directory, f"{sample.name}.deepvariant.flt.vcf.gz")

            hificnv_directory = os.path.join(self.output_dirs["SVariants_directory"], sample.name, "hificnv")
            output_prefix = os.path.join(hificnv_directory, "hificnv")

            hificnv_out = os.path.join(hificnv_directory, f"hificnv.{sample.name}.vcf.gz")
            hificnv_filtered = os.path.join(hificnv_directory, f"{sample.name}.filt.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(hificnv_directory),
                        hificnv.run(
                            in_bam,
                            output_prefix,
                            sample.name,
                            in_maf
                        ),
                        bcftools.view(
                            hificnv_out,
                            hificnv_filtered,
                            "-f PASS -Oz"
                        ),
                        htslib.tabix(
                            hificnv_filtered,
                            "-pvcf"
                        )
                    ],
                    name=f"hificnv.{sample.name}"
                )
            )

        return jobs
    
    def trgt_genotyping(self):
        """
        Call tandem repeats for pathogenic and full repeats with TRGT.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            in_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            trgt_directory = os.path.join(self.output_dirs["SVariants_directory"], sample.name, "trgt")

            pathogenic_prefix = os.path.join(trgt_directory, f"{sample.name}.pathogenic_repeats")
            full_prefix = os.path.join(trgt_directory, f"{sample.name}.full_repeats")

            pathogenic_repeats = global_conf.global_get("trgt_genotyping", 'pathogenic_repeats', required=True)
            full_repeats = global_conf.global_get("trgt_genotyping", 'full_repeat_catalog', required=True)

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(trgt_directory),
                        trgt.genotype(
                            in_bam,
                            pathogenic_repeats,
                            pathogenic_prefix
                        ),
                        bcftools.sort(
                            f"{pathogenic_prefix}.vcf.gz",
                            f"{pathogenic_prefix}.sorted.vcf.gz",
                            "-Oz"
                        ),
                        bcftools.index(
                            f"{pathogenic_prefix}.sorted.vcf.gz"
                        ),
                        samtools.sort(
                            f"{pathogenic_prefix}.spanning.bam",
                            f"{pathogenic_prefix}.spanning.sorted.bam"
                        ),
                        samtools.index(
                            f"{pathogenic_prefix}.spanning.sorted.bam"
                        )
                    ],
                    name=f"trgt_genotyping.pathogenic.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    removable_files=[
                        f"{pathogenic_prefix}.vcf.gz",
                        f"{pathogenic_prefix}.spanning.bam"
                    ]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(trgt_directory),
                        trgt.genotype(
                            in_bam,
                            full_repeats,
                            full_prefix
                        ),
                        bcftools.sort(
                            f"{full_prefix}.vcf.gz",
                            f"{full_prefix}.sorted.vcf.gz",
                            "-Oz"
                        ),
                        bcftools.index(
                            f"{full_prefix}.sorted.vcf.gz"
                        ),
                        samtools.sort(
                            f"{full_prefix}.spanning.bam",
                            f"{full_prefix}.spanning.sorted.bam"
                        ),
                        samtools.index(
                            f"{full_prefix}.spanning.sorted.bam"
                        )
                    ],
                    name=f"trgt_genotyping.full.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    removable_files=[
                        f"{full_prefix}.vcf.gz",
                        f"{full_prefix}.spanning.bam"
                    ]
                )
            )

        return jobs
    
    def annotSV(self):
        """
        Annotate and rank structural variants with AnnotSV.
        """
        jobs =[]

        for sample in self.samples:
            annotsv_directory = os.path.join(self.output_dirs["annotsv_directory"], sample.name)
            svariants_dir = os.path.join(self.output_dirs["SVariants_directory"], sample.name)
            hificnv_vcf = os.path.join(svariants_dir, "hificnv", f"{sample.name}.filt.vcf.gz")
            sawfish_vcf = os.path.join(svariants_dir, "sawfish", f"{sample.name}.sawfish.flt.vcf.gz")
            deepvariant_vcf = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant", f"{sample.name}.deepvariant.flt.vcf.gz")

            hificnv_dir = os.path.join(annotsv_directory, sample.name, "hificnv")
            sawfish_dir = os.path.join(annotsv_directory, sample.name, "sawfish")
            hificnv_annot = os.path.join(hificnv_dir, f"{sample.name}.hificnv.annotsv.tsv")
            sawfish_annot = os.path.join(sawfish_dir, f"{sample.name}.sawfish.annotsv.tsv")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(hificnv_dir),
                        annotsv.annotate(
                            hificnv_vcf,
                            hificnv_annot,
                            deepvariant_vcf
                        ),
                        annotsv.html(
                            hificnv_annot,
                            hificnv_dir,
                            f"{sample.name}.hificnv.annotsv"
                        ),
                        annotsv.excel(
                            hificnv_annot,
                            hificnv_dir,
                            f"{sample.name}.hificnv.annotsv"
                        )
                    ],
                    name=f"annotsv.hificnv.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(sawfish_dir),
                        annotsv.annotate(
                            sawfish_vcf,
                            sawfish_annot,
                            deepvariant_vcf
                        ),
                        annotsv.html(
                            sawfish_annot,
                            sawfish_dir,
                            f"{sample.name}.sawfish.annotsv"
                        ),
                        annotsv.excel(
                            sawfish_annot,
                            sawfish_dir,
                            f"{sample.name}.sawfish.annotsv"
                        )
                    ],
                    name=f"annotsv.sawfish.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
                )
            )

        return jobs
    
    def hiphase(self):
        """
        Phase variant calls with HiPhase.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            variants_dir = os.path.join(self.output_dirs["variants_directory"], sample.name)
            svariants_dir = os.path.join(self.output_dirs["SVariants_directory"], sample.name)
            deepvariant_vcf = os.path.join(variants_dir, "deepvariant", f"{sample.name}.deepvariant.flt.vcf.gz")
            sawfish_vcf = os.path.join(svariants_dir, "sawfish", f"{sample.name}.sawfish.flt.vcf.gz")
            trgt_vcf = os.path.join(svariants_dir, "trgt", f"{sample.name}.pathogenic_repeats.sorted.vcf.gz")

            hiphase_directory = os.path.join(self.output_dirs["hiphase_directory"], sample.name)
            stats_out = os.path.join(hiphase_directory, f"{sample.name}.stats.csv")
            blocks_out = os.path.join(hiphase_directory, f"{sample.name}.blocks.tsv")
            summary_out = os.path.join(hiphase_directory, f"{sample.name}.summary.tsv")

            deepvariant_out = os.path.join(hiphase_directory, f"{sample.name}.deepvariant.hiphase.vcf.gz")
            sawfish_out = os.path.join(hiphase_directory, f"{sample.name}.sawfish.hiphase.vcf.gz")
            trgt_out = os.path.join(hiphase_directory, f"{sample.name}.trgt.pathogenic_repeats.hiphase.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(hiphase_directory),
                        hiphase.run(
                            input_bam,
                            stats_out,
                            blocks_out,
                            summary_out,
                            deepvariant_vcf,
                            deepvariant_out,
                            sawfish_vcf,
                            sawfish_out,
                            trgt_vcf,
                            trgt_out
                        )
                    ],
                    name=f"hiphase.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
                )
            )

        return jobs
    
    def report_cpsr(self):
        """
        Creates a cpsr germline report (https://sigven.github.io/cpsr/)
        input: annotated/filter vcf
        output: html report and addtional flat files
        Returns:
            list: A list of cpsr report jobs.
        """
        jobs = []

        for sample in self.samples:
            hiphase_directory = os.path.join(self.output_dirs["hiphase_directory"], sample.name)
            deepvariant_phased = os.path.join(hiphase_directory, f"{sample.name}.deepvariant.hiphase.vcf.gz")
            cpsr_directory = os.path.join(self.output_dirs["report_directory"], sample.name, "cpsr")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(cpsr_directory),
                        cpsr.report(
                            deepvariant_phased,
                            cpsr_directory,
                            sample.name,
                            "report_cpsr"
                        )
                    ],
                    name=f"report_cpsr.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
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

        output = os.path.join(self.output_dirs['report_directory'], f"LongRead_DnaSeq.{self.protocol}.multiqc")

        job = concat_jobs(
            [
                bash.mkdir(os.path.join(self.output_dirs['report_directory'])),
                multiqc.run(
                    self.multiqc_inputs,
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
        return { 'nanopore': [
                self.blastqc,
                self.minimap2_align,
                self.pycoqc,
                self.picard_merge_sam_files,
                self.svim
            ], 'revio':
            [
                self.metrics_nanoplot,
                self.pbmm2_align,
                self.picard_merge_sam_files,
                self.metrics_mosdepth,
                self.set_deepvariant_regions,
                self.deepvariant,
                self.merge_filter_deepvariant,
                self.hificnv,
                self.trgt_genotyping,
                self.sawfish,
                self.annotSV,
                self.hiphase,
                self.report_cpsr,
                self.multiqc
            ]
        }

def main(parsed_args):
    """
    The function that will call this pipeline!
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    json_pt = parsed_args.json_pt
    force = parsed_args.force
    force_mem_per_cpu = parsed_args.force_mem_per_cpu
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    protocol = parsed_args.protocol
    design_file = parsed_args.design_file

    pipeline = LongReadDnaSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, protocol=protocol, design_file=design_file, json_pt=json_pt, container=container)

    pipeline.submit_jobs()
