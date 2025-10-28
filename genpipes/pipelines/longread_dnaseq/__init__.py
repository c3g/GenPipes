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
import argparse
import os
import logging
import re
import shtab

# GenPipes Modules
from ...core.config import global_conf, SanitycheckError, _raise
from ...core.job import Job, concat_jobs, pipe_jobs
from ...core.sample_tumor_pairs import parse_tumor_pair_file
from .. import common

from ...bfx import (
    annotsv,
    bash_cmd as bash,
    bcftools,
    bvatools,
    clair3,
    clairS,
    cnvkit,
    cpsr,
    deepvariant,
    djerba,
    dysgu,
    gatk4,
    hificnv,
    hiphase,
    htslib,
    job2json_project_tracking,
    minimap2,
    modkit,
    mosdepth,
    multiqc,
    nanoplot,
    pbmm2,
    pcgr,
    pycoqc,
    samtools,
    savana,
    sawfish,
    svim,
    tools,
    trgt,
    whatshap
    )

log = logging.getLogger(__name__)

class LongReadDnaSeq(common.LongRead):
    """
LongRead DNA-Seq Pipeline
==============

The LongRead Pipeline is used to analyse long reads produced by the Oxford Nanopore Technologies (ONT) 
and PacBio Revio sequencers. The protocols used are nanopore and revio, respectively. For nanopore reads, 
there is now an additional protocol called nanopore_paired_somatic that is used to analyze data from a matched 
tumor-normal pair. 

Currently, the nanopore protocol of the pipeline uses minimap2 to align reads to the reference genome.
Additionally, it produces a QC report that includes an interactive dashboard with data from the basecalling
summary file as well as the alignment. A step aligning random reads to the NCBI nt database and reporting 
the species of the highest hits is also done as QC.

Once the QC and alignments have been produced, Picard is used to merge readsets coming from the same
sample. Nanoplot and Mosdepth provide metrics at the sample level. Variant calling is performed by Clair3 and 
results are phased with Whatshap. Dysgu and SVIM are used to detect Structural Variants (SV) including deletions, 
insertions and translocations. For a full summary of the types of SVs detected by SVIM, please consult the 
following [site](https://github.com/eldariont/svim#background-on-structural-variants-and-long-reads).

The SV calls produced by SVIM and Dysgu are saved as VCFs for each sample, which can then be used in downstream
analyses. No filtering is performed on the SV calls from SVIM.

For epigenomic datasets, Modkit can be run as an optional final step. 

This pipeline currently does not perform base calling and requires both FASTQ and a sequencing_summary
file produced by a ONT supported basecaller (we recommend Guppy). Additionally, the testing and
development of the pipeline were focused on genomics applications, and functionality has not been tested
for transcriptomics datasets.

For more information on using ONT data for structural variant detection, as well as an alternative
approach, please consult [this GitHub repository](https://github.com/nanoporetech/pipeline-structural-variation).

For the nanopore_paired_somatic protocol, alignment and metrics generation follow the same steps for both the normal 
and the tumor sample. Variant calling for each sample is done with ClairS, followed by detection of somatic structural 
variants with SAVANA. Finally, CPSR and PCGR reports are created for germline and somatic variants, respectively. 

The Revio protocol uses pbmm2 to align reads to the reference genome, followed by variant calling with DeepVariant
and structural variant calling with HiFiCNV, TRGT, and Sawfish. Variants are annotated with AnnotSV and phased
with HiPhase. A CPSR report can be produced from the phased variants. Metrics on the raw and mapped reads are
collected with NanoPlot and mosdepth, respectively. 

Both protocols require as input a readset file, which provides sample metadata and paths to input data (FASTQ, FAST5 or BAM).
For information on the structure and contents of the LongRead readset file, please consult [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/readset_file.html).
    """

    def __init__(self, *args, pairs_file=None, protocol='nanopore', **kwargs):
        self.pairs = pairs_file
        self.protocol = protocol
        super(LongReadDnaSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument(
            "-t",
            "--type",
            help="Type of protocol (default nanopore)",
            dest='protocol',
            choices=["nanopore", "nanopore_paired_somatic", "revio"],
            default="nanopore"
            )
        cls._argparser.add_argument(
            "-p",
            "--pairs",
            help="pairs file",
            type=argparse.FileType('r')
            ).complete = shtab.FILE
        return cls._argparser

    @property
    def tumor_pairs(self):
        """
        The tumor_pairs attribute is a property that returns a parsed tumor pairs file.
        Returns:
            dict: A dictionary of tumor pairs.
        """
        # Create property only if somatic protocols called
        if 'somatic' in self.protocol:
            if not hasattr(self, "_tumor_pairs"):
                self._tumor_pairs = parse_tumor_pair_file(
                    self.pairs.name,
                    self.samples
                )
            return self._tumor_pairs
    
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
            'methylation_directory': os.path.relpath(os.path.join(self.output_dir, 'methylation'), self.output_dir)
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
            self._multiqc_inputs = {}
            for sample in self.samples:
                self._multiqc_inputs[sample.name] = []

            if 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
                for tumor_pair in self.tumor_pairs.values():
                    self._multiqc_inputs[tumor_pair.name] = []

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

            is_directory = False

            if readset.summary_file:
                input_summary = os.path.join(nanoplot_directory, f"{readset.name}.sequencing_summary.txt")
                link_job = bash.ln(
                    os.path.abspath(readset.summary_file),
                    input_summary,
                    readset.summary_file
                    )
                input_fastq = None
                input_bam = None
            elif readset.fastq_files:
                if os.path.isdir(readset.fastq_files):
                    is_directory = True
                input_fastq = os.path.join(nanoplot_directory, f"{readset.name}_fastq_pass")
                link_job = concat_jobs(
                    [
                        bash.rm(input_fastq),
                        bash.ln(
                            os.path.abspath(readset.fastq_files),
                            input_fastq,
                            readset.fastq_files
                        )
                    ]
                )
                input_bam = None
                input_summary = None
            elif readset.bam:
                input_bam = os.path.join(nanoplot_directory, f"{readset.name}.{os.path.basename(readset.bam)}")
                link_job = bash.ln(
                    os.path.abspath(readset.bam),
                    input_bam,
                    readset.bam
                    )
                input_fastq = None
                input_summary = None
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
                            input_fastq,
                            input_summary,
                            is_directory=is_directory
                        )
                    ],
                    name=f"nanoplot.{readset.name}",
                    samples=[readset.sample],
                    readsets=[readset],
                    input_dependency=[readset.bam, readset.fastq_files, readset.summary_file]
                )
            )

            self.multiqc_inputs[readset.sample.name].append(
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

            read_group = "'@RG" + \
                         "\\tID:" + readset.name + \
                         "\\tSM:" + readset.sample.name + \
                         "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                         ("\\tPU:run" + readset.run if readset.run else "") + \
                         "\\tPL:Nanopore" + \
                         "'"
            
            if readset.fastq_files:
                minimap2_input = readset.fastq_files
                input_dependency = readset.fastq_files
                bam2fq_job = None
            elif readset.bam_files:
                minimap2_input = None
                input_dependency = readset.bam_files
                bam2fq_job = samtools.fastq(
                        readset.bam_files,
                        "-TMM,ML"
                        )

            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))


            job = concat_jobs(
                [
                    bash.mkdir(os.path.dirname(out_bam)),
                    pipe_jobs(
                        [
                            bam2fq_job,
                            minimap2.minimap2_ont(
                                minimap2_input,
                                read_group,
                                ini_section= "minimap2_align"
                            ),
                            samtools.sort(
                                "-",
                                out_bam,
                                ini_section = 'minimap2_align'
                            )
                        ]
                    ),
                    samtools.quickcheck(
                        out_bam
                    ),
                    samtools.index(
                        out_bam,
                    )
                ],
                name="minimap2_align." + readset.name,
                samples=[readset.sample],
                input_dependency=[input_dependency]
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

    def samtools_merge_bam_files(self):
        """
        BAM readset files are merged into one file per sample.
        Merge is done using [Samtools](https://www.htslib.org/doc/samtools-merge.html).

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
                        samtools.merge(
                            sample_bam,
                            readset_bams,
                            ini_section='samtools_merge_bam_files'
                            ),
                            samtools.index(
                                sample_bam
                            )
                    ],
                    samples=[sample],
                    name=f"samtools_merge_bam_files.{sample.name}"
                )

            jobs.append(job)

        return jobs
    
    def metrics_nanoplot_aligned(self):
        """
        Collect QC metrics on aligned bam file with nanoplot.
        """
        jobs =[]

        for sample in self.samples:
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'], sample.name)
            nanoplot_directory = os.path.join(metrics_directory, "nanoplot")
            nanoplot_prefix = f"{sample.name}.aligned."

            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(nanoplot_directory),
                        nanoplot.qc(
                            nanoplot_directory,
                            nanoplot_prefix,
                            input_bam,
                            aligned=True
                        )
                    ],
                    name=f"nanoplot.aligned.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
                )
            )

            self.multiqc_inputs[sample.name].append(
                os.path.join(nanoplot_directory, f"{nanoplot_prefix}NanoStats.txt")
                )

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
            self.multiqc_inputs[sample.name].extend(
                [
                    os.path.join(mosdepth_directory, os.path.basename(output_dist)),
                    os.path.join(mosdepth_directory, os.path.basename(output_summary))
                ]
            )

        return jobs
    
    def set_variant_calling_regions(self):
        """
        Create an interval list with ScatterIntervalsByNs from GATK: [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360041416072-ScatterIntervalsByNs-Picard).
        Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. 
        By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering 
        the variants with the resulting interval list will be the same as calling with one large region.
        Returns:
            list: A list of set interval list jobs.
        """
        jobs = []

        if self.protocol == "revio":
            caller = "deepvariant"
        elif self.protocol == "nanopore_paired_somatic":
            caller = "clairS"
        else:
            caller = "clair3"

        reference = global_conf.global_get(caller, 'genome_fasta', param_type='filepath')
        scatter_jobs = global_conf.global_get(caller, 'nb_jobs', param_type='posint')

        for sample in self.samples:
            interval_directory = os.path.join(self.output_dirs["variants_directory"], sample.name, caller, "regions")
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
    
    def clair3(self):
        """
        Call germline small variants with clair3.
        """

        jobs = []

        nb_jobs = global_conf.global_get('clair3', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            align_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(align_directory, f"{sample.name}.sorted.bam")
            clair3_dir = os.path.join(self.output_dirs["variants_directory"], sample.name, "clair3")
            region_directory = os.path.join(clair3_dir, "regions")
            region_param = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            if coverage_bed:
                region = coverage_bed
                region_param = f"--bed_fn={region}"
            elif nb_jobs == 1:
                region = global_conf.global_get('clair3', 'region') if global_conf.global_get('clair3', 'region') else None
                if region:
                    if os.path.isfile(region):
                        region_param = f"--bed_fn={region}"
                    else:
                        region_param = f"--ctg_name={region}"
            
            if nb_jobs == 1 or coverage_bed:
                                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(clair3_dir),
                            clair3.run(
                                input_bam,
                                clair3_dir,
                                sample.name,
                                "ont",
                                region_param
                            )
                        ],
                        name=f"clair3.{sample.name}",
                        samples=[sample],
                        readsets=[*list(sample.readsets)]
                    )
                )
            else:
                regions = [os.path.join(region_directory, f"{idx:04d}-region.bed") for idx in range(nb_jobs)]

                for idx, region in enumerate(regions):

                    output_dir = os.path.join(clair3_dir, str(idx))

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(output_dir),
                                clair3.run(
                                    input_bam,
                                    output_dir,
                                    sample.name,
                                    "ont",
                                    f"--bed_fn={region}"
                                )
                            ],
                            name=f"clair3.{sample.name}.{str(idx)}",
                            input_dependency=[input_bam, region],
                            samples=[sample],
                            readsets=[*list(sample.readsets)]
                        )
                    )

        return jobs
    
    def merge_filter_clair3(self):
        """
        Merge clair3 outputs, if applicable, and filter vcf.
        """
        jobs = []

        nb_jobs = global_conf.global_get('clair3', 'nb_jobs', param_type='posint')

        for sample in self.samples:
            
            clair3_dir = os.path.join(self.output_dirs["variants_directory"], sample.name, "clair3")
            clair3_vcf = os.path.join(clair3_dir, "phased_merge_output.vcf.gz")
            clair3_filtered = os.path.join(clair3_dir, f"{sample.name}.clair3.phased.flt.vcf.gz")

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            job = concat_jobs(
                [
                    bcftools.view(
                        clair3_vcf,
                        clair3_filtered,
                        "-f PASS -Oz"
                    ),
                    htslib.tabix(
                        clair3_filtered,
                        "-f -pvcf"
                    )
                ]
            )

            if nb_jobs == 1 or coverage_bed:
                job.name = f"merge_filter_clair3.{sample.name}"
                job.samples = [sample]
                job.readsets = [*list(sample.readsets)]
                jobs.append(job)

            else:
                vcfs_to_merge = [os.path.join(clair3_dir, str(idx), "phased_merge_output.vcf.gz") for idx in range(nb_jobs)]
                jobs.append(
                    concat_jobs(
                        [
                            bcftools.concat(
                                vcfs_to_merge,
                                clair3_vcf,
                                "-a -oZ"
                            ),
                            htslib.tabix(
                                clair3_vcf,
                                "-f -pvcf"
                            ),
                            job
                        ],
                        name = f"merge_filter_clair3.{sample.name}",
                        samples = [sample],
                        readsets = [*list(sample.readsets)],
                        removable_files=vcfs_to_merge
                    )
                )

        return jobs
    
    def clairS(self):
        """
        Call somatic small variants with clairS.
        """

        jobs = []

        nb_jobs = global_conf.global_get('clairS', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            normal_align_directory = os.path.join(self.output_dirs["alignment_directory"], tumor_pair.normal.name)
            tumor_align_directory = os.path.join(self.output_dirs["alignment_directory"], tumor_pair.tumor.name)
            normal_bam = os.path.join(normal_align_directory, f"{tumor_pair.normal.name}.sorted.bam")
            tumor_bam = os.path.join(tumor_align_directory, f"{tumor_pair.tumor.name}.sorted.bam")
            clairS_dir = os.path.join(self.output_dirs["variants_directory"], tumor_pair.name, "clairS")
            region_directory = os.path.join(self.output_dirs["variants_directory"], tumor_pair.normal.name, "clairS", "regions")
            region_param = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if coverage_bed:
                region = coverage_bed
                region_param = f"--bed_fn={region}"
            elif nb_jobs == 1:
                region = global_conf.global_get('clairS', 'region') if global_conf.global_get('clairS', 'region') else None
                if region:
                    if os.path.isfile(region):
                        region_param = f"--bed_fn={region}"
                    else:
                        region_param = f"--ctg_name={region}"
            
            if nb_jobs == 1 or coverage_bed:
                                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(clairS_dir),
                            clairS.run(
                                normal_bam,
                                tumor_bam,
                                clairS_dir,
                                tumor_pair.tumor.name,
                                "ont_r10_dorado_sup_5khz_ssrs",
                                region_param
                            )
                        ],
                        name=f"clairS.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.readsets)]
                    )
                )
            else:
                regions = [os.path.join(region_directory, f"{idx:04d}-region.bed") for idx in range(nb_jobs)]

                for idx, region in enumerate(regions):

                    output_dir = os.path.join(clairS_dir, str(idx))

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(output_dir),
                                clairS.run(
                                    normal_bam,
                                    tumor_bam,
                                    output_dir,
                                    tumor_pair.tumor.name,
                                    "ont",
                                    f"--bed_fn={region}"
                                )
                            ],
                            name=f"clairS.{tumor_pair.name}.{str(idx)}",
                            input_dependency=[normal_bam, tumor_bam, region],
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs
    
    def merge_filter_clairS(self):
        """
        Merge clairS outputs and filter vcf.
        Germline and somatic outputs are merged for downstream use in CPSR/PCGR, respectively.
        """
        jobs = []

        nb_jobs = global_conf.global_get('clairS', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            
            clairS_dir = os.path.join(self.output_dirs["variants_directory"], tumor_pair.name, "clairS")
            clairS_germline_vcf = os.path.join(clairS_dir, "clair3_normal_germline_output.vcf.gz")
            clairS_germline_filtered = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.germline.flt.vcf.gz")
            clairS_somatic_vcf = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.somatic.vcf.gz")
            clairS_somatic_filtered = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.somatic.flt.vcf.gz")

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if nb_jobs == 1 or coverage_bed:

                clairS_indel = os.path.join(clairS_dir, "indel.vcf.gz")
                clairS_snv = os.path.join(clairS_dir, "snv.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bcftools.reheader(
                                        clairS_germline_vcf,
                                        None,
                                        f"-n {tumor_pair.normal.name}",
                                        ini_section="merge_filter_clairS"
                                    ),
                                    bcftools.view(
                                        None,
                                        clairS_germline_filtered,
                                        "-f PASS -Oz",
                                        ini_section="merge_filter_clairS"
                                    )
                                ]
                            ),
                            htslib.tabix(
                                clairS_germline_filtered,
                                "-f -pvcf"
                            ),
                        ],
                        name=f"merge_filter_clairS.germline.{tumor_pair.name}",
                        samples=[tumor_pair.normal],
                        readsets=[*list(tumor_pair.normal.readsets)]
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bcftools.concat(
                                [clairS_indel, clairS_snv],
                                clairS_somatic_vcf
                            ),
                            htslib.tabix(
                                clairS_somatic_vcf,
                                "-f -pvcf"
                            ),
                            bcftools.view(
                                clairS_somatic_vcf,
                                clairS_somatic_filtered,
                                "-f PASS -Oz"
                            ),
                            htslib.tabix(
                                clairS_somatic_filtered,
                                "-f -pvcf"
                            )
                        ],
                        name = f"merge_filter_clairS.somatic.{tumor_pair.name}",
                        samples=[tumor_pair.tumor],
                        readsets=[*list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                germline_vcfs_to_merge = [os.path.join(clairS_dir, str(idx), "clair3_normal_germline_output.vcf.gz") for idx in range(nb_jobs)]
                somatic_vcfs_to_merge = [os.path.join(clairS_dir, str(idx), "indel.vcf.gz") for idx in range(nb_jobs)]
                somatic_vcfs_to_merge.extend(
                    [os.path.join(clairS_dir, str(idx), "snv.vcf.gz") for idx in range(nb_jobs)]
                    )

                jobs.append(
                    concat_jobs(
                        [
                            bcftools.concat(
                                germline_vcfs_to_merge,
                                clairS_germline_vcf,
                                "-a -oZ"
                            ),
                            htslib.tabix(
                                clairS_germline_vcf,
                                "-f -pvcf"
                            ),
                            pipe_jobs(
                                [
                                    bcftools.reheader(
                                        clairS_germline_vcf,
                                        None,
                                        f"-n {tumor_pair.normal.name}"
                                    ),
                                    bcftools.view(
                                        None,
                                        clairS_germline_filtered,
                                        "-f PASS -Oz"
                                    )
                                ]
                            ),
                            htslib.tabix(
                                clairS_germline_filtered,
                                "-f -pvcf"
                            )
                        ],
                        name = f"merge_filter_clairS.germline.{tumor_pair.name}",
                        samples = [tumor_pair.normal],
                        readsets = [*list(tumor_pair.normal.readsets)],
                        removable_files=germline_vcfs_to_merge
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bcftools.concat(
                                somatic_vcfs_to_merge,
                                clairS_somatic_vcf,
                                "-a -oZ"
                            ),
                            htslib.tabix(
                                clairS_somatic_vcf,
                                "-f -pvcf"
                            ),
                            bcftools.view(
                                clairS_somatic_vcf,
                                clairS_somatic_filtered,
                                "-f PASS -Oz"
                            ),
                            htslib.tabix(
                                clairS_somatic_filtered,
                                "-f -pvcf"
                            )
                        ],
                        name = f"merge_filter_clairS.somatic.{tumor_pair.name}",
                        samples = [tumor_pair.tumor],
                        readsets = [*list(tumor_pair.tumor.readsets)],
                        removable_files=somatic_vcfs_to_merge
                    )
                )

        return jobs
    
    def savana(self):
        """
        Call somatic structural variants and copy number aberrations with Savana.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            normal_align_directory = os.path.join(self.output_dirs["alignment_directory"], tumor_pair.normal.name)
            tumor_align_directory = os.path.join(self.output_dirs["alignment_directory"], tumor_pair.tumor.name)
            normal_bam = os.path.join(normal_align_directory, f"{tumor_pair.normal.name}.sorted.bam")
            tumor_bam = os.path.join(tumor_align_directory, f"{tumor_pair.tumor.name}.sorted.bam")
            clairS_dir = os.path.join(self.output_dirs["variants_directory"], tumor_pair.name, "clairS")
            clairS_germline_vcf = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.germline.flt.vcf.gz")
            output_directory = os.path.join(self.output_dirs['SVariants_directory'], tumor_pair.name, 'savana')

            jobs.append(
                concat_jobs(
                    [
                        bash.rm(output_directory),
                        bash.mkdir(output_directory),
                        savana.run(
                            normal_bam,
                            tumor_bam,
                            output_directory,
                            tumor_pair.name,
                            clairS_germline_vcf
                        )
                    ],
                    name=f"savana.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[normal_bam, tumor_bam, clairS_germline_vcf]
                )
            )

        return jobs

    def whatshap(self):
        """
        Create a haplotagged file using Whatshap.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            clair3_dir = os.path.join(self.output_dirs["variants_directory"], sample.name, "clair3")
            clair3_filtered = os.path.join(clair3_dir, f"{sample.name}.clair3.phased.flt.vcf.gz")
            output_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.haplotag.bam")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        whatshap.haplotag(
                            input_bam,
                            clair3_filtered,
                            output_bam
                        )
                    ],
                    name=f"whatshap.{sample.name}",
                    samples=[sample],
                    readsets = [*list(sample.readsets)]
                )
            )

        return jobs
    
    def qdnaseq(self):
        """
        Run QDNAseq R script.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            output_dir = os.path.join(self.output_dirs["SVariants_directory"], sample.name, "QDNAseq")
            bin_size = global_conf.global_get("qdnaseq", "bin_size")
            reference = global_conf.global_get("qdnaseq", "reference")
            qdnaseq_vcf = os.path.join(output_dir, f"{sample.name}.CNV_calls.{bin_size}k.{reference}.vcf")
            qdnaseq_filtered =  os.path.join(output_dir, f"{sample.name}.CNV_calls.{bin_size}k.{reference}.filtered.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_dir),
                        tools.r_qdna_seq(
                            input_bam,
                            output_dir,
                            sample.name,
                            bin_size
                        ),
                        bcftools.view(
                            qdnaseq_vcf,
                            qdnaseq_filtered,
                        "-f PASS -Oz"
                        ),
                        htslib.tabix(
                            qdnaseq_filtered,
                            "-f -pvcf"
                        )
                    ],
                    name=f"qdnaseq.{sample.name}",
                    samples=[sample],
                    readsets = [*list(sample.readsets)],
                    removable_files=[qdnaseq_vcf]
                )
            )

        return jobs
    
    def dysgu(self):
        """
        Call structural variants with dysgu.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.haplotag.bam")
            output_dir = os.path.join(self.output_dirs["SVariants_directory"], sample.name, "dysgu")
            output_vcf = os.path.join(output_dir, f"{sample.name}.dysgu.vcf")
            dysgu_filtered = os.path.join(output_dir, f"{sample.name}.dysgu.filtered.vcf.gz")
            
            region = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                region = coverage_bed
            elif global_conf.global_get('dysgu', 'region'):
                region = global_conf.global_get('dysgu', 'region')

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_dir),
                        dysgu.call(
                            input_bam,
                            output_vcf,
                            region
                        ),
                        bcftools.view(
                            output_vcf,
                            dysgu_filtered,
                            "-f PASS -Oz"
                        ),
                        htslib.tabix(
                            dysgu_filtered,
                            "-f -pvcf"
                        )
                    ],
                    name=f"dysgu.{sample.name}",
                    samples=[sample],
                    readsets = [*list(sample.readsets)],
                    removable_files=[output_vcf]
                )
            )

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
            hificnv_filtered = os.path.join(hificnv_directory, f"{sample.name}.hificnv.filt.vcf.gz")

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
            hificnv_vcf = os.path.join(svariants_dir, "hificnv", f"{sample.name}.hificnv.filt.vcf.gz")
            sawfish_vcf = os.path.join(svariants_dir, "sawfish", f"{sample.name}.sawfish.flt.vcf.gz")
            deepvariant_vcf = os.path.join(self.output_dirs["variants_directory"], sample.name, "deepvariant", f"{sample.name}.deepvariant.flt.vcf.gz")

            hificnv_dir = os.path.join(annotsv_directory, "hificnv")
            sawfish_dir = os.path.join(annotsv_directory, "sawfish")
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

        # Set directory, ini_section, job and sample name for nanopore tumor only protocol
        if 'tumor_only' in self.protocol:
            output_directory = os.path.join(self.output_dirs['variants_directory'], "split")
            ini_section = 'report_cpsr_tumor_only'

            for sample in self.samples:
                job_name = f"report_cpsr_tumor_only.{sample.name}"
                samples = [sample]

                input_file = os.path.join(
                    output_directory,
                    sample.name,
                    f"{sample.name}.annot.vcf.gz"
                )
                cpsr_directory = os.path.join(output_directory, sample.name, "cpsr")
                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cpsr_directory,
                            ),
                            cpsr.report(
                                input_file,
                                cpsr_directory,
                                sample.name,
                                ini_section=ini_section
                            )
                        ],
                        name=job_name,
                        samples=samples,
                        readsets=[*list(sample.readsets)],
                    )
                )

        elif 'somatic' in self.protocol:
            for tumor_pair in self.tumor_pairs.values():
                clairS_dir = os.path.join(self.output_dirs["variants_directory"], tumor_pair.name, "clairS")
                input_file = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.germline.flt.vcf.gz")
                cpsr_directory = os.path.join(self.output_dirs["report_directory"], tumor_pair.name, "cpsr")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cpsr_directory,
                            ),
                            cpsr.report(
                                input_file,
                                cpsr_directory,
                                tumor_pair.normal.name,
                                "report_cpsr"
                            )
                        ],
                        name=f"report_cpsr.{tumor_pair.name}",
                        samples=[tumor_pair.normal],
                        readsets=[*list(tumor_pair.normal.readsets)],
                    )
                )
                
        #Set directory, ini_section, job and sample name for revio protocol
        else:
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
    
    def report_pcgr(self):
        """
        Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
        input: filtered somatic vcf
        output: html report and addtional flat files
        Returns:
            list: A list of pcgr report jobs.
        """
        jobs = []

        # Set directory, ini_section, job and sample name for dnaseq tumor only protocol
        if 'tumor_only' in self.protocol:
            output_directory = os.path.join(self.output_dirs['variants_directory'], "split")
            ini_section = 'report_pcgr_tumor_only'
            assembly = global_conf.global_get(ini_section, 'assembly')

            for sample in self.samples:
                cpsr_directory = os.path.join(
                    output_directory,
                    sample.name,
                    "cpsr"
                )
                input_cpsr = os.path.join(
                    cpsr_directory,
                    f"{sample.name}.cpsr.{assembly}.json.gz"
                )
                input_file = os.path.join(
                    output_directory,
                    sample.name,
                    f"{sample.name}.annot.vcf.gz"
                )

                input_cna = os.path.join(
                    self.output_dirs['sv_variants_directory'],
                    f"{sample.name}.cnvkit.cna.tsv"
                )

                pcgr_directory = os.path.join(
                    output_directory,
                    sample.name,
                    "pcgr"
                )

                input_cpsr = os.path.join(
                        cpsr_directory,
                        f"{sample.name}.cpsr.{assembly}"
                    )
                output = os.path.join(
                        pcgr_directory,
                        f"{sample.name}.pcgr.{assembly}.html"
                    )
                input_dependencies = [input_file, input_cpsr + ".classification.tsv.gz", input_cpsr + ".conf.yaml", input_cna]
                job_name = f"report_pcgr_tumor_only.{sample.name}"

                pcgr_job = concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory,
                        ),
                        pcgr.report(
                            input_file,
                            input_cpsr,
                            pcgr_directory,
                            sample.name,
                            input_cna,
                            ini_section=ini_section
                        ),
                        bash.ls(output)
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    input_dependency=input_dependencies,
                    output_dependency=[output]
                )

                if self.project_tracking_json:
                    samples = [sample]
                    pcgr_output_file = os.path.join(self.output_dir, "job_output", "report_pcgr", f"{job_name}_$TIMESTAMP.o")
                    jobs.append(
                        concat_jobs(
                            [
                                pcgr_job,
                                pcgr.parse_pcgr_passed_variants_pt(pcgr_output_file),
                                job2json_project_tracking.run(
                                    input_file=pcgr_output_file,
                                    samples=",".join([sample.name for sample in samples]),
                                    readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                                    job_name=job_name,
                                    metrics="pcgr_passed_variants=$pcgr_passed_variants"
                                )
                            ],
                            name=job_name,
                            samples=[sample],
                            readsets=[*list(sample.readsets)],
                            input_dependency=input_dependencies,
                            output_dependency=[output]
                        )
                    )
                else:
                    jobs.append(pcgr_job)

        else:
            for tumor_pair in self.tumor_pairs.values():
                # Set directory, ini_section, job and sample name for nanopore somatic protocol
 
                ini_section = 'report_pcgr'
                assembly = global_conf.global_get(ini_section, 'assembly')
                job_name = f"report_pcgr.{tumor_pair.name}"

                cpsr_directory = os.path.join(self.output_dirs["report_directory"], tumor_pair.name, "cpsr")
                input_cpsr = os.path.join(cpsr_directory, f"{tumor_pair.normal.name}.cpsr.{assembly}")
                pcgr_directory = os.path.join(self.output_dirs["report_directory"], tumor_pair.name, "pcgr")

                savana_directory = os.path.join(self.output_dirs['SVariants_directory'], tumor_pair.name, 'savana')
                output_savana = os.path.join(savana_directory, f"{tumor_pair.name}_segmented_absolute_copy_number.tsv")
                input_cna = os.path.join(savana_directory, f"{tumor_pair.name}.savana.cna.tsv")
                clairS_dir = os.path.join(self.output_dirs["variants_directory"], tumor_pair.name, "clairS")
                output_clairS = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.somatic.flt.vcf.gz")
                input_vcf = os.path.join(clairS_dir, f"{tumor_pair.name}.clairS.somatic.flt.pcgr.vcf.gz")

                output_report = os.path.join(pcgr_directory, f"{tumor_pair.name}.pcgr.{assembly}.html")
                output_maf = os.path.join(pcgr_directory, f"{tumor_pair.name}.pcgr.{assembly}.maf")
                input_dependencies = [input_vcf, input_cpsr + ".classification.tsv.gz", input_cpsr + ".conf.yaml", input_cna]

                format_savana_job = concat_jobs(
                    [
                        tools.savana2cnvkit(
                            output_savana,
                            input_cna
                        ),
                        cnvkit.file_check(
                                input_cna,
                                f"{input_cna}.pass"
                            )
                    ],
                    name = f"savana2cnvkit.{tumor_pair.name}",
                    samples = [tumor_pair.normal, tumor_pair.tumor],
                    readsets = [*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )

                jobs.append(format_savana_job)
                
                format_clairS_job = concat_jobs(
                    [
                        tools.format2pcgr(
                            output_clairS,
                            input_vcf,
                            1,
                            "somatic",
                            tumor_pair.tumor.name,
                            ini_section="format2pcgr"
                        ),
                        htslib.tabix(
                            input_vcf,
                            "-f -pvcf"
                        )
                    ],
                    name=f"format2pcgr.{tumor_pair.name}",
                    samples = [tumor_pair.normal, tumor_pair.tumor],
                    readsets = [*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],

                )

                jobs.append(format_clairS_job)
                
                pcgr_job = concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory
                        ),
                        pcgr.report(
                            input_vcf,
                            input_cpsr,
                            pcgr_directory,
                            tumor_pair.name,
                            input_cna,
                            ini_section=ini_section
                        ),
                        bash.ls(output_report)
                    ],
                    name=job_name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=input_dependencies,
                    output_dependency=[output_report, output_maf]
                )

                samples = [tumor_pair.normal, tumor_pair.tumor]

                if self.project_tracking_json:
                    pcgr_output_file = os.path.join(self.output_dir, "job_output", "report_pcgr", f"{job_name}_$TIMESTAMP.o")
                    jobs.append(
                        concat_jobs(
                            [
                                pcgr_job,
                                pcgr.parse_pcgr_passed_variants_pt(pcgr_output_file),
                                job2json_project_tracking.run(
                                    input_file=pcgr_output_file,
                                    samples=",".join([sample.name for sample in samples]),
                                    readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                                    job_name=job_name,
                                    metrics="pcgr_passed_variants=$pcgr_passed_variants"
                                )
                            ],
                            name=job_name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            input_dependency=input_dependencies,
                            output_dependency=[output_report, output_maf]
                        )
                    )
                else:
                    jobs.append(pcgr_job)

        return jobs
    
    def report_djerba(self):
        """
        Produce Djerba report.
        """
        jobs = []
        
        token = global_conf.global_get('report_djerba', 'oncokb_token', param_type='filepath', required=False)

        if token:
            assembly = global_conf.global_get('report_pcgr', 'assembly')
        
            for tumor_pair in self.tumor_pairs.values():
                djerba_dir = os.path.join(self.output_dirs['report_directory'], tumor_pair.name, "djerba")
                pcgr_directory = os.path.join(self.output_dirs["report_directory"], tumor_pair.name, "pcgr")
                input_maf = os.path.join(pcgr_directory, tumor_pair.name + ".pcgr." + assembly + ".maf")
                clean_maf =  os.path.join(djerba_dir, tumor_pair.name + ".pcgr." + assembly + ".clean.maf")
                config_file = os.path.join(djerba_dir, tumor_pair.name + ".djerba.ini")
                djerba_script = os.path.join(djerba_dir, "djerba_report." + tumor_pair.name + ".sh")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(djerba_dir),
                            djerba.clean_maf(
                                input_maf,
                                clean_maf
                                ),
                            djerba.make_config(
                                config_file,
                                tumor_pair.name,
                                tumor_pair.tumor.name,
                                tumor_pair.normal.name,
                                clean_maf + ".gz",
                                None,
                                "WGS"
                                ),
                            # djerba report requires internet connection. Script is produced but must be executed locally.
                            djerba.make_script(
                                config_file,
                                djerba_dir,
                                djerba_script
                                )
                        ],
                        name="report_djerba." + tumor_pair.name,
                        samples=[tumor_pair.tumor],
                        readsets=list(tumor_pair.tumor.readsets),
                        input_dependency=[input_maf],
                        output_dependency=[config_file, djerba_script]
                        )
                    )

        else:
            log.debug("No OncoKB token provided in config file, skipping djerba report step.")

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
        multiqc_files_paths = [item for sample in self.samples for item in self.multiqc_inputs[sample.name]]

        job = concat_jobs(
            [
                bash.mkdir(os.path.join(self.output_dirs['report_directory'])),
                multiqc.run(
                    multiqc_files_paths,
                    output
                )
            ]
        )
        job.name = "multiqc"
        job.input_files = multiqc_files_paths
        jobs.append(job)

        return jobs
    
    def modkit(self):
        """
        Methylation analysis for nanopore data.
        """
        
        jobs = []

        for sample in self.samples:
            alignment_dir = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_bam = os.path.join(alignment_dir, f"{sample.name}.sorted.bam")
            output_dir = os.path.join(self.output_dirs['methylation_directory'], sample.name)
            output_bed = os.path.join(output_dir, f"{sample.name}.cpg.pileup.bed")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_dir),
                        modkit.pileup(
                            input_bam,
                            output_bed
                        )
                    ],
                    name=f"modkit.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return { 'nanopore': [
                self.blastqc,
                self.metrics_nanoplot,
                self.minimap2_align,
                self.pycoqc,
                self.samtools_merge_bam_files,
                self.metrics_nanoplot_aligned,
                self.metrics_mosdepth,
                self.set_variant_calling_regions,
                self.clair3,
                self.merge_filter_clair3,
                self.whatshap,
                self.qdnaseq,
                self.dysgu,
                self.svim,
                self.multiqc,
                self.modkit
            ], 'nanopore_paired_somatic': [
                self.blastqc,
                self.metrics_nanoplot,
                self.minimap2_align,
                #self.pycoqc,
                self.samtools_merge_bam_files,
                self.metrics_nanoplot_aligned,
                self.metrics_mosdepth,
                self.set_variant_calling_regions,
                self.clairS,
                self.merge_filter_clairS,
                self.savana,
                self.report_cpsr,
                self.report_pcgr,
                self.report_djerba,
                self.multiqc
            ], 'revio':
            [
                self.metrics_nanoplot,
                self.pbmm2_align,
                self.samtools_merge_bam_files,
                self.metrics_nanoplot_aligned,
                self.metrics_mosdepth,
                self.set_variant_calling_regions,
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
    pairs_file = parsed_args.pairs

    pipeline = LongReadDnaSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, protocol=protocol, design_file=design_file, json_pt=json_pt, container=container, pairs_file=pairs_file)

    pipeline.submit_jobs()
