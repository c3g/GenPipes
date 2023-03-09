#!/usr/bin/env python

######
# Copyright etc
#####

# Python Standard Modules
import logging
import os
import re
import sys
import argparse

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
import utils.utils
from pipelines import common

from bfx import agent
from bfx import bwa
from bfx import samtools
from bfx import vardict
from bfx import hmm
from bfx import ichorCNA

from bfx import bash_cmd as bash
from core.config import config
from core.sample_tumor_pairs import parse_tumor_pair_file # using tumor pair system for ichorCNA step for now

log = logging.getLogger(__name__)

class DOvEE_gene(common.Illumina):
    """
    DOvEE Gene Pipeline
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments?
        self.argparser.add_argument("-p", "--pairs", help="File with sample pairing information", type=argparse.FileType('r')) # only needed for copy number protocol ichorCNA step
        self.argparser.add_argument("-t", "--type", help="Type of pipeline (default vardict)", choices=["vardict", "copy-number"], default="vardict") # FINAL NAMES TBD
        super(DOvEE_gene, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
                'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
                'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
                'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
                'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
                'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
                'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
                'wig_directory' : os.path.relpath(os.path.join(self.output_dir, 'wig'), self.output_dir), #temp name
                'cna_directory' : os.path.relpath(os.path.join(self.output_dir, 'cna'), self.output_dir) #temp name
                }
        return dirs

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(
                    self.args.pairs.name,
                    self.samples,
                    )
        return self._tumor_pairs

    def trimmer(self):
        """
        Read trimming with AGeNT Trimmer.
        """

        jobs = []

        for readset in self.readsets:
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim")
            fastq1 = readset.fastq1
            fastq2 = readset.fastq2
            
            job = agent.trimmer(
                fastq1,
                fastq2,
                trim_file_prefix
                )

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(trim_directory),
                            job,
                            bash.mv(trim_file_prefix + "_R1.fastq.gz", trim_file_prefix + ".pair1.fastq.gz"), # trimmer only allows specifying outfile prefix to which it adds _R1/2.fastq.gz
                            bash.mv(trim_file_prefix + "_R2.fastq.gz", trim_file_prefix + ".pair2.fastq.gz"), # can either rename here or change following to expect that naming
                        ], name="agent_trimmer." + readset.name
                        )
                 )
        return jobs

    def bwa_mem_samtools_sort(self):
        """
        The trimmed reads are aligned to the reference genome with bwa mem, followed by sorting with samtools.
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name + ".sorted")
            
            fastq1 = ""
            fastq2 = ""
            # Find input readset FASTQs first from previous trimmer job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        pipe_jobs(
                            [
                                bwa.mem(
                                    fastq1,
                                    fastq2,
                                    read_group="'@RG" + \
                                        "\\tID:" + readset.name + \
                                        "\\tSM:" + readset.sample.name + \
                                        "\\tLB:" + readset.library + \
                                        "\\tPU:" + readset.library + \
                                        ("\\tCN:" + config.param('bwa_mem', 'sequencing_center')) + \
                                        ("\\tPL:" + config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem', 'sequencing_technology', required=False) else "\\tPL:Illumina") + \
                                        "'",
                                    ini_section="bwa_mem"
                                ),
                                samtools.view(
                                    "-",
                                    None,
                                    "-b "
                                ),
                                samtools.sort(
                                    "-",
                                    readset_bam,
                                )
                            ]
                        )
                    ],
                    name="bwa_mem." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs

    def samtools_merge(self):
        """
        Merges sorted bam files for each sample, if there are multiple readsets for the sample. If only a single readset file is found, a symlink to the bam is created. 
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            readset_bams = [os.path.join(alignment_directory, readset.name + ".sorted.bam") for readset in sample.readsets]
            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]

                job = bash.ln(
                       os.path.relpath(readset_bam, os.path.dirname(sample_bam)),
                       sample_bam,
                       input=readset_bam
                        )
                job.name = "symlink_readset_sample_bam." + sample.name
                job.samples = [sample]

                jobs.append(job)

            elif len(sample.readsets) > 1:
                job = samtools.merge(
                        sample_bam,
                        readset_bams
                        )
                job.name = "samtools_merge." + sample.name
                job.samples = [sample]

                jobs.append(job)

        return jobs
                        
    def locatit_dedup_bam(self):
        """
        Deduplicated bam files with AGeNT locatit in hybrid or duplex mode.
        """

        jobs = []

        for sample in self.samples:
           alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
           input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

           if 'saliva' in sample.name:
               output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
               output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
               covered_bed = config.param('agent_locatit', 'covered_bedv8', param_type='filepath')

               jobs.append(
                       concat_jobs(
                           [
                               agent.locatit(
                                   input_bam,
                                   output_duplex,
                                   covered_bed,
                                   "duplex"
                                   ),
                               agent.locatit(
                                   input_bam,
                                   output_hybrid,
                                   covered_bed,
                                   "hybrid"
                                   )
                            ],
                           name='agent_locatit_dedup.' + sample.name
                           )
                       )
           elif 'brush' in sample.name:
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                covered_bed = config.param('agent_locatit', 'covered_bedv7', param_type='filepath')

                job = agent.locatit(
                        input_bam,
                        output_duplex,
                        covered_bed,
                        "duplex"
                        )
                job.name='agent_locatit_dedup.' + sample.name
                jobs.append(job)
                
        return jobs
    

    def samtools_sort(self):
        """
        Sort deduplicated bams by coordinate with samtools.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            
            if 'saliva' in sample.name: #temporary solution
                input_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                input_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
                output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.sorted")
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted")

                jobs.append(
                        concat_jobs(
                            [
                                samtools.sort(
                                    input_duplex,
                                    output_duplex
                                    ),
                                samtools.sort(
                                    input_hybrid,
                                    output_hybrid
                                    )
                            ],
                            name = "samtools.sort." + sample.name
                        )
                    )
            elif 'brush' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted")

                job = samtools.sort(
                        input_duplex,
                        output_duplex
                        )
                job.name = "samtools.sort." + sample.name
                jobs.append(job)

        return jobs

    def samtools_index(self):
        """
        Index deduplicated and sorted bams with samtools.
        """

        jobs=[]

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            if 'saliva' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")
                input_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.sorted.bam")

                jobs.append(
                        concat_jobs(
                            [
                                samtools.index(
                                    input_duplex
                                    ),
                                samtools.index(
                                    input_hybrid
                                    )
                            ],
                            name = "samtools.index." + sample.name
                            )
                        )
            elif 'brush' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")

                job = samtools.index(
                            input_duplex
                        )
                job.name = "samtools.index." + sample.name
                jobs.append(job)

        return jobs

    def vardict_single(self):
        """
        Variant calling with vardict.
        """

        jobs = []
        
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            variants_directory = os.path.join(self.output_dirs['variants_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam") # which bam should be used here? For saliva two bams exist
            output = os.path.join(variants_directory, sample.name + "out.vcf") # temp name

            if 'brush' in sample.name:
                freq=0.001
                region=config.param('vardict_single', 'target_file', param_type='filepath')  #target bed file
                nosv=True

            elif 'saliva' in sample.name:
                freq=0.1
                region=config.param('vardict_single', 'target_filev7', param_type='filepath')  #target bed file
                nosv=False

            vardict_job = vardict.single_java(
                    input_bam,
                    sample.name,
                    None,
                    nosv=nosv,
                    freq=freq,
                    region=region
                    )

            teststrandbias_job = vardict.teststrandbias(
                    None,
                    None
                    )

            var2vcf_job = vardict.var2vcf_valid(
                    output,
                    sample.name,
                    freq,
                    None
                    )

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(variants_directory),
                            pipe_jobs(
                                [
                                    vardict_job,
                                    teststrandbias_job,
                                    var2vcf_job
                                    ]
                                )
                            ], name = "vardict." + sample.name
                        )
                    )

        return jobs

    def hmm_readCounter(self):
        """
        Counting number of reads in non-overlapping windows of fixed width directly from BAM files with HMM Copy Utils readCounter:
        https://github.com/shahcompbio/hmmcopy_utils
        """
        # Which bam file is used as input here? Need to clarify, whether duplex or hybrid dedup, etc.
        # Run on both saliva and brush, as both wigs needed as input for IchorCNA

        jobs = []

        for sample in self.samples: # or iterate over sample_pair/tumor_pair name?
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            wig_directory = os.path.join(self.output_dirs['wig_directory'], sample.name) 
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam") # how should bam be generated here? For locatit, two dedup options exist. Sorting needed?
            output = os.path.join(wig_directory, sample.name + "out.wig") # temp name

            jobs.append(
                    concat_jobs(
                    [
                        bash.mkdir(wig_directory),
                        hmm.readCounter(
                            input_bam,
                            output
                            )
                    ], name = "hmm_readCounter." + sample.name
                )
            )

        return jobs

    def run_ichorCNA(self):
        """
        https://github.com/broadinstitute/ichorCNA
        'Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.'
        """
        # Run once per patient using paired saliva and brush sample wigs, requires knowing which saliva and brush came from same patient.

        jobs = []
        
        for sample_pair in self.tumor_pairs.values(): # or another id that allows us to pair brush and saliva, using tumor pair system for now?
            wig_directory = self.output_dirs['wig_directory']
            input_brush = os.path.join(wig_directory, sample_pair.tumor.name, sample_pair.tumor.name + "out.wig") # temp name 
            input_saliva = os.path.join(wig_directory, sample_pair.normal.name, sample_pair.normal.name + "out.wig") # temp name 
            output_dir = os.path.join(self.output_dirs['cna_directory'], sample_pair.name)

            jobs.append(
                    concat_jobs(
                        [
                        bash.mkdir(output_dir),
                        ichorCNA.run_ichorCNA(
                            input_brush,
                            input_saliva,
                            sample_pair.name,
                            output_dir
                            )
                        ], 
                        name = "run_ichorCNA." + sample_pair.name
                    )
                )

        return jobs

    @property
    def steps(self): # what kind of metrics and reports to add?
        return [
                [
                    self.trimmer,
                    self.bwa_mem_samtools_sort,
                    self.samtools_merge,
                    self.locatit_dedup_bam,
                    self.samtools_sort,
                    self.samtools_index,
                    self.vardict_single,
                ],
                [
                    self.trimmer, # same trimming, mapping, dedup steps for copy number protocol? 
                    self.bwa_mem_samtools_sort,
                    self.samtools_merge,
                    self.locatit_dedup_bam,
                    self.samtools_sort,
                    self.samtools_index,
                    self.hmm_readCounter,
                    self.run_ichorCNA
                 ]
            ]


if __name__ == "__main__":
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        DOvEE_gene(protocol=['vardict', 'copy-number'])
