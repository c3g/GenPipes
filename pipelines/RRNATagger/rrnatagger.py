#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(sys.argv[0]))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.design import *
from bio.readset import *

from bio import rrna_amplicons
from bio import microbial_ecology

from pipelines.illumina import illumina

log = logging.getLogger(__name__)

class RRNATagger(illumina.Illumina):

    def merge_barcodes(self):
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.

    def remove_contam(self):

    def split_barcodes(self):

    def qscores(self):

    def qc(self):

    def qscores_qced(self):

    def generate_clusters(self):

    def classify(self): # Add PacBio's blast here!
    
    def generate_otu_table(self):

    def filter_otu_table(self):

    def summarize_taxonomy(self):

    def plot_taxonomy(self):

    def deliverables(self):

    def cleanup(self):




    @property
    def contrasts(self):
        if not hasattr(self, "_contrasts"):
            self._contrasts = parse_design_file(self.args.design.name, self.samples)
        return self._contrasts

    def trim_metrics(self):
        # Transform pipeline 'PAIRED_END' or 'SINGLE_END' run_type to 'paired' or 'single' parameter
        run_type = re.sub("_END$", "", self.run_type).lower()

        job = metrics.merge_trimmomatic_stats("trim.stats.csv", "trim", os.path.join("metrics", "trimming.stats"), run_type)
        job.input_files = [os.path.join("trim", readset.sample.name, readset.name + ".trim.stats.csv") for readset in self.readsets]
        job.name = "trim_metrics"
        return [job]

    def tophat(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_alignment_directory = os.path.join(alignment_directory, readset.name)

            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job = tophat.tophat(
                fastq1,
                fastq2,
                readset_alignment_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library,
                rg_platform_unit=readset.run + "_" + readset.lane,
                rg_platform=config.param('tophat', 'platform'),
                rg_center=config.param('tophat', 'TBInstitution'),
            )

            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            if len(readset.sample.readsets) == 1:
                sample_bam = os.path.join(alignment_directory, readset.sample.name + ".sorted.bam")
                job.command += " && \\\nln -s " + os.path.join(readset_alignment_directory, "accepted_hits.bam") + " " + sample_bam
                job.output_files.append(sample_bam)

            job.name = "tophat." + readset.name
            jobs.append(job)

        return jobs

    def picard_merge_sam_files(self):
        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name + ".sorted.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def picard_reorder_sam(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)

            job = picard.reorder_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".merged.karyotypic.bam"
            )
            job.name = "picard_reorder_sam." + sample.name
            jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".merged.")

            job = picard.mark_duplicates(
                [alignment_file_prefix + "karyotypic.bam"],
                alignment_file_prefix + "mdup.bam",
                alignment_file_prefix + "mdup.metrics"
            )
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def rnaseqc(self):

        sample_file = os.path.join("alignment", "rnaseqc.samples.txt")
        output_directory = os.path.join("metrics", "rnaseqRep")

        job = metrics.rnaseqc(sample_file, output_directory, self.run_type == "SINGLE_END")

        sample_file_job = Job([os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam") for sample in self.samples], [sample_file])
        project_name = config.param('DEFAULT', 'project_name')
        job.command = """\
mkdir -p {output_directory} && \\
echo \\"Sample\tBamFile\tNote
{input_bams}\\" \
> {sample_file} && \\
{job.command} && \\
zip -r {output_directory}.zip {output_directory}""".format(
            output_directory=output_directory,
            input_bams=" \\\n".join(["\t".join([sample.name, os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam"), project_name]) for sample in self.samples]),
            sample_file=sample_file,
            job=job
        )

        job.input_files.extend([os.path.join("alignment", sample.name, sample.name + ".merged.mdup.bam") for sample in self.samples])
        job.output_files.extend([sample_file, output_directory + ".zip"])
        job.name = "rnaseqc"
        return [job]

    def wiggle(self):
        jobs = []

        for sample in self.samples:
            bam_file_prefix = os.path.join("alignment", sample.name, sample.name + ".merged.mdup.")
            input_bam = bam_file_prefix + "bam"
            bed_graph_prefix = os.path.join("tracks", sample.name, sample.name)
            big_wig_prefix = os.path.join("tracks", "bigWig", sample.name)

            if config.param('tophat', 'strandInfo') != 'fr-unstranded':
                input_bam_f1 = bam_file_prefix + "tmp1.forward.bam"
                input_bam_f2 = bam_file_prefix + "tmp2.forward.bam"
                input_bam_r1 = bam_file_prefix + "tmp1.reverse.bam"
                input_bam_r2 = bam_file_prefix + "tmp2.reverse.bam"
                output_bam_f = bam_file_prefix + "forward.bam"
                output_bam_r = bam_file_prefix + "reverse.bam"

                bam_f_job = concat_jobs([
                    samtools.view(input_bam, input_bam_f1, "-bh -F 256 -f 81"),
                    samtools.view(input_bam, input_bam_f2, "-bh -F 256 -f 161"),
                    picard.merge_sam_files([input_bam_f1, input_bam_f2], output_bam_f),
                    Job(command="rm " + input_bam_f1 + " " + input_bam_f2)
                ])
                bam_f_job.name = "wiggle." + sample.name + ".forward_strandspec"

                bam_r_job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks", "bigWig")),
                    samtools.view(input_bam, input_bam_r1, "-bh -F 256 -f 97"),
                    samtools.view(input_bam, input_bam_r2, "-bh -F 256 -f 145"),
                    picard.merge_sam_files([input_bam_r1, input_bam_r2], output_bam_r),
                    Job(command="rm " + input_bam_r1 + " " + input_bam_r2)
                ])
                bam_r_job.name = "wiggle." + sample.name + ".reverse_strandspec"

                jobs.extend([bam_f_job, bam_r_job])

                outputs = [
                    [bed_graph_prefix + ".forward.bedGraph", big_wig_prefix + ".forward.bw"],
                    [bed_graph_prefix + ".reverse.bedGraph", big_wig_prefix + ".reverse.bw"],
                ]
            else:
                outputs = [[bed_graph_prefix + ".bedGraph", big_wig_prefix + ".bw"]]

            for bed_graph_output, big_wig_output in outputs:
                job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks", "bigWig")),
                    bedtools.graph(input_bam, bed_graph_output, big_wig_output)
                ])
                job.name = "wiggle." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output))
                jobs.append(job)

        return jobs

    def raw_counts(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)
            input_bam = alignment_file_prefix + ".merged.mdup.bam"
            sorted_bam = alignment_file_prefix + ".queryNameSorted.bam"
            sort_order = "queryname"

            # Sort BAM by query
            job = picard.sort_sam(input_bam, sorted_bam, sort_order)
            job.name = "sort_sam.qnsort." + sample.name
            jobs.append(job)

            # Count reads
            output_count = os.path.join("raw_counts", sample.name + ".readcounts.csv")
            stranded = "no" if config.param('align', 'strandInfo') == "fr-unstranded" else "reverse"
            job = concat_jobs([
                Job(command="mkdir -p raw_counts"),
                htseq.htseq_count(
                    sorted_bam,
                    config.param('htseq', 'referenceGtf', type='filepath'),
                    output_count,
                    config.param('htseq', 'options'),
                    stranded
                )
            ])
            job.name = "htseq_count." + sample.name
            jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.trim_metrics,
            self.tophat,
            self.picard_merge_sam_files,
            self.picard_reorder_sam,
            self.picard_mark_duplicates,
            self.rnaseqc,
            self.wiggle,
            self.raw_counts
        ]

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)

        super(RnaSeq, self).__init__()
        
RnaSeq().submit_jobs()

# Steps from old Perl pipeline:

# MERGE
#1 mergeBarcodes

# REMOVE CONTAM
#2 duk_wrapper_contam
#3 duk_wrapper_phix

# SPLIT BY BARCODES
#4 barcodes

# QSCORES
#5 qscore_sheet_R1_raw
#6 qscore_plot_R1_raw

# QC
#7 itags_QC

# QSCORE QCED
#8 qscore_sheet_R1QCed
#9 qscore_plot_R1QCed

# CLUSTER SEQUENCED
#10 clustering_method_3

# CLASSIFY
#11 RDP 

# GENERATE RAW OTU TABLE
#12 add_taxonomy

# FILTER OTU TABLE
#13 filter_obs
#14 split_otu_table
#15 convert_otu_to_biom
#16 convert_otu_to_biom2
#17 convert_otu_to_biom3
#18 rarefy
#19 rarefy
#20 filter_obs_table2X2
#21 filter_obs_table2X2
#22 filter_obs_table1X1
#23 filter_obs_table1X1

# SUMMARIZE TAXONOMY
#24 summarize_taxonomy_absolute_L1
#25 summarize_taxonomy_absolute_L2
#26 summarize_taxonomy_absolute_L3
#27 summarize_taxonomy_absolute_L4
#28 summarize_taxonomy_absolute_L5
#29 summarize_taxonomy_absolute_L6
#30 summarize_taxonomy_absolute_L7
#31 summarize_taxonomy_relative_L1
#32 summarize_taxonomy_relative_L2
#33 summarize_taxonomy_relative_L3
#34 summarize_taxonomy_relative_L4
#35 summarize_taxonomy_relative_L5
#36 summarize_taxonomy_relative_L6
#37 summarize_taxonomy_relative_L7
#38 summarize_taxonomy_absolute_raw_L1
#39 summarize_taxonomy_absolute_raw_L2
#40 summarize_taxonomy_absolute_raw_L3
#41 summarize_taxonomy_absolute_raw_L4
#42 summarize_taxonomy_absolute_raw_L5
#43 summarize_taxonomy_absolute_raw_L6
#44 summarize_taxonomy_absolute_raw_L7
#45 summarize_taxonomy_relative_raw_L1
#46 summarize_taxonomy_relative_raw_L2
#47 summarize_taxonomy_relative_raw_L3
#48 summarize_taxonomy_relative_raw_L4
#49 summarize_taxonomy_relative_raw_L5
#50 summarize_taxonomy_relative_raw_L6
#51 summarize_taxonomy_relative_raw_L7

# PLOT TAXONOMY
#52 plot_taxonomy_absolute
#53 plot_taxonomy_absolute_raw
#54 phylum_barplot_all
#55 phylum_barplot_all_rel
#56 phylum_barplot_bacteria_rarefied

# BLAST (PACBIO ONLY)
#57 blast raw OTUs

# DELIVERABLES
#58 count_report
#59 txtToPdf
#60 merge_pdf
#61 Nozzle deliverables

# CLEANUP
#62 cleanup
