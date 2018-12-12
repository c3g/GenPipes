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
from os.path import basename

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.sequence_dictionary import *

from pipelines import common
from bfx import tools
from bfx import flash
from bfx import qiime
from bfx import vsearch
from bfx import krona

log = logging.getLogger(__name__)

class AmpliconSeq(common.Illumina):
    """
    Amplicon-Seq Pipeline
    ================

    """
  
    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        super(AmpliconSeq, self).__init__(protocol)


    def flash(self):
        """
        Merge paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/).
        """
        jobs = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            merge_directory = os.path.join("merge", readset.sample.name)
            merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")
            merge_file_prefix_log = os.path.join(merge_directory, readset.name + ".log")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            else:
                raise Exception("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!")

            job = flash.flash(
                fastq1,
                fastq2,
                merge_directory,
                merge_file_prefix,
                readset.name,
                merge_file_prefix_log
            )
            job.samples = [readset.sample]

            jobs.append(concat_jobs([
                # FLASh does not create output directory by default
                Job(command="mkdir -p " + merge_directory),
                job
            ], name="flash." + readset.sample.name))

        return jobs

    def merge_flash_stats(self):
        """
        The paired end merge statistics per readset are merged at this step.
        """

        readset_merge_flash_stats = os.path.join("metrics", "mergeReadsetTable.tsv")
        job = concat_jobs([
            Job(command="mkdir -p metrics"),
            Job(command="echo 'Sample\tReadset\tTrim Paired Reads #\tMerged Paired Reads #\tMerged Paired Reads %' > " + readset_merge_flash_stats)
        ])

        for readset in self.readsets:
            flash_log = os.path.join("merge", readset.sample.name, readset.name + ".log")

            job = concat_jobs([
                job,
                Job(
                    command="""\
printf '{sample}\t{readset}\t' \\
  >> {stats}""".format(
                        sample=readset.sample.name,
                        readset=readset.name,
                        stats=readset_merge_flash_stats,
                        ),
                    samples=[readset.sample]
                )
            ])

            # Retrieve merge statistics using re search in python.
            python_command = """\
module load {module_python}; \\
python -c 'import re; \\
import sys; \\
log_file = open("{flash_log}","r"); \\
merge_stat=[]; \\
merge_stat.append([i.split()[3] for i in log_file if re.search("Total pairs",i)][0]); \\
log_file.seek(0); \\
merge_stat.append([i.split()[3] for i in log_file if re.search("Combined pairs",i)][0]); \\
log_file.seek(0); \\
merge_stat.append([i.split()[3] for i in log_file if re.search("Percent combined",i)][0][:-1]); \\
log_file.close(); \\
print "\t".join(merge_stat)'""".format(
                module_python=config.param('DEFAULT', 'module_python'),
                flash_log=flash_log
            )

            job = concat_jobs([
                job,
                Job(
                    [flash_log],
                    [readset_merge_flash_stats],
                    # Create readset merging stats TSV file with paired read count using python.
                    command="""\
{python_command} \\
  >> {readset_merge_flash_stats}""".format(
                        python_command=python_command,
                        readset_merge_flash_stats=readset_merge_flash_stats
                    )
                )
            ])

        sample_merge_flash_stats = os.path.join("metrics", "mergeSampleTable.tsv")
        report_file = os.path.join("report", "Illumina.flash_stats.md")
        return [concat_jobs([
            job,
            Job(
                [readset_merge_flash_stats],
                [sample_merge_flash_stats],
                # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                command="""\
cut -f1,3- {readset_merge_flash_stats} | \\
awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Trim Reads #", "Merged Reads #", "Merged %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_flash_stats}""".format(
                    readset_merge_flash_stats=readset_merge_flash_stats,
                    sample_merge_flash_stats=sample_merge_flash_stats
                )
            ),
            Job(
                [sample_merge_flash_stats],
                [report_file],
                [['flash', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {readset_merge_flash_stats} {sample_merge_flash_stats} report/ && \\
merge_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_flash_stats}` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable min_overlap="{min_overlap}" \\
  --variable max_overlap="{max_overlap}" \\
  --variable read_type="{read_type}" \\
  --variable merge_readset_table="$merge_readset_table_md" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    min_overlap=config.param('flash', 'min_overlap', type='int'),
                    max_overlap=config.param('flash', 'max_overlap', type='int'),
                    read_type="Paired",
                    report_template_dir=self.report_template_dir,
                    readset_merge_flash_stats=readset_merge_flash_stats,
                    sample_merge_flash_stats=sample_merge_flash_stats,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file]
            )], name="merge_flash_stats", samples=self.samples)]

    def catenate(self):

        """
        Catenate all the reads in one file for further analysis.

        This step takes as input files:

        1. Merged FASTQ files from previous step flash.
        """

        jobs = []
        input_files = []
        sample_name = []
        catenate_fasta = os.path.join("catenate", "seqs.fna")

        for readset in self.readsets:
            merge_directory = os.path.join("merge", readset.sample.name)
            merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")

            # Find input readset FASTQs first from previous FLASh job,
            input_files.append(merge_file_prefix)
            sample_name.append(str(readset.sample.name).replace("_", "."))

        if config.param('qiime_catenate', 'map_file'):
            job = qiime.catenate(
                input_files,
                sample_name,
                catenate_fasta
            )
            job.name = "catenate"
            job.samples = self.samples
            jobs.append(job)
        else:
            mapbuild_job = tools.py_ampliconSeq(
                [],
                [],
                'map_build',
                """-s {samples}""".format(
                    samples=','.join(sample_name)
                )
            )
            catenate_job = qiime.catenate(
                input_files,
                sample_name,
                catenate_fasta
            )
            catenate_job.samples = self.samples
            jobs.append(concat_jobs([
                mapbuild_job,
                catenate_job
            ], name="catenate"))

        return jobs

    def uchime(self):
        """
        Reference based chimera detection is performed using [vsearch](https://github.com/torognes/vsearch)

        This step takes as input files:

        1. Catenated FASTA file from previous step catenate.
        """

        jobs = []

        cat_sequence_fasta = os.path.join("catenate", "seqs.fna")

        filter_directory = "catenate_without_chimeras"
        filter_fasta = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")
        filter_log = os.path.join(filter_directory, "seqs_chimeras_filtered.log")

        uchime_job = vsearch.uchime(
            cat_sequence_fasta,
            filter_fasta
        )
        uchime_job.samples = self.samples

        job_log = tools.py_ampliconSeq(
            [filter_fasta],
            [filter_log],
            'catenate_stat',
            """\
  -i {filter_fasta} \\
  -j {filter_log}""".format(
                filter_fasta=filter_fasta,
                filter_log=filter_log)
            )

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + filter_directory),
            uchime_job,
            job_log
        ], name="uchime"))

        return jobs

    def merge_uchime_stats(self):
        """
        The chimeric sequences filtered out statistics per readset are merged at this step.
        """

        readset_merge_uchime_stats = os.path.join("metrics", "uchimeReadsetTable.tsv")
        job = concat_jobs([
            Job(command="mkdir -p metrics"),
            Job(command="echo 'Sample\tReadset\tMerged Paired Reads #\tFiltered Paired Reads #\tFiltered Paired Reads %' > " + readset_merge_uchime_stats)
        ])

        # Database
        if config.param('uchime', 'name') == 'gold':
            chimera_db = 'GOLD'
            chimera_ref = 'gold'
        else:
            chimera_db = 'UNITE'
            chimera_ref = 'unite'

        filter_directory = "catenate_without_chimeras"
        filter_log = os.path.join(filter_directory, "seqs_chimeras_filtered.log")

        for readset in self.readsets:
            flash_log = os.path.join("merge", readset.sample.name, readset.name + ".log")

            job = concat_jobs([
                job,
                Job(
                    command="""\
printf '{sample}\t{readset}\t' \\
  >> {stats}""".format(
                        sample=readset.sample.name,
                        readset=readset.name,
                        stats=readset_merge_uchime_stats
                    ),
                    samples=[readset.sample]
            )])

            job = concat_jobs([
                job,
                tools.py_ampliconSeq(
                    [filter_log, flash_log],
                    [readset_merge_uchime_stats],
                    'uchime',
                    """\
  -i {filter_log} \\
  -j {flash_log} \\
  -s {sample} \\
  >> {stats}""".format(
                        filter_log=filter_log,
                        flash_log=flash_log,
                        sample=str(readset.sample.name),
                        stats=readset_merge_uchime_stats
                    )
                )

            ])

        sample_merge_uchime_stats = os.path.join("metrics", "uchimeSampleTable.tsv")
        report_file = os.path.join("report", "AmpliconSeq.uchime.md")

        return [concat_jobs([
            job,
            Job(
                [],
                [sample_merge_uchime_stats],
                # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                command="""\
cut -f1,3- {readset_merge_uchime_stats} | awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Merged Reads #", "Chimera filtered out Reads #", "Chimera filtered out %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_uchime_stats}""".format(
                    readset_merge_uchime_stats=readset_merge_uchime_stats,
                    sample_merge_uchime_stats=sample_merge_uchime_stats
                )
            ),
            Job(
                [sample_merge_uchime_stats],
                [report_file],
                [['uchime', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {readset_merge_uchime_stats} {sample_merge_uchime_stats} report/ && \\
uchime_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_uchime_stats}` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable read_type="{read_type}" \\
  --variable sequence_max_n="{sequence_max_n}" \\
  --variable chimera_db="{chimera_db}" \\
  --variable chimera_ref="{chimera_ref}" \\
  --variable uchime_readset_table="$uchime_readset_table_md" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    read_type="Paired",
                    report_template_dir=self.report_template_dir,
                    sequence_max_n=config.param('qiime_catenate', 'sequence_max_n'),
                    chimera_db=chimera_db,
                    chimera_ref=chimera_ref,
                    readset_merge_uchime_stats=readset_merge_uchime_stats,
                    sample_merge_uchime_stats=sample_merge_uchime_stats,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file]
            )], name="merge_uchime_stats", samples=self.samples)]

    def otu_picking(self):
        """
        The OTU picking step (de novo & close_ref) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

        This step takes as input file:

        1. Catenated and filtered FASTA file from previous step.

        """

        jobs = []

        method = config.param('qiime_otu_picking', 'method')
        otu_picking_method = "otu_" + method + "_picking"

        filter_directory = "catenate_without_chimeras"
        filter_fastq = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")

        output_directory = method + "_otus/pick_otus"
        otu_file = os.path.join(output_directory, "seqs_chimeras_filtered_otus.txt")

        job = eval("qiime." + otu_picking_method)(
            filter_fastq,
            output_directory,
            otu_file
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + method + "_otus/"),
            job
        ], name="qiime_otu_picking." + method))

        return jobs

    def otu_rep_picking(self):
        """
        After picking OTUs, this step pick a representative sequence for each OTU.

        This step takes as input files:

        1. OTU file from previous step
        2. Catenated and filtered FASTA file from filter_chimeras step.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_file = self.select_input_files([os.path.join(otu_directory, "pick_otus", "seqs_chimeras_filtered_otus.txt")] for otu_directory in otu_directories)[0]

        filter_directory = "catenate_without_chimeras"
        filter_fasta = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")

        otu_directory = os.path.dirname(os.path.dirname(otu_file))
        output_directory = os.path.join(otu_directory, "pick_rep_set")
        otu_rep_file = os.path.join(output_directory, "rep_set.fna")

        job = qiime.otu_rep_picking(
            otu_file,
            filter_fasta,
            otu_rep_file
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + output_directory),
            job
        ], name="qiime_rep_picking." + re.sub("_otus", "", otu_directory)))

        return jobs

    def otu_assigning(self):
        """
        Given a set of OTUS, this step attempts to assign the taxonomy of each OTU using [Uclust] (http://drive5.com/usearch/manual/uclust_algo.html).

        This step takes as input files:

        1. OTU representative sequence file from previous step.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_rep_picking_fasta = self.select_input_files([os.path.join(otu_directory, "pick_rep_set", "rep_set.fna")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(os.path.dirname(otu_rep_picking_fasta))
        output_directory = os.path.join(otu_directory, "taxonomy_assignment")
        tax_assign_file = os.path.join(output_directory, "rep_set_tax_assignments.txt")
 
        job = qiime.otu_assigning(
            otu_rep_picking_fasta,
            output_directory,
            tax_assign_file
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + output_directory),
            job
        ], name="qiime_otu_assigning." + re.sub("_otus", "", otu_directory)))

        return jobs

    def otu_table(self):
        """
        This step make a consensus OTU table in biom format. It tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU.

        This step takes as input files:

        1. OTU picking file.
        2. Taxonomy assignment for each OTU from the previous step.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_file = self.select_input_files([os.path.join(otu_directory, "pick_otus", "seqs_chimeras_filtered_otus.txt")] for otu_directory in otu_directories)[0]
        tax_assign_file = self.select_input_files([os.path.join(otu_directory, "taxonomy_assignment", "rep_set_tax_assignments.txt")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(os.path.dirname(otu_file))
        otu_table_file = os.path.join(otu_directory, "otu_table_prefiltered.biom")
        otu_table_sample_file = os.path.join(otu_directory, "otu_table_prefiltered_sample.biom")
        otu_table_final = os.path.join(otu_directory, "otu_table.biom")
        otu_table_summary = os.path.join(otu_directory, "otu_table.sum")

        otu_sample_directory = os.path.join(otu_directory, "sample")
        sample_name_control = os.path.join(otu_sample_directory, "done.txt")

        job = qiime.otu_table(
            otu_file,
            tax_assign_file,
            otu_table_file,
            otu_table_summary
        )
        job.samples = self.samples

        # Remove singleton
        job_filter = Job(
            [otu_table_file],
            [otu_table_sample_file],
            command="""\
$QIIME_HOME/filter_otus_from_otu_table.py \\
  -i {otu_table_file} \\
  -n 2 \\
  -o {otu_table_sample_file}""".format(
                otu_table_file=otu_table_file,
                otu_table_sample_file=otu_table_sample_file
            )
        )

        # Remove samples with fewer than 2 observations.
        job_filter2 = Job(
            [otu_table_sample_file],
            [otu_table_final],
            command="""\
$QIIME_HOME/filter_samples_from_otu_table.py \\
  -i {otu_table_sample_file} \\
  -n 2 \\
  -o {otu_table_final}""".format(
                otu_table_sample_file=otu_table_sample_file,
                otu_table_final=otu_table_final
            )
        )

        # Sample remained after filtering.
        job_sample = tools.py_ampliconSeq(
            [otu_table_summary],
            [sample_name_control],
            'sample_name',
        """\
  -i {otu_table_summary} \\
  -j {otu_sample_directory}""".format(
                otu_table_summary=otu_table_summary,
                otu_sample_directory=otu_sample_directory
            )
        )

        jobs.append(concat_jobs([
            job,
            job_filter,
            job_filter2,
            Job(command="""\
$QIIME_HOME/biom summarize-table \\
  -i {otu_table_final} \\
  > {otu_table_summary}""".format(
                otu_table_final=otu_table_final,
                otu_table_summary=otu_table_summary
            )),
            Job(command="mkdir -p " + otu_sample_directory),
            job_sample
        ], name="qiime_otu_table." + re.sub("_otus", "", otu_directory)))

        return jobs

    def otu_alignment(self):
        """
        Align OTU representative sequences using [PyNAST] (http://biocore.github.io/pynast/).

        This step takes as input file:

        1. OTU representative sequence file.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_rep_picking_fasta = self.select_input_files([os.path.join(otu_directory, "pick_rep_set", "rep_set.fna")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(os.path.dirname(otu_rep_picking_fasta))
        output_directory = os.path.join(otu_directory, "align_seq")
        align_seq_fasta= os.path.join(output_directory, "rep_set_aligned.fasta")

        job = qiime.otu_alignment(
            otu_rep_picking_fasta,
            output_directory,
            align_seq_fasta
        )
        job.samples = self.samples

        job.name = "qiime_otu_alignment." + re.sub("_otus", "", otu_directory)
        jobs.append(job)

        return jobs

    def filter_alignment(self):
        """
        Filter the alignment by removing positions which are gaps in every sequence.

        This step takes as input file:

        1. Alignment sequence file.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        align_seq_fasta = self.select_input_files([os.path.join(otu_directory, "align_seq", "rep_set_aligned.fasta")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(os.path.dirname(align_seq_fasta))
        output_directory = os.path.join(otu_directory, "filter_alignment")
        filter_align_fasta = os.path.join(output_directory, "rep_set_aligned_pfiltered.fasta")

        job = qiime.filter_alignment(
            align_seq_fasta,
            output_directory,
            filter_align_fasta
        )
        job.samples = self.samples

        job.name = "qiime_filter_alignment." + re.sub("_otus", "", otu_directory)
        jobs.append(job)

        return jobs

    def phylogeny(self):
        """
        Build a phylogenetic tree from a multiple sequence alignment using [FastTree] (http://www.microbesonline.org/fasttree/).

        This step takes as input file:

        1. Filtered alignment sequence file from previous step.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        filter_align_fasta = self.select_input_files([os.path.join(otu_directory, "filter_alignment", "rep_set_aligned_pfiltered.fasta")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(os.path.dirname(filter_align_fasta))
        output_directory = os.path.join(otu_directory, "phylogenetic_tree")
        phylo_file = os.path.join(output_directory, "rep_phylo.tre")

        job = qiime.phylogeny(
            filter_align_fasta,
            phylo_file
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + output_directory),
            job
        ], name="qiime_phylogeny." + re.sub("_otus", "", otu_directory)))

        return jobs

    def qiime_report(self):
        """
        1st part report for taxonomic affiliation.
        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_table = self.select_input_files([os.path.join(otu_directory, "otu_table.biom")] for otu_directory in otu_directories)[0]
        phylo_file = self.select_input_files([os.path.join(otu_directory, "phylogenetic_tree", "rep_phylo.tre")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(otu_table)

        report_file = os.path.join("report", "AmpliconSeq.qiime.md")

        if config.param('qiime', 'amplicon_type') == '16S':
            amp_db = 'Greengenes'
        elif config.param('qiime', 'amplicon_type') == '18S':
            amp_db = 'Silva'
        elif config.param('qiime', 'amplicon_type') == 'ITS':
            amp_db = 'UNITE'
        else:
            amp_db = 'Unknown'

        jobs.append(Job(
            [otu_table, phylo_file],
            [report_file],
            [['qiime', 'module_pandoc']],
            command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable amplicon_type="{amplicon_type}" \\
  --variable similarity="{similarity}" \\
  --variable amplicon_db="{amplicon_db}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                amplicon_type=config.param('qiime', 'amplicon_type'),
                similarity=config.param('qiime_otu_picking', 'similarity'),
                amplicon_db=amp_db,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            ),
            report_files=[report_file],
            name="qiime_report." + re.sub("_otus", "", otu_directory),
            samples=self.samples
        ))

        return jobs

    def multiple_rarefaction(self):
        """
        1st step (/4) for rarefaction plot.
        Rarefies OTU table by random sampling (without replacement) at different depth in order to perform rarefaction analysis.
        You need to provide the minimum/maximum number of sequences per samples and the size of each steps between the min/max of seqs/sample.

        This step takes as input files:

        1. OTU non rarefied table in biom format.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_table = self.select_input_files([os.path.join(otu_directory, "otu_table.biom")] for otu_directory in otu_directories)[0]

        otus_input = [otu_table]

        otu_directory = os.path.dirname(otu_table)
        alpha_directory = re.sub('otus', 'alpha_diversity', otu_directory)
        rarefied_otu_directory = os.path.join(alpha_directory, "rarefied_otu_tables")

        job = qiime.multiple_rarefaction(
            otus_input,
            rarefied_otu_directory
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + alpha_directory),
            job
        ], name="qiime_multiple_rarefaction." + re.sub("_otus", "", otu_directory)))

        return jobs

    def alpha_diversity(self):
        """
        2nd step (/4) for rarefaction plot.
        Calculate alpha diversity on each sample using a variety of alpha diversity metrics (chao1, shannon, observed otus).

        This step takes as input files:

        1. Multiple OTU rarefied table in biom format from previous step.

        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']
        rarefied_otu_directory = self.select_input_files([os.path.join(alpha_directory, "rarefied_otu_tables")] for alpha_directory in alpha_directories)[0]

        alpha_directory = os.path.dirname(rarefied_otu_directory)
        alpha_diversity_directory = os.path.join(alpha_directory, "alpha_diversity_compute")

        job = qiime.alpha_diversity(
            rarefied_otu_directory,
            alpha_diversity_directory
        )
        job.samples = self.samples
        job.name = "qiime_alpha_diversity." + re.sub("_alpha_diversity", "", alpha_directory)

        jobs.append(job)

        return jobs

    def collate_alpha(self):
        """
        3rd step (/4) for rarefaction plot.
        Merge all the alpha diversity computed in the previous step.
        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']
        alpha_diversity_directory = self.select_input_files([os.path.join(alpha_directory, "alpha_diversity_compute")] for alpha_directory in alpha_directories)[0]

        alpha_directory = os.path.dirname(alpha_diversity_directory)
        alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
        alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
        chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
        observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
        shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")

        job = qiime.collate_alpha(
            alpha_diversity_directory,
            alpha_diversity_collated_merge_directory,
            chao1_stat,
            observed_species_stat,
            shannon_stat
        )
        job.samples = self.samples

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + alpha_diversity_collated_directory),
            job
        ], name="qiime_collate_alpha." + re.sub("_alpha_diversity", "", alpha_directory)))

        return jobs

    def sample_rarefaction_plot(self):
        """
        Last step for rarefaction plot.
        Plot the rarefaction curve for each sample
        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']
        chao1_stat = self.select_input_files([os.path.join(alpha_directory, "alpha_diversity_collated", "merge_samples", "chao1.txt")] for alpha_directory in alpha_directories)[0]

        alpha_directory = os.path.dirname(os.path.dirname(os.path.dirname(chao1_stat)))
        alpha_diversity_collated_directory = os.path.dirname(os.path.dirname(chao1_stat))

        curve_sample=[]
        alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
        alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")

        observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
        shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")

        for readset in self.readsets:

            sample_collated_general_directory = os.path.join(alpha_diversity_collated_directory, readset.sample.name)
            sample_map = os.path.join(sample_collated_general_directory, "map.txt")
            sample_collated_directory = os.path.join(sample_collated_general_directory, "stat")
            sample_rarefaction_directory = os.path.join(alpha_diversity_rarefaction_directory, readset.sample.name)
            chao1_dir = os.path.join(sample_collated_directory, "chao1.txt")
            observed_species_dir = os.path.join(sample_collated_directory, "observed_species.txt")
            shannon_dir = os.path.join(sample_collated_directory, "shannon.txt")
            observed_species_file = """{}/average_plots/observed_species{}.png""".format(readset.sample.name,str(readset.sample.name).replace('_','.'))
            curve_sample.append(os.path.join(alpha_diversity_rarefaction_directory,observed_species_file))

            job = qiime.sample_rarefaction_plot(
                chao1_stat,
                observed_species_stat,
                shannon_stat,
                sample_collated_directory,
                sample_map,
                sample_rarefaction_directory,
                curve_sample,
            )
            job.samples = [readset.sample]

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + sample_collated_directory),
                tools.py_ampliconSeq(
                    [],
                    [sample_map],
                    'map_per_sample',
                    """\
  -s {sample} \\
  -j {sample_map}""".format(
                        sample=readset.sample.name,
                        sample_map=sample_map
                    )
                ),
                tools.py_ampliconSeq(
                    [chao1_stat],
                    [chao1_dir],
                    'sample_rarefaction',
                    """\
  -i {cstat} \\
  -j {cdir} \\
  -s {readset}""".format(
                        cstat=chao1_stat,
                        cdir=chao1_dir,
                        readset=readset.sample.name
                    )
                ),
                tools.py_ampliconSeq(
                    [observed_species_stat],
                    [observed_species_dir],
                    'sample_rarefaction',
                    """\
  -i {ostat} \\
  -j {odir} \\
  -s {readset}""".format(
                        ostat=observed_species_stat,
                        odir=observed_species_dir,
                        readset=readset.sample.name
                    )
                ),
                tools.py_ampliconSeq(
                    [shannon_stat],
                    [shannon_dir],
                    'sample_rarefaction',
                    """\
  -i {sstat} \\
  -j {sdir} \\
  -s {readset}""".format(
                        sstat=shannon_stat,
                        sdir=shannon_dir,
                        readset=readset.sample.name
                    )
                ),
                job
            ], name="qiime_sample_rarefaction." + re.sub("_alpha_diversity", ".", alpha_directory) + readset.sample.name))

        return jobs

    def qiime_report2(self):
        """
        2nd part report for taxonomic affiliation. Plot rarefaction curve for each sample.
        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']

        report_file = os.path.join("report", "AmpliconSeq.plot_curve_no_rar.md")
        num_sample = 0

        inputs = []
        curve_sample = []
        for readset in self.readsets:
            inputs.append(str(readset.sample.name).replace('_','.'))
            observed_species_file="""\
{sample}/average_plots/observed_species{readset}.png""".format(
                sample=readset.sample.name,
                readset=str(readset.sample.name).replace('_','.')
            )
            curve_sample.append(self.select_input_files([os.path.join(alpha_directory, "alpha_rarefaction", observed_species_file)] for alpha_directory in alpha_directories)[0])
            num_sample+=1

        alpha_directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(curve_sample[0]))))

        fig_src_link = "fig/" + alpha_directory + "/alpha.rarefaction_sample.png"

        jobs.append(Job(
            curve_sample,
            [report_file],
            [['qiime', 'module_pandoc']],
            command="""\
mkdir -p report/fig/{alpha_directory}/ && \\
montage -mode concatenate -tile {plot_dimension} {curve_sample} report/fig/{alpha_directory}/alpha.rarefaction_sample.png  && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable fig_src_link="{fig_src_link}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                alpha_directory=alpha_directory,
                plot_dimension=str(math.sqrt(num_sample)+1)+"x",  # "3x"+str((num_sample/3)+1)
                curve_sample=' '.join(curve_sample),
                fig_src_link=fig_src_link,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            ),
            report_files=[report_file],
            samples=self.samples,
            name="qiime_report2." + re.sub("_alpha_diversity", "", alpha_directory))
        )

        return jobs

    def single_rarefaction(self):
        """
        This step is recommended. It subsamples (rarefy) all the samples to an equal number of sequences for further comparaison.
        You have to provide the number of sequences to subsample per sample in the configuration file (single_rarefaction_depth).

        This step takes as input files:

        1. OTU table in biom format.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_table = self.select_input_files([os.path.join(otu_directory, "otu_table.biom")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(otu_table)

        otu_normalized_directory = re.sub('otus', 'otu_rarefaction_normalized', otu_directory)
        otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
        normalization_method = os.path.join(otu_normalized_directory,"rarefaction.txt")

        alpha_directory = re.sub('otus', 'alpha_diversity', otu_directory)
        alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
        alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
        chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
        observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
        shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")

        alpha_diversity_collated_merge_rarefied_directory = os.path.join(alpha_diversity_collated_directory, "rarefaction/merge_samples_rarefied")
        chao1_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "chao1.txt")
        observed_species_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "observed_species.txt")
        shannon_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "shannon.txt")

        job = qiime.single_rarefaction(
            otu_table,
            chao1_rarefied_stat,
            observed_species_rarefied_stat,
            shannon_rarefied_stat,
            otu_normalized_table,
            normalization_method
        )
        job.samples = self.samples

        job_chao1 = tools.py_ampliconSeq(
            [chao1_stat],
            [chao1_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {depth}""".format(
                stat=chao1_stat,
                rarefied_stat=chao1_rarefied_stat,
                depth=config.param('qiime_single_rarefaction', 'single_rarefaction_depth')
            )
        )

        job_observed_species = tools.py_ampliconSeq(
            [observed_species_stat],
            [observed_species_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {depth}""".format(
                stat=observed_species_stat,
                rarefied_stat=observed_species_rarefied_stat,
                depth=config.param('qiime_single_rarefaction', 'single_rarefaction_depth')
            )
        )

        job_shannon = tools.py_ampliconSeq(
            [shannon_stat],
            [shannon_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {depth}""".format(
                stat=shannon_stat,
                rarefied_stat=shannon_rarefied_stat,
                depth=config.param('qiime_single_rarefaction', 'single_rarefaction_depth')
            )
        )

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + otu_normalized_directory),
            Job(command="mkdir -p " + alpha_diversity_collated_merge_rarefied_directory),
            Job(command="touch " + normalization_method),
            job_chao1,
            job_observed_species,
            job_shannon,
            job
        ], name="qiime_single_rarefaction." + re.sub("_otus", "", otu_directory)))

        return jobs

    def css_normalization(self):
        """
        This step is recommended. Alternative method for normalization to rarefaction.
        Performs the CSS Matrix normalization.

        This step takes as input files:

        1. OTU table in biom format.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']
        otu_table = self.select_input_files([os.path.join(otu_directory, "otu_table.biom")] for otu_directory in otu_directories)[0]

        otu_directory = os.path.dirname(otu_table)

        otu_normalized_directory = re.sub('otus', 'otu_css_normalized', otu_directory)
        otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
        normalization_method = os.path.join(otu_normalized_directory,"css.txt")

        alpha_directory = re.sub('otus', 'alpha_diversity', otu_directory)
        alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
        alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
        chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
        observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
        shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")

        alpha_diversity_collated_merge_rarefied_directory = os.path.join(alpha_diversity_collated_directory, "css/merge_samples_rarefied")
        chao1_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "chao1.txt")
        observed_species_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "observed_species.txt")
        shannon_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "shannon.txt")

        job = qiime.css_normalization(
            otu_table,
            chao1_rarefied_stat,
            observed_species_rarefied_stat,
            shannon_rarefied_stat,
            otu_normalized_table,
            normalization_method
        )
        job.samples = self.samples

        job_chao1 = tools.py_ampliconSeq(
            [chao1_stat],
            [chao1_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {raremax}""".format(
                stat=chao1_stat,
                rarefied_stat=chao1_rarefied_stat,
                raremax=config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max')
            )
        )

        job_observed_species = tools.py_ampliconSeq(
            [observed_species_stat],
            [observed_species_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {raremax}""".format(
                stat=observed_species_stat,
                rarefied_stat=observed_species_rarefied_stat,
                raremax=config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max')
            )
        )

        job_shannon = tools.py_ampliconSeq(
            [shannon_stat],
            [shannon_rarefied_stat],
            'single_rarefaction',
            """\
  -i {stat} \\
  -j {rarefied_stat} \\
  -s {raremax}""".format(
                stat=shannon_stat,
                rarefied_stat=shannon_rarefied_stat,
                raremax=config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max')
            )
        )

        jobs.append(concat_jobs([
            # Create an output directory
            Job(command="mkdir -p " + otu_normalized_directory),
            Job(command="mkdir -p " + alpha_diversity_collated_merge_rarefied_directory),
            Job(command="touch " + normalization_method),
            job_chao1,
            job_observed_species,
            job_shannon,
            job
        ], name="qiime_css_normalization." + re.sub("_otus", "", otu_directory)))

        return jobs

    def rarefaction_plot(self):
        """
        Last step for rarefaction plot.
        Rarefaction curve for each sample on the same plot.
        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']

        for method in ['css', 'rarefaction']:

            chao1_stat = self.select_input_files([os.path.join(alpha_directory, "alpha_diversity_collated", method, "merge_samples_rarefied", "chao1.txt")] for alpha_directory in alpha_directories)[0]
            observed_species_stat = self.select_input_files([os.path.join(alpha_directory, "alpha_diversity_collated", method, "merge_samples_rarefied", "observed_species.txt")] for alpha_directory in alpha_directories)[0]
            shannon_stat = self.select_input_files([os.path.join(alpha_directory, "alpha_diversity_collated", method, "merge_samples_rarefied", "shannon.txt")] for alpha_directory in alpha_directories)[0]

            alpha_directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(chao1_stat))))
            alpha_diversity_collated_merge_rarefied_directory = os.path.dirname(chao1_stat)

            alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
            alpha_diversity_rarefaction_rarefied_directory = os.path.join(alpha_diversity_rarefaction_directory, method, "merge_samples_rarefied")
            alpha_diversity_rarefaction_file = os.path.join(alpha_diversity_rarefaction_rarefied_directory, "rarefaction_plots.html")

            if config.param('rarefaction_plot', 'map_file'):
                map_file = config.param('rarefaction_plot', 'map_file')
            else:
                map_file = "map.txt"

            job = qiime.rarefaction_plot(
                alpha_diversity_collated_merge_rarefied_directory,
                chao1_stat,
                observed_species_stat,
                shannon_stat,
                map_file,
                alpha_diversity_rarefaction_file,
                alpha_diversity_rarefaction_rarefied_directory
            )
            job.samples = self.samples

            job.name = "qiime_rarefaction_plot." + re.sub("_alpha_diversity", ".", alpha_directory[0] + method)
            jobs.append(job)

        return jobs

    def summarize_taxa(self):
        """
        1st step (/3) for taxonomic affiliation plot.
        Summarize information of taxonomic groups within each sample at different taxonomic level.

        This step takes as input files:

        1. OTU rarefied table in biom format if available.
        2. Else, OTU non rarefied table in biom format.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']

        for method in ['css', 'rarefaction']:

            otu_normalized_table = self.select_input_files([os.path.join(re.sub('otus', 'otu_' + method + '_normalized', otu_directory),"otu_normalized_table.biom")] for otu_directory in otu_directories)[0]

            alpha_directory = re.sub('otu_' + method + '_normalized', 'alpha_diversity', os.path.dirname(otu_normalized_table))
            taxonomic_directory = os.path.join(alpha_directory, method, "taxonomic_affiliation/")
            taxonomic_phylum = os.path.join(taxonomic_directory, "otu_normalized_table_L2.txt")
            taxonomic_class = os.path.join(taxonomic_directory, "otu_normalized_table_L3.txt")
            taxonomic_order = os.path.join(taxonomic_directory, "otu_normalized_table_L4.txt")
            taxonomic_family = os.path.join(taxonomic_directory, "otu_normalized_table_L5.txt")
            taxonomic_genus = os.path.join(taxonomic_directory, "otu_normalized_table_L6.txt")

            job = qiime.summarize_taxa(
                otu_normalized_table,
                taxonomic_directory,
                taxonomic_phylum,
                taxonomic_class,
                taxonomic_order,
                taxonomic_family,
                taxonomic_genus
            )
            job.samples = self.samples

            jobs.append(concat_jobs([
                # Create an output directory
                Job(command="mkdir -p " + alpha_directory),
                job
            ], name="qiime_summarize_taxa." + re.sub("_alpha_diversity", ".", alpha_directory) + method))

        return jobs

    def plot_taxa(self):
        """
        2nd step (/3) for taxonomic affiliation plot.
        Make taxonomy summary bar plots based on taxonomy assignment.

        This step takes as input files:

        1. Summarized information from previous step.

        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']

        for method in ['css', 'rarefaction']:

            taxonomic_phylum = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "otu_normalized_table_L2.txt")] for alpha_directory in alpha_directories)[0]
            taxonomic_class = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "otu_normalized_table_L3.txt")] for alpha_directory in alpha_directories)[0]
            taxonomic_order = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "otu_normalized_table_L4.txt")] for alpha_directory in alpha_directories)[0]
            taxonomic_family = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "otu_normalized_table_L5.txt")] for alpha_directory in alpha_directories)[0]
            taxonomic_genus = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "otu_normalized_table_L6.txt")] for alpha_directory in alpha_directories)[0]

            taxonomic_directory = os.path.dirname(taxonomic_phylum)
            alpha_directory = os.path.dirname(os.path.dirname(taxonomic_directory))

            taxonomic_input = [taxonomic_phylum, taxonomic_class, taxonomic_order, taxonomic_family, taxonomic_genus]

            alpha_diversity_taxonomy_bar_plot = os.path.join(taxonomic_directory, "bar_charts.html")

            job = qiime.plot_taxa(
                taxonomic_input,
                alpha_diversity_taxonomy_bar_plot,
                taxonomic_directory
            )
            job.samples = self.samples

            job.name = "qiime_plot_taxa." + re.sub("_alpha_diversity", ".", alpha_directory) + method
            jobs.append(job)

        return jobs

    def plot_heatmap(self):
        """
        Last step for taxonomic affiliation plot.
        Make heatmap at phylum level.

        This step takes as input files:

        1. Summarized information from previous step.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']

        for method in ['css', 'rarefaction']:

            otu_normalized_table = self.select_input_files([os.path.join(re.sub('otus', 'otu_' + method + '_normalized', otu_directory), "otu_normalized_table.biom")] for otu_directory in otu_directories)[0]
            taxonomic_phylum = self.select_input_files([os.path.join(re.sub('otus', 'alpha_diversity', otu_directory), method, "taxonomic_affiliation", "otu_normalized_table_L2.txt")] for otu_directory in otu_directories)[0]

            otu_directory = re.sub('otu_' + method + '_normalized', 'otus', os.path.dirname(otu_normalized_table))

            beta_directory = re.sub('otus', 'beta_diversity', otu_directory)
            heatmap_directory = os.path.join(beta_directory, method, "heatmap")

            heatmap_script = os.path.join(heatmap_directory, "OTU_Phylum_to_R.R")
            heatmap_chart = os.path.join(heatmap_directory, "otu_heatmap.png")

            heatmap_otu_data_R = os.path.join(heatmap_directory, "OTU_data.txt")
            heatmap_otu_name_R = os.path.join(heatmap_directory, "OTU_name.txt")
            heatmap_otu_tax_R = os.path.join(heatmap_directory, "OTU_tax_final.txt")

            heatmap_otu_table = os.path.join(heatmap_directory, "otumat.tsv")
            heatmap_tax_table = os.path.join(heatmap_directory, "taxmat.tsv")

            job = tools.py_ampliconSeq(
                [taxonomic_phylum],
                [heatmap_script, heatmap_otu_data_R, heatmap_otu_name_R, heatmap_otu_tax_R],
                'plot_heatmap',
                """\
  -i {taxonomic_phylum} \\
  -j {heatmap_directory} \\
  -s {taxon_lvl}""".format(
                    taxonomic_phylum=taxonomic_phylum,
                    heatmap_directory=heatmap_directory,
                    taxon_lvl=1
                )
            )
            job.samples = self.samples

            # Create a job that cleans the generated OTU_data.txt i.e. removes the lines with characters
            jobClean = tools.clean_otu(heatmap_otu_data_R)

            jobR = Job(
                [heatmap_script, heatmap_otu_data_R, heatmap_otu_name_R, heatmap_otu_tax_R],
                [heatmap_chart, heatmap_otu_table, heatmap_tax_table],
                module_entries=[['plot_heatmap', 'module_R']],
                command="./" +  heatmap_script
            )

            jobs.append(concat_jobs([
                # Create an output directory
                Job(command="mkdir -p " + heatmap_directory),
                job,
                Job(command="chmod +x " + heatmap_directory + "/OTU_Phylum_to_R.R"),
                jobClean,
                jobR
            ], name="plot_heatmap." + re.sub("_otus", ".", otu_directory) + method))

        return jobs

    def krona(self):
        """
        Plot Krona chart for taxonomic affiliation
        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']

        for method in ['css', 'rarefaction']:

            sample_name = []

            otu_normalized_table = self.select_input_files([os.path.join(re.sub('otus', 'otu_' + method + '_normalized', otu_directory), "otu_normalized_table.biom")] for otu_directory in otu_directories)[0]
            otu_directory = re.sub('otu_' + method + '_normalized', 'otus', os.path.dirname(otu_normalized_table))

            alpha_directory = re.sub('otus', 'alpha_diversity', otu_directory)
            alpha_diversity_krona_directory = os.path.join(alpha_directory, method, "krona_chart")
            alpha_diversity_krona_file = os.path.join(alpha_diversity_krona_directory, "krona_chart.html")

            for readset in self.readsets:
                sample_name.append(alpha_diversity_krona_directory+'/'+str(readset.sample.name).replace("_", ".")+'.txt')

            job = krona.krona(
                otu_normalized_table,
                sample_name,
                alpha_diversity_krona_file,
            )
            job.samples = self.samples

            jobs.append(concat_jobs([
                # Create an output directory
                Job(command="mkdir -p " + alpha_diversity_krona_directory),
                Job(command="""\
$QIIME_HOME/biom convert -i {otu_normalized_table} \\
  -o {alpha_directory}/{method}/table_tax.txt \\
  --table-type="OTU table" \\
  --to-tsv \\
  --header-key taxonomy""".format(
                    method=method,
                    alpha_directory=alpha_directory,
                    otu_normalized_table=otu_normalized_table
                )),
                tools.py_ampliconSeq(
                    [otu_normalized_table],
                    [],
                    'krona',
                    """\
  -i {alpha_directory}/{method}/table_tax.txt \\
  -j {alpha_diversity_krona_directory}""".format(
                    method=method,
                    alpha_directory=alpha_directory,
                    alpha_diversity_krona_directory=alpha_diversity_krona_directory
                )),
                job
            ], name="krona." + re.sub("_otus", ".", otu_directory) + method))
        return jobs

    def plot_to_alpha(self):
        """
        Final report 1st part for the Amplicon-Seq pipeline. Display results (taxonomy, heatmap and alpha diversity).
        """

        jobs = []

        alpha_directories = ['open_ref_alpha_diversity', 'denovo_alpha_diversity', 'closed_ref_alpha_diversity']

        for method in ['css', 'rarefaction']:

            alpha_diversity_taxonomy_bar_plot = self.select_input_files([os.path.join(alpha_directory, method, "taxonomic_affiliation", "bar_charts.html")] for alpha_directory in alpha_directories)[0]
            alpha_diversity_krona_file = self.select_input_files([os.path.join(alpha_directory, method, "krona_chart", "krona_chart.html")] for alpha_directory in alpha_directories)[0]
            alpha_diversity_rarefaction_file = self.select_input_files([os.path.join(alpha_directory, "alpha_rarefaction", method, "merge_samples_rarefied/rarefaction_plots.html")] for alpha_directory in alpha_directories)[0]

            beta_diversity_heatmap_plot = self.select_input_files([os.path.join(re.sub('alpha_diversity', 'beta_diversity', alpha_directory), method, "heatmap", "otu_heatmap.png")] for alpha_directory in alpha_directories)[0]
            beta_diversity_heatmap_otumat = self.select_input_files([os.path.join(re.sub('alpha_diversity', 'beta_diversity', alpha_directory), method, "heatmap", "otumat.tsv")] for alpha_directory in alpha_directories)[0]
            beta_diversity_heatmap_taxmat = self.select_input_files([os.path.join(re.sub('alpha_diversity', 'beta_diversity', alpha_directory), method, "heatmap", "taxmat.tsv")] for alpha_directory in alpha_directories)[0]

            normalization_method = self.select_input_files([os.path.join(re.sub('alpha_diversity', 'otu_' + method + '_normalized', alpha_directory), method + ".txt")] for alpha_directory in alpha_directories)[0]

            inputs = [alpha_diversity_taxonomy_bar_plot, alpha_diversity_krona_file, alpha_diversity_rarefaction_file, beta_diversity_heatmap_plot, beta_diversity_heatmap_otumat, beta_diversity_heatmap_taxmat, normalization_method]

            alpha_directory = re.sub('otu_' + method + '_normalized', 'alpha_diversity', os.path.dirname(normalization_method))
            beta_directory = re.sub('alpha_diversity', 'beta_diversity', alpha_directory)
            report_file = os.path.join("report", "AmpliconSeq.plot_to_alpha_" + method[:3] + ".md")

            if method == 'rarefaction':
                method_title = "from rarefaction"
                method_link = "These results have been generated after a rarefaction step. All the samples have been rarefied to **" + config.param('qiime_single_rarefaction', 'single_rarefaction_depth') + "** sequences."
            else:
                method_title = "from CSS normalization"
                method_link = "These results have been generated after the [CSS]\ [@css] normalization method."

            bar_plot_link = "[Interactive html plots available here](fig/" + alpha_directory + "/" + method + "/taxonomic_affiliation/bar_charts.html)"
            krona_chart_link = "[available here](fig/" + alpha_directory + "/" + method + "/krona_chart/krona_chart.html)"
            heatmap_otu_link = "[download OTU table](fig/" + beta_directory + "/" + method + "/heatmap/otumat.tsv)"
            heatmap_taxon_link = "[download taxon table](fig/" + beta_directory + "/" + method + "/heatmap/taxmat.tsv))](fig/" + beta_directory + "/" + method + "/heatmap/otu_heatmap.png"
            alpha_plots_link = "[Interactive html plots for alpha diversity available here](fig/" + alpha_directory + "/alpha_rarefaction/" + method + "/rarefaction_plots.html)"

            jobs.append(Job(
                inputs,
                [report_file],
                [['plot_to_alpha', 'module_pandoc']],
                command="""\
mkdir -p report/fig/{alpha_directory}/{method}/ && \\
mkdir -p report/fig/{alpha_directory}/alpha_rarefaction/ && \\
mkdir -p report/fig/{beta_directory}/{method}/heatmap/ && \\
cp -r {alpha_directory}/{method}/taxonomic_affiliation/ report/fig/{alpha_directory}/{method}/taxonomic_affiliation/ && \\
cp -r {alpha_directory}/{method}/krona_chart/ report/fig/{alpha_directory}/{method}/krona_chart/ && \\
cp -r {alpha_directory}/alpha_rarefaction/{method}/merge_samples_rarefied/ report/fig/{alpha_directory}/alpha_rarefaction/{method}/ && \\
cp {beta_directory}/{method}/heatmap/otu_heatmap.png report/fig/{beta_directory}/{method}/heatmap/otu_heatmap.png && \\
cp {beta_directory}/{method}/heatmap/otumat.tsv report/fig/{beta_directory}/{method}/heatmap/otumat.tsv && \\
cp {beta_directory}/{method}/heatmap/taxmat.tsv report/fig/{beta_directory}/{method}/heatmap/taxmat.tsv && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable method_title="{method_title}" \\
  --variable method_link="{method_link}" \\
  --variable bar_plot_link="{bar_plot_link}" \\
  --variable krona_chart_link="{krona_chart_link}" \\
  --variable heatmap_otu_link="{heatmap_otu_link}" \\
  --variable heatmap_taxon_link="{heatmap_taxon_link}" \\
  --variable alpha_plots_link="{alpha_plots_link}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    alpha_directory=alpha_directory,
                    beta_directory=beta_directory,
                    method=method,
                    method_title=method_title,
                    method_link=method_link,
                    bar_plot_link=bar_plot_link,
                    krona_chart_link=krona_chart_link,
                    heatmap_otu_link=heatmap_otu_link,
                    heatmap_taxon_link=heatmap_taxon_link,
                    alpha_plots_link=alpha_plots_link,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=re.sub('_'+method[:3], '', os.path.basename(report_file)),
                    report_file=report_file
                ),
                samples=self.samples,
                name="plot_to_alpha." + re.sub("_alpha_diversity", ".", alpha_directory) + method
            ))

        return jobs

    def beta_diversity(self):
        """
        1st step (/3) for 2D PCoA plot.
        Calculate beta diversity (pairwise sample dissimilarity) on OTU table. The OTU table has to be normalized.
        Only works with >= 4 samples

        This step takes as input files:

        1. OTU rarefied table in biom format.
        2. Tree file.

        """

        jobs = []

        otu_directories = ['open_ref_otus', 'denovo_otus', 'closed_ref_otus']

        for method in ['css', 'rarefaction']:

            otu_normalized_table = self.select_input_files([os.path.join(re.sub('otus', 'otu_' + method + '_normalized', otu_directory), "otu_normalized_table.biom")] for otu_directory in otu_directories)[0]

            otu_directory = re.sub('otu_' + method + '_normalized', 'otus', os.path.dirname(otu_normalized_table))

            phylogenetic_tree_directory = os.path.join(otu_directory, "phylogenetic_tree")
            phylogenetic_tree_file = os.path.join(phylogenetic_tree_directory, "rep_phylo.tre")

            beta_directory = re.sub('otus', 'beta_diversity', otu_directory)
            dm_directory = os.path.join(beta_directory, method, "dissimilarity_matrix")
            dm_unweighted_file = os.path.join(dm_directory, "unweighted_unifrac_otu_normalized_table.txt")
            dm_weighted_file = os.path.join(dm_directory, "weighted_unifrac_otu_normalized_table.txt")
            dm_euclidean_file = os.path.join(dm_directory, "euclidean_otu_normalized_table.txt")

            job = qiime.beta_diversity(
                config.param('beta_diversity', 'beta_diversity_metric'),
                otu_normalized_table,
                phylogenetic_tree_file,
                dm_directory,
                dm_unweighted_file,
                dm_weighted_file,
                dm_euclidean_file
            )
            job.samples = self.samples

            jobs.append(concat_jobs([
                # Create an output directory
                Job(command="mkdir -p " + dm_directory),
                job
            ], name="beta_diversity." + re.sub("_otus", ".", otu_directory) + method))

        return jobs

    def pcoa(self):
        """
        2nd step (/3) for 2D PCoA plot.
        Compute coordinates pour PCoA

        This step takes as input file:

        1. Matrix produced in the previous step.

        """

        jobs = []

        beta_directories = ['open_ref_beta_diversity', 'denovo_beta_diversity', 'closed_ref_beta_diversity']

        for method in ['css', 'rarefaction']:

            if config.param('pcoa', 'beta_diversity_metric') == 'unifrac':
                dm_unweighted_file = self.select_input_files([os.path.join(beta_directory, method, "dissimilarity_matrix", "unweighted_unifrac_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                dm_weighted_file = self.select_input_files([os.path.join(beta_directory, method, "dissimilarity_matrix", "weighted_unifrac_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                dm_euclidean_file = ''
                dm_directory = os.path.dirname(dm_unweighted_file)
            else:
                dm_unweighted_file = ''
                dm_weighted_file = ''
                dm_euclidean_file = self.select_input_files([os.path.join(beta_directory, method, "dissimilarity_matrix", "euclidean_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                dm_directory = os.path.dirname(dm_euclidean_file)

            beta_directory = os.path.dirname(os.path.dirname(dm_directory))

            pcoa_directory = os.path.join(beta_directory, method, "principal_coordinates")
            pcoa_unweighted_file = os.path.join(pcoa_directory, "pcoa_unweighted_unifrac_otu_normalized_table.txt")
            pcoa_weighted_file = os.path.join(pcoa_directory, "pcoa_weighted_unifrac_otu_normalized_table.txt")
            pcoa_euclidean_file = os.path.join(pcoa_directory, "pcoa_euclidean_otu_normalized_table.txt")

            job = qiime.pcoa(
                config.param('pcoa', 'beta_diversity_metric'),
                dm_unweighted_file,
                dm_weighted_file,
                dm_euclidean_file,
                dm_directory,
                pcoa_directory,
                pcoa_unweighted_file,
                pcoa_weighted_file,
                pcoa_euclidean_file
            )
            job.samples = self.samples

            jobs.append(concat_jobs([
                # Create an output directory
                Job(command="mkdir -p " + pcoa_directory),
                job
            ], name="pcoa." + re.sub("_beta_diversity", ".", beta_directory) + method))

        return jobs

    def pcoa_plot(self):
        """
        Last step for 2D PCoA plot.

        This step takes as input file:

        1. PCoA from the previous step.

        """

        jobs = []

        beta_directories = ['open_ref_beta_diversity', 'denovo_beta_diversity', 'closed_ref_beta_diversity']

        for method in ['css', 'rarefaction']:

            if config.param('pcoa_plot', 'beta_diversity_metric') == 'unifrac':
                pcoa_unweighted_file = self.select_input_files([os.path.join(beta_directory, method, "principal_coordinates", "pcoa_unweighted_unifrac_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                pcoa_weighted_file = self.select_input_files([os.path.join(beta_directory, method, "principal_coordinates", "pcoa_weighted_unifrac_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                pcoa_euclidean_file = ''
                pcoa_directory = os.path.dirname(pcoa_unweighted_file)
            else:
                pcoa_unweighted_file = ''
                pcoa_weighted_file = ''
                pcoa_euclidean_file = self.select_input_files([os.path.join(beta_directory, method, "principal_coordinates", "pcoa_euclidean_otu_normalized_table.txt")] for beta_directory in beta_directories)[0]
                pcoa_directory = os.path.dirname(pcoa_euclidean_file)

            beta_directory = os.path.dirname(os.path.dirname(pcoa_directory))

            pcoa_plot_directory = os.path.join(beta_directory, method, "2d_plots")
            beta_diversity_pcoa_unweighted = os.path.join(pcoa_plot_directory, "pcoa_unweighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
            beta_diversity_pcoa_weighted = os.path.join(pcoa_plot_directory, "pcoa_weighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
            beta_diversity_pcoa_euclidean = os.path.join(pcoa_plot_directory, "pcoa_euclidean_otu_normalized_table_2D_PCoA_plots.html")

            if config.param('pcoa_plot', 'map_file'):
                map_file = config.param('pcoa_plot', 'map_file')
            else:
                map_file = "map.txt"

            if config.param('pcoa_plot', 'beta_diversity_metric') == 'unifrac':

                job1 = qiime.pcoa_plot(
                    pcoa_unweighted_file,
                    map_file,
                    beta_diversity_pcoa_unweighted,
                    pcoa_plot_directory
                )
                job1.samples = self.samples

                job2 = qiime.pcoa_plot(
                    pcoa_weighted_file,
                    map_file,
                    beta_diversity_pcoa_weighted,
                    pcoa_plot_directory
                )

                jobs.append(concat_jobs([
                    # Create an output directory
                    Job(command="mkdir -p " + pcoa_plot_directory),
                    job1,
                    job2
                ], name="pcoa_plot." + re.sub("_beta_diversity", ".", beta_directory) + method))

            else:

                job = qiime.pcoa_plot(
                    pcoa_euclidean_file,
                    map_file,
                    beta_diversity_pcoa_euclidean,
                    pcoa_plot_directory
                )
                job.samples = self.samples

                jobs.append(concat_jobs([
                    # Create an output directory
                    Job(command="mkdir -p " + pcoa_plot_directory),
                    job
                ], name="pcoa_plot." + re.sub("_beta_diversity", ".", beta_directory) + method))

        return jobs

    def plot_to_beta(self):
        """
        Final report's 2nd part for the Amplicon-Seq pipeline. Display results (beta diversity PCoA plots).
        """

        jobs = []

        beta_directories = ['open_ref_beta_diversity', 'denovo_beta_diversity', 'closed_ref_beta_diversity']

        for method in ['css', 'rarefaction']:

            if config.param('plot_to_beta', 'beta_diversity_metric') == 'unifrac':
                beta_diversity_pcoa_unweighted = self.select_input_files([os.path.join(beta_directory, method, "2d_plots", "pcoa_unweighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")] for beta_directory in beta_directories)[0]
                beta_diversity_pcoa_weighted = self.select_input_files([os.path.join(beta_directory, method, "2d_plots", "pcoa_weighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")] for beta_directory in beta_directories)[0]

                inputs = [beta_diversity_pcoa_unweighted, beta_diversity_pcoa_weighted]
                beta_directory = os.path.dirname(os.path.dirname(os.path.dirname(inputs[0])))

                description_metric_t = "Beta diversity is a measure of diversity between samples. Using phylogenetic information, the Amplicon-Seq pipeline provides unweighted and weighted UNIFRAC distance metrics PCoA plots."
                link_metric1_t = "Unweighted UNIFRAC distance ([Interactive html plots available here](fig/" + beta_directory + "/" + method + "/2d_plots/pcoa_unweighted_unifrac_otu_normalized_table_2D_PCoA_plots.html))"
                link_metric2_t = "Weighted UNIFRAC distance ([Interactive html plots available here](fig/" + beta_directory + "/" + method + "/2d_plots/pcoa_weighted_unifrac_otu_normalized_table_2D_PCoA_plots.html))"

            else:
                beta_diversity_pcoa_euclidean = self.select_input_files([os.path.join(beta_directory, method, "2d_plots", "pcoa_euclidean_otu_normalized_table_2D_PCoA_plots.html")] for beta_directory in beta_directories)[0]

                inputs = [beta_diversity_pcoa_euclidean]
                beta_directory = os.path.dirname(os.path.dirname(os.path.dirname(inputs[0])))

                description_metric_t = "Beta diversity is a measure of diversity between samples. The Amplicon-Seq pipeline provides euclidean distance metrics PCoA plots."
                link_metric1_t = "Euclidean distance ([Interactive html plots available here](fig/" + beta_directory + "/" + method + "/2d_plots/pcoa_euclidean_otu_normalized_table_2D_PCoA_plots.html))"
                link_metric2_t = " "

            report_file_beta = os.path.join("report", "AmpliconSeq.plot_to_beta_" + method[:3] + ".md")
            report_file_alpha = re.sub("beta", "alpha", report_file_beta)
            report_file = re.sub('plot_to_beta', 'taxonomic_affiliation', report_file_beta)

            inputs.append(report_file_alpha)

            jobs.append(Job(
                inputs,
                [report_file],
                [['plot_to_beta', 'module_pandoc']],
                command="""\
mkdir -p report/fig/{beta_directory}/{method} && \\
cp -r {beta_directory}/{method}/2d_plots/ report/fig/{beta_directory}/{method}/2d_plots/ && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file_beta} \\
  --variable description_metric="{description_metric}" \\
  --variable link_metric1="{link_metric1}" \\
  --variable link_metric2="{link_metric2}" \\
  {report_template_dir}/{basename_report_file_beta} \\
  > {report_file_beta} && \\
cat {report_file_alpha} {report_file_beta} > {report_file}""".format(
                    beta_directory=beta_directory,
                    method=method,
                    description_metric=description_metric_t,
                    link_metric1=link_metric1_t,
                    link_metric2=link_metric2_t,
                    report_template_dir=self.report_template_dir,
                    basename_report_file_beta=re.sub('_'+method[:3], '', os.path.basename(report_file_beta)),
                    report_file_alpha=report_file_alpha,
                    report_file_beta=report_file_beta,
                    report_file=report_file
                ),
                report_files=[report_file],
                samples=self.samples,
                name="plot_to_beta." + re.sub("_beta_diversity", ".", beta_directory) + method
            ))

        return jobs

    @property
    def steps(self):
        return [
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.flash,
            self.merge_flash_stats,
            self.catenate,                  #5
            self.uchime,
            self.merge_uchime_stats,
            self.otu_picking,
            self.otu_rep_picking,
            self.otu_assigning,             #10
            self.otu_table,
            self.otu_alignment,
            self.filter_alignment,
            self.phylogeny,
            self.qiime_report,              #15
            self.multiple_rarefaction,
            self.alpha_diversity,
            self.collate_alpha,
            self.sample_rarefaction_plot,
            self.qiime_report2,             #20
            self.single_rarefaction,
            self.css_normalization,
            self.rarefaction_plot,
            self.summarize_taxa,
            self.plot_taxa,                 #25
            self.plot_heatmap,
            self.krona,
            self.plot_to_alpha,
            self.beta_diversity,
            self.pcoa,                      #30
            self.pcoa_plot,
            self.plot_to_beta
        ]

if __name__ == '__main__':
    AmpliconSeq()
