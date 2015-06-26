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
from bfx.design import *

from bfx import gq_seq_utils
from bfx import picard
from bfx import samtools
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class ChipSeq(dnaseq.DnaSeq):
    """
    ChIP-Seq Pipeline
    =================

    ChIP-Seq experiments allows the Isolation and sequencing of genomic DNA bound by a specific transcription factor,
    covalently modified histone, or other nuclear protein. The pipeline starts by trimming adaptors and
    low quality bases and mapping the reads (single end or paired end ) to a reference genome using bwa.
    Reads are filtered by mapping quality and duplicate reads are marked. Then, Homer quality control routines
    are used to provide information and feedback about the quality of the experiment. Peak calls is executed by MACS
    and annotation and motif discovery for narrow peaks are executed using Homer. Statistics of annotated peaks
    are produced for narrow peaks and a standard report is generated.

    An example of the ChIP-Seq report for an analysis on public ENCODE data is available for illustration purpose only:
    [ChIP-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/chipReport/index.html).

    [Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_ChIP-Seq.pptx)
    is more information about ChIP-Seq pipeline that you may find interesting.
    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(ChipSeq, self).__init__()

    @property
    def contrasts(self):
        contrasts = super(ChipSeq, self).contrasts

        # Parse contrasts to retrieve name and type
        for contrast in contrasts:
            if re.search("^\w[\w.-]*,[BN]$", contrast.name):
                contrast.real_name = contrast.name.split(",")[0]
                if contrast.name.split(",")[1] == 'B':
                    contrast.type = 'broad'
                elif contrast.name.split(",")[1] == 'N':
                    contrast.type = 'narrow'
            else:
                raise Exception("Error: contrast name \"" + contrast.name + "\" is invalid (should be <contrast>,B for broad or <contrast>,N for narrow)!")

        return contrasts

    def mappable_genome_size(self):
        genome_index = csv.reader(open(config.param('DEFAULT', 'genome_fasta', type='filepath') + ".fai", 'rb'), delimiter='\t')
        # 2nd column of genome index contains chromosome length
        # HOMER and MACS2 mappable genome size (without repetitive features) is about 80 % of total size
        return sum([int(chromosome[1]) for chromosome in genome_index]) * 0.8

    def samtools_view_filter(self):
        """
        Filter unique reads by mapping quality using [Samtools](http://www.htslib.org/).
        """

        jobs = []
        for readset in self.readsets:
            readset_bam_prefix = os.path.join("alignment", readset.sample.name, readset.name, readset.name + ".sorted.")
            readset_bam = readset_bam_prefix + "bam"
            filtered_readset_bam = readset_bam_prefix + "filtered.bam"

            job = samtools.view(readset_bam, filtered_readset_bam, "-b -F4 -q " + str(config.param('samtools_view_filter', 'min_mapq', type='int')))
            job.name = "samtools_view_filter." + readset.name
            jobs.append(job)

        report_file = os.path.join("report", "ChipSeq.samtools_view_filter.md")
        jobs.append(
            Job(
                [os.path.join("alignment", readset.sample.name, readset.name, readset.name + ".sorted.filtered.bam") for readset in self.readsets],
                [report_file],
                [['samtools_view_filter', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable min_mapq="{min_mapq}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    min_mapq=config.param('samtools_view_filter', 'min_mapq', type='int'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="samtools_view_filter_report")
        )

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
            readset_bams = [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.filtered.bam") for readset in sample.readsets]
            sample_bam = os.path.join(alignment_directory, sample.name + ".merged.bam")

            mkdir_job = Job(command="mkdir -p " + os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, alignment_directory)

                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + target_readset_bam + " " + sample_bam, removable_files=[sample_bam]),
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:
                job = concat_jobs([
                    mkdir_job,
                    picard.merge_sam_files(readset_bams, sample_bam)
                ])
                job.name = "picard_merge_sam_files." + sample.name

            jobs.append(job)

        return jobs

    def picard_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "merged.bam"
            output = alignment_file_prefix + "sorted.dup.bam"
            metrics_file = alignment_file_prefix + "sorted.dup.metrics"

            job = picard.mark_duplicates([input], output, metrics_file)
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)

        report_file = os.path.join("report", "ChipSeq.picard_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in self.samples],
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
                name="picard_mark_duplicates_report")
        )

        return jobs

    def metrics(self):
        """
        The number of raw/filtered and aligned reads per sample are computed at this stage.
        """

        jobs = []
        jobs.append(concat_jobs([samtools.flagstat(os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam"), os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam.flagstat")) for sample in self.samples], name="metrics.flagstat"))

        trim_metrics_file = os.path.join("metrics", "trimSampleTable.tsv")
        metrics_file = os.path.join("metrics", "SampleMetrics.stats")
        report_metrics_file = os.path.join("report", "trimMemSampleTable.tsv")
        report_file = os.path.join("report", "ChipSeq.metrics.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam.flagstat") for sample in self.samples],
                [report_metrics_file],
                [['metrics', 'module_pandoc']],
                # Retrieve number of aligned and duplicate reads from sample flagstat files
                # Merge trimming stats per sample with aligned and duplicate stats using ugly awk
                # Format merge stats into markdown table using ugly awk (knitr may do this better)
                command="""\
for sample in {samples}
do
  flagstat_file=alignment/$sample/$sample.sorted.dup.bam.flagstat
  echo -e "$sample\t`grep -P '^\d+ \+ \d+ mapped' $flagstat_file | grep -Po '^\d+'`\t`grep -P '^\d+ \+ \d+ duplicate' $flagstat_file | grep -Po '^\d+'`"
done | \\
awk -F"\t" '{{OFS="\t"; print $0, $3 / $2 * 100}}' | sed '1iSample\tAligned Filtered Reads\tDuplicate Reads\tDuplicate %' \\
  > {metrics_file} && \\
mkdir -p report && \\
if [[ -f {trim_metrics_file} ]]
then
  awk -F "\t" 'FNR==NR{{trim_line[$1]=$0; surviving[$1]=$3; next}}{{OFS="\t"; if ($1=="Sample") {{print trim_line[$1], $2, "Aligned Filtered %", $3, $4}} else {{print trim_line[$1], $2, $2 / surviving[$1] * 100, $3, $4}}}}' {trim_metrics_file} {metrics_file} \\
  > {report_metrics_file}
else
  cp {metrics_file} {report_metrics_file}
fi && \\
trim_mem_sample_table=`if [[ -f {trim_metrics_file} ]] ; then LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4), sprintf("%\\47d", $5), sprintf("%.1f", $6), sprintf("%\\47d", $7), sprintf("%.1f", $8)}}}}' {report_metrics_file} ; else LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4)}}}}' {report_metrics_file} ; fi` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trim_mem_sample_table="$trim_mem_sample_table" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}
""".format(
                    samples=" ".join([sample.name for sample in self.samples]),
                    trim_metrics_file=trim_metrics_file,
                    metrics_file=metrics_file,
                    report_metrics_file=report_metrics_file,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                name="metrics_report",
                removable_files=[report_metrics_file],
                report_files=[report_file]
            )
        )
        return jobs

    def homer_make_tag_directory(self):
        """
        The Homer Tag directories, used to check for quality metrics, are computed at this step. 
        """

        jobs = []
        for sample in self.samples:
            alignment_file = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam")
            output_dir = os.path.join("tags", sample.name)

            jobs.append(Job(
                [alignment_file],
                [os.path.join(output_dir, "tagInfo.txt")],
                [['homer_make_tag_directory', 'module_samtools'], ['homer_make_tag_directory', 'module_homer']],
                command="""\
makeTagDirectory \\
  {output_dir} \\
  {alignment_file} \\
  -checkGC -genome {genome}""".format(
                    output_dir=output_dir,
                    alignment_file=alignment_file,
                    genome=config.param('homer_make_tag_directory', 'assembly')
                ),
                name="homer_make_tag_directory." + sample.name,
                removable_files=[output_dir]
            ))

        return jobs

    def qc_metrics(self):
        """
        Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated.
        """

         # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)

        report_file = os.path.join("report", "ChipSeq.qc_metrics.md")
        output_files = [os.path.join("graphs", sample.name + "_QC_Metrics.ps") for sample in self.samples] + [report_file]

        return [Job(
            [os.path.join("tags", sample.name, "tagInfo.txt") for sample in self.samples],
            output_files,
            [
                ['qc_plots_R', 'module_mugqic_tools'],
                ['qc_plots_R', 'module_R']
            ],
            command="""\
mkdir -p graphs && \\
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \\
  {design_file} \\
  {output_dir} && \\
cp {report_template_dir}/{basename_report_file} {report_file} && \\
for sample in {samples}
do
  cp --parents graphs/${{sample}}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${{sample}}_QC_Metrics.ps report/graphs/${{sample}}_QC_Metrics.png
  echo -e "----\n\n![QC Metrics for Sample $sample ([download high-res image](graphs/${{sample}}_QC_Metrics.ps))](graphs/${{sample}}_QC_Metrics.png)\n" \\
  >> {report_file}
done""".format(
                samples=" ".join([sample.name for sample in self.samples]),
                design_file=design_file,
                output_dir=self.output_dir,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            ),
            name="qc_plots_R",
            removable_files=output_files,
            report_files=[report_file]
        )]

    def homer_make_ucsc_file(self):
        """
        Wiggle Track Format files are generated from the aligned reads using Homer.
        The resulting files can be loaded in browsers like IGV or UCSC.
        """

        jobs = []
        for sample in self.samples:
            tag_dir = os.path.join("tags", sample.name)
            bedgraph_dir = os.path.join("tracks", sample.name)
            bedgraph_file = os.path.join(bedgraph_dir, sample.name + ".ucsc.bedGraph.gz")

            jobs.append(Job(
                [os.path.join(tag_dir, "tagInfo.txt")],
                [bedgraph_file],
                [['homer_make_ucsc_files', 'module_homer']],
                command="""\
mkdir -p {bedgraph_dir} && \\
makeUCSCfile \\
  {tag_dir} | \\
gzip -1 -c > {bedgraph_file}""".format(
                    tag_dir=tag_dir,
                    bedgraph_dir=bedgraph_dir,
                    bedgraph_file=bedgraph_file
                ),
                name="homer_make_ucsc_file." + sample.name,
                removable_files=[bedgraph_dir]
            ))

        report_file = os.path.join("report", "ChipSeq.homer_make_ucsc_file.md")
        jobs.append(
            Job(
                [os.path.join("tracks", sample.name, sample.name + ".ucsc.bedGraph.gz")],
                [report_file],
                command="""\
mkdir -p report && \\
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \\
cp {report_template_dir}/{basename_report_file} report/""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="homer_make_ucsc_file_report")
        )

        return jobs

    def macs2_callpeak(self):
        """
        Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
        The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
        The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
        The default mfold parameter of MACS2 is [10,30].
        """

        jobs = []

        for contrast in self.contrasts:
            if contrast.treatments:
                treatment_files = [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in contrast.treatments]
                control_files = [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in contrast.controls]
                output_dir = os.path.join("peak_call", contrast.real_name)

                if contrast.type == 'broad':  # Broad region
                  other_options = " --broad --nomodel"
                else:  # Narrow region
                    if control_files:
                        other_options = " --nomodel"
                    else:
                        other_options = " --fix-bimodal"

                jobs.append(Job(
                    treatment_files + control_files,
                    [os.path.join(output_dir, contrast.real_name + "_peaks." + contrast.type + "Peak")],
                    [['macs2_callpeak', 'module_python'], ['macs2_callpeak', 'module_macs2']],
                    command="""\
mkdir -p {output_dir} && \\
macs2 callpeak {format}{other_options} \\
  --gsize {genome_size} \\
  --treatment \\
  {treatment_files}{control_files} \\
  --name {output_prefix_name} \\
  >& {output_prefix_name}.diag.macs.out""".format(
                        output_dir=output_dir,
                        format="--format " + ("BAMPE" if self.run_type == "PAIRED_END" else "BAM"),
                        other_options=other_options,
                        genome_size=self.mappable_genome_size(),
                        treatment_files=" \\\n  ".join(treatment_files),
                        control_files=" \\\n  --control \\\n  " + " \\\n  ".join(control_files) if control_files else " \\\n  --nolambda",
                        output_prefix_name=os.path.join(output_dir, contrast.real_name)
                    ),
                    name="macs2_callpeak." + contrast.real_name,
                    removable_files=[output_dir]
                ))
            else:
                log.warning("No treatment found for contrast " + contrast.name + "... skipping")

        report_file = os.path.join("report", "ChipSeq.macs2_callpeak.md")
        jobs.append(
            Job(
                [os.path.join("peak_call", contrast.real_name, contrast.real_name + "_peaks." + contrast.type + "Peak") for contrast in self.contrasts],
                [report_file],
                command="""\
mkdir -p report && \\
cp {report_template_dir}/{basename_report_file} report/ && \\
for contrast in {contrasts}
do
  cp -a --parents peak_call/$contrast/ report/ && \\
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${{contrast}}_peaks.xls)" \\
  >> {report_file}
done""".format(
                    contrasts=" ".join([contrast.real_name for contrast in self.contrasts]),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="macs2_callpeak_report")
        )

        return jobs

    def homer_annotate_peaks(self):
        """
        The peaks called previously are annotated with HOMER using RefSeq annotations for the reference genome.
        Gene ontology and genome ontology analysis are also performed at this stage.
        """

        jobs = []

        for contrast in self.contrasts:
            if contrast.treatments:
                peak_file = os.path.join("peak_call", contrast.real_name, contrast.real_name + "_peaks." + contrast.type + "Peak")
                output_prefix = os.path.join("annotation", contrast.real_name, contrast.real_name)
                annotation_file = output_prefix + ".annotated.csv"

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + output_prefix),
                    Job(
                        [peak_file],
                        [
                            annotation_file,
                            os.path.join(output_prefix, "geneOntology.html"),
                            os.path.join(output_prefix, "GenomeOntology.html")
                        ],
                        [['homer_annotate_peaks', 'module_perl'], ['homer_annotate_peaks', 'module_homer']],
                        command="""\
annotatePeaks.pl \\
  {peak_file} \\
  {genome} \\
  -gsize {genome} \\
  -cons -CpG \\
  -go {output_prefix} \\
  -genomeOntology {output_prefix} \\
  > {annotation_file}""".format(
                            peak_file=peak_file,
                            genome=config.param('homer_annotate_peaks', 'assembly'),
                            genome_size=config.param('homer_annotate_peaks', 'assembly'),
                            output_prefix=output_prefix,
                            annotation_file=annotation_file
                        )
                    ),
                    Job(
                        [annotation_file],
                        [
                            output_prefix + ".tss.stats.csv",
                            output_prefix + ".exon.stats.csv",
                            output_prefix + ".intron.stats.csv",
                            output_prefix + ".tss.distance.csv"
                        ],
                        [['homer_annotate_peaks', 'module_perl'], ['homer_annotate_peaks', 'module_mugqic_tools']],
                        command="""\
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "{annotation_file}",
  "{output_prefix}",
  {proximal_distance},
  {distal_distance},
  {distance5d_lower},
  {distance5d_upper},
  {gene_desert_size}
)'""".format(
                            annotation_file=annotation_file,
                            output_prefix=output_prefix,
                            proximal_distance=config.param('homer_annotate_peaks', 'proximal_distance', type='int'),
                            distal_distance=config.param('homer_annotate_peaks', 'distal_distance', type='int'),
                            distance5d_lower=config.param('homer_annotate_peaks', 'distance5d_lower', type='int'),
                            distance5d_upper=config.param('homer_annotate_peaks', 'distance5d_upper', type='int'),
                            gene_desert_size=config.param('homer_annotate_peaks', 'gene_desert_size', type='int')
                        ),
                        removable_files=[os.path.join("annotation", contrast.real_name)]
                    )
                ], name="homer_annotate_peaks." + contrast.real_name))
            else:
                log.warning("No treatment found for contrast " + contrast.name + "... skipping")

        report_file = os.path.join("report", "ChipSeq.homer_annotate_peaks.md")
        jobs.append(
            Job(
                [os.path.join("annotation", contrast.real_name, contrast.real_name + ".annotated.csv") for contrast in self.contrasts],
                [report_file],
                command="""\
mkdir -p report/annotation/ && \\
cp {report_template_dir}/{basename_report_file} report/ && \\
for contrast in {contrasts}
do
  rsync -avP annotation/$contrast report/annotation/ && \\
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${{contrast}}.annotated.csv)\n* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)\n* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \\
  >> {report_file}
done""".format(
                    contrasts=" ".join([contrast.real_name for contrast in self.contrasts]),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="homer_annotate_peaks_report")
        )

        return jobs

    def homer_find_motifs_genome(self):
        """
        De novo and known motif analysis per design are performed using HOMER.
        """

        jobs = []

        for contrast in self.contrasts:
            # Don't find motifs for broad peaks
            if contrast.type == 'narrow' and contrast.treatments:
                peak_file = os.path.join("peak_call", contrast.real_name, contrast.real_name + "_peaks." + contrast.type + "Peak")
                output_dir = os.path.join("annotation", contrast.real_name, contrast.real_name)

                jobs.append(Job(
                    [peak_file],
                    [
                        os.path.join(output_dir, "homerResults.html"),
                        os.path.join(output_dir, "knownResults.html")
                    ],
                    [
                        ['homer_find_motifs_genome', 'module_perl'],
                        ['homer_find_motifs_genome', 'module_weblogo'],
                        ['homer_find_motifs_genome', 'module_homer']
                    ],
                    command="""\
mkdir -p {output_dir} && \\
findMotifsGenome.pl \\
  {peak_file} \\
  {genome} \\
  {output_dir} \\
  -preparsedDir {output_dir}/preparsed \\
  -p {threads}""".format(
                        peak_file=peak_file,
                        genome=config.param('homer_find_motifs_genome', 'assembly'),
                        output_dir=output_dir,
                        threads=config.param('homer_find_motifs_genome', 'threads', type='posint')
                    ),
                    name="homer_find_motifs_genome." + contrast.real_name,
                    removable_files=[os.path.join("annotation", contrast.real_name)]
                ))
            else:
                log.warning("No treatment found for contrast " + contrast.name + "... skipping")

        report_file = os.path.join("report", "ChipSeq.homer_find_motifs_genome.md")
        jobs.append(
            Job(
                [os.path.join("annotation", contrast.real_name, contrast.real_name, "homerResults.html") for contrast in self.contrasts if contrast.type == 'narrow' and contrast.treatments] +
                [os.path.join("annotation", contrast.real_name, contrast.real_name, "knownResults.html") for contrast in self.contrasts if contrast.type == 'narrow' and contrast.treatments],
                [report_file],
                command="""\
mkdir -p report/annotation/ && \\
cp {report_template_dir}/{basename_report_file} report/ && \\
for contrast in {contrasts}
do
  rsync -avP annotation/$contrast report/annotation/ && \\
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)\n* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \\
  >> {report_file}
done""".format(
                    contrasts=" ".join([contrast.real_name for contrast in self.contrasts if contrast.type == 'narrow' and contrast.treatments]),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="homer_find_motifs_genome_report")
        )

        return jobs

    def annotation_graphs(self):
        """
        The peak location statistics. The following peak location statistics are generated per design:
        proportions of the genomic locations of the peaks. The locations are: Gene (exon or intron),
        Proximal ([0;2] kb upstream of a transcription start site), Distal ([2;10] kb upstream
        of a transcription start site), 5d ([10;100] kb upstream of a transcription start site),
        Gene desert (>= 100 kb upstream or downstream of a transcription start site), Other (anything
        not included in the above categories); The distribution of peaks found within exons and introns;
        The distribution of peak distance relative to the transcription start sites (TSS);
        the Location of peaks per design.
        """

         # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)

        input_files = []
        output_files = []
        for contrast in self.contrasts:
            annotation_prefix = os.path.join("annotation", contrast.real_name, contrast.real_name)
            input_files.append(annotation_prefix + ".tss.stats.csv")
            input_files.append(annotation_prefix + ".exon.stats.csv")
            input_files.append(annotation_prefix + ".intron.stats.csv")
            input_files.append(annotation_prefix + ".tss.distance.csv")

            output_files.append(os.path.join("graphs", contrast.real_name + "_Misc_Graphs.ps"))

        peak_stats_file = os.path.join("annotation", "peak_stats.csv")
        output_files.append(peak_stats_file)
        report_file = os.path.join("report", "ChipSeq.annotation_graphs.md")
        output_files.append(report_file)

        return [Job(
            input_files,
            output_files,
            [
                ['annotation_graphs', 'module_mugqic_tools'],
                ['annotation_graphs', 'module_R'],
                ['annotation_graphs', 'module_pandoc']
            ],
            command="""\
mkdir -p graphs && \\
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \\
  {design_file} \\
  {output_dir} && \\
mkdir -p report/annotation/ && \\
if [[ -f {peak_stats_file} ]]
then
  cp {peak_stats_file} report/annotation/
peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, $2,  sprintf("%\\47d", $3), $4, sprintf("%\\47.1f", $5), sprintf("%\\47.1f", $6), sprintf("%\\47.1f", $7), sprintf("%\\47.1f", $8)}}}}' {peak_stats_file}`
else
  peak_stats_table=""
fi
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable peak_stats_table="$peak_stats_table" \\
  --variable proximal_distance="{proximal_distance}" \\
  --variable distal_distance="{distal_distance}" \\
  --variable distance5d_lower="{distance5d_lower}" \\
  --variable distance5d_upper="{distance5d_upper}" \\
  --variable gene_desert_size="{gene_desert_size}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file} && \\
for contrast in {contrasts}
do
  cp --parents graphs/${{contrast}}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${{contrast}}_Misc_Graphs.ps report/graphs/${{contrast}}_Misc_Graphs.png
  echo -e "----\n\n![Annotation Statistics for Design $contrast ([download high-res image](graphs/${{contrast}}_Misc_Graphs.ps))](graphs/${{contrast}}_Misc_Graphs.png)\n" \\
  >> {report_file}
done""".format(
                design_file=design_file,
                output_dir=self.output_dir,
                peak_stats_file=peak_stats_file,
                contrasts=" ".join([contrast.real_name for contrast in self.contrasts if contrast.type == 'narrow' and contrast.treatments]),
                proximal_distance=config.param('homer_annotate_peaks', 'proximal_distance', type='int') / -1000,
                distal_distance=config.param('homer_annotate_peaks', 'distal_distance', type='int') / -1000,
                distance5d_lower=config.param('homer_annotate_peaks', 'distance5d_lower', type='int') / -1000,
                distance5d_upper=config.param('homer_annotate_peaks', 'distance5d_upper', type='int') / -1000,
                gene_desert_size=config.param('homer_annotate_peaks', 'gene_desert_size', type='int') / 1000,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            ),
            name="annotation_graphs",
            report_files=[report_file],
            removable_files=output_files
        )]

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.bwa_mem_picard_sort_sam,
            self.samtools_view_filter,
            self.picard_merge_sam_files,
            self.picard_mark_duplicates,
            self.metrics,
            self.homer_make_tag_directory,
            self.qc_metrics,
            self.homer_make_ucsc_file,
            self.macs2_callpeak,
            self.homer_annotate_peaks,
            self.homer_find_motifs_genome,
            self.annotation_graphs
        ]

if __name__ == '__main__': 
    ChipSeq()
