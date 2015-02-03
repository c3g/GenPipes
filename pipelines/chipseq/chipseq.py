#!/usr/bin/env python

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
        return jobs

    # Generic function called in metrics step
    def merge_metrics(self, name, alignment_file,flagstat_file, metrics_file):
        return concat_jobs([
                samtools.flagstat(alignment_file, flagstat_file),
                Job(
                    [flagstat_file],
                    [metrics_file],
                    [['metrics', 'module_perl'], ['metrics', 'module_mugqic_tools']],
                    command="""\
perl -MReadMetrics -e 'ReadMetrics::mergeStats(
  "{name}",
  "{metrics_file}",
  ReadMetrics::parseFlagstats(
    "{name}",
    "{flagstat_file}"
  )
)'""".format(
                        name=name,
                        metrics_file=metrics_file,
                        flagstat_file=flagstat_file
                    ),
                    removable_files=[metrics_file]
                )
            ])

    def metrics(self):
        """
        The number of raw/filtered and aligned reads per sample are computed at this stage.
        """

        jobs = []
        trimming_stats = os.path.join("metrics", "trimming.stats")
        for readset in self.readsets:
            alignment_file = os.path.join("alignment", readset.sample.name, readset.name, readset.name + ".sorted.bam")
            flagstat_file = alignment_file + ".flagstat"
            metrics_file = os.path.join("metrics", readset.name + ".readstats.csv")

            job = concat_jobs([
                self.merge_metrics(readset.sample.name + " " + readset.name, alignment_file, flagstat_file, metrics_file + ".tmp"),
                Job(
                    [trimming_stats, metrics_file + ".tmp"],
                    [metrics_file],
                    command="""\
grep "{sample_name}\t{readset_name}" {trimming_stats} | cut -f3- | sed 's/\t/,/g' | sed '1iRaw Fragments,Fragment Surviving,Single Surviving' | paste -d, <(cut -f1 -d, {metrics_file}.tmp) - <(cut -f2- -d, {metrics_file}.tmp) \\
  > {metrics_file} && \\
rm {metrics_file}.tmp
""".format(
                        sample_name=readset.sample.name,
                        readset_name=readset.name,
                        trimming_stats=trimming_stats,
                        metrics_file=metrics_file
                    ),
                    removable_files=[metrics_file]
                )
            ], name = "metrics." + readset.name)
            # Remove uncessary temporary file otherwise job is never up to date
            job.output_files.remove(metrics_file + ".tmp")
            jobs.append(job)

        for sample in self.samples:
            alignment_file = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam")
            flagstat_file = alignment_file + ".flagstat"
            metrics_file = os.path.join("metrics", sample.name + ".memstats.csv")

            job = self.merge_metrics(sample.name, alignment_file, flagstat_file, metrics_file)
            job.name = "metrics." + sample.name
            jobs.append(job)
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

        return jobs

    def qc_metrics(self):
        """
        Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated.
        """

         # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)

        output_files = [os.path.join("graphs", sample.name + "_QC_Metrics.ps") for sample in self.samples]

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
  {output_dir}""".format(
                design_file=design_file,
                output_dir=self.output_dir
            ),
            name="qc_plots_R",
            removable_files=output_files
        )]

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

        output_files.append(os.path.join("annotation", "peak_stats.csv"))

        return [Job(
            input_files,
            output_files,
            [
                ['annotation_graphs', 'module_mugqic_tools'],
                ['annotation_graphs', 'module_R']
            ],
            command="""\
mkdir -p graphs && \\
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \\
  {design_file} \\
  {output_dir}""".format(
                design_file=design_file,
                output_dir=self.output_dir
            ),
            name="annotation_graphs",
            removable_files=output_files
        )]

    def gq_seq_utils_report(self):
        """
        Generate the standard report. A summary html report is automatically generated by the pipeline.
        This report contains description of the sequencing experiment as well as a detailed presentation
        of the pipeline steps and results. Various Quality Control (QC) summary statistics are included
        in the report and additional QC analysis is accessible for download directly through the report.
        The report includes also the main references of the software and methods used during the analysis,
        together with the full list of parameters passed to the pipeline main script.
        """

        job = gq_seq_utils.report(
            [config_file.name for config_file in self.args.config],
            self.output_dir,
            "CHIPseq",
            self.output_dir
        )

        # Mammouth compute nodes don't have 'convert' command by default. ImageMagick module must be loaded explicitely
        if config.param('gq_seq_utils_report', 'module_imagemagick', required=False):
            job.modules.append(config.param('gq_seq_utils_report', 'module_imagemagick', required=False))

        job.input_files = \
            [os.path.join("metrics", readset.name + ".readstats.csv") for readset in self.readsets] + \
            [os.path.join("metrics", sample.name + ".memstats.csv") for sample in self.samples] + \
            [os.path.join("graphs", sample.name + "_QC_Metrics.ps") for sample in self.samples] + \
            [os.path.join("annotation", contrast.real_name, contrast.real_name + ".annotated.csv") for contrast in self.contrasts if contrast.type == 'narrow'] + \
            [os.path.join("annotation", contrast.real_name, contrast.real_name, "GenomeOntology.html") for contrast in self.contrasts if contrast.type == 'narrow'] + \
            [os.path.join("annotation", contrast.real_name, contrast.real_name, "geneOntology.html") for contrast in self.contrasts if contrast.type == 'narrow'] + \
            [os.path.join("annotation", contrast.real_name, contrast.real_name, "homerResults.html") for contrast in self.contrasts if contrast.type == 'narrow'] + \
            [os.path.join("annotation", contrast.real_name, contrast.real_name, "knownResults.html") for contrast in self.contrasts if contrast.type == 'narrow'] + \
            [os.path.join("graphs", contrast.real_name + "_Misc_Graphs.ps") for contrast in self.contrasts] + \
            [os.path.join("annotation", "peak_stats.csv")] if [contrast for contrast in self.contrasts if contrast.type == 'narrow'] else []
        job.name = "gq_seq_utils_report"
        job.removable_files = [os.path.join("graphs", sample.name + "_QC_Metrics.png") for sample in self.samples] + \
            [os.path.join("graphs", contrast.real_name + "_Misc_Graphs.png") for contrast in self.contrasts] + \
            [os.path.join("graphs", contrast.real_name + "_Misc_Graphs.pdf") for contrast in self.contrasts]
        return [job]

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
            self.homer_make_ucsc_file,
            self.qc_metrics,
            self.macs2_callpeak,
            self.homer_annotate_peaks,
            self.homer_find_motifs_genome,
            self.annotation_graphs,
            self.gq_seq_utils_report
        ]

if __name__ == '__main__': 
    ChipSeq()
