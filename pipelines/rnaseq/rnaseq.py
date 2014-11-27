#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

from bfx import bedtools
from bfx import cufflinks
from bfx import differential_expression
from bfx import htseq
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import star
from bfx import gq_seq_utils
from pipelines import common
import utils

log = logging.getLogger(__name__)

class RnaSeq(common.Illumina):
    
    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(RnaSeq, self).__init__()

    def star(self):
        jobs = []
        project_index_directory = "reference.Merged"
        project_junction_file =  os.path.join("alignment_1stPass", "AllSample.SJ.out.tab")
        individual_junction_list=[]
        ######
        #pass 1 -alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_1stPass_directory = os.path.join("alignment_1stPass", readset.sample.name)
            individual_junction_list.append(os.path.join(alignment_1stPass_directory,"SJ.out.tab"))

            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_1stPass_directory,
                genome_index_folder=None,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else ""
            )
            job.name = "star_align.1." + readset.name
            jobs.append(job)
        
        ######
        #pass 1 - contatenate junction
        job = star.concatenate_junction(
            input_junction_files_list=individual_junction_list,
            output_junction_file=project_junction_file
        )
        job.name = "star_concat.AllSample"
        jobs.append(job)
        
        ######
        #pass 1 - genome indexing
        job = star.index(
            genome_index_folder=project_index_directory,
            junction_file=project_junction_file
        )
        job.name = "star_index.AllSample"
        jobs.append(job)
        
        ######
        #Pass 2 - alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_2ndPass_directory = os.path.join("alignment", readset.sample.name)

            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_2ndPass_directory,
                genome_index_folder=project_index_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else "",
                create_wiggle_track=True,
                search_chimeres=True,
                cuff_follow=True,
                sort_bam=True
            )
            job.input_files.append(os.path.join(project_index_directory, "SAindex"))
 
            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            # remove older symlink before otherwise it raise an error if the link already exist (in case of redo)
            if len(readset.sample.readsets) == 1:
                sample_bam = os.path.join("alignment",readset.sample.name ,readset.sample.name + ".sorted.bam")
                job.command += " && \\\n rm -f " + sample_bam + " && \\\n ln -s " + os.path.join(self.output_dir,"alignment",readset.name, "Aligned.sortedByCoord.out.bam") + " " + sample_bam
                job.output_files.append(sample_bam)

            job.name = "star_align.2." + readset.name
            jobs.append(job)

        return jobs

    def picard_merge_sam_files(self):
        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name, "Aligned.sortedByCoord.out.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def picard_sort_sam(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)

            job = picard.sort_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".QueryNameSorted.bam",
                "queryname"
            )
            job.name = "picard_sort_sam." + sample.name
            jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.")

            job = picard.mark_duplicates(
                [alignment_file_prefix + "bam"],
                alignment_file_prefix + "mdup.bam",
                alignment_file_prefix + "mdup.metrics"
            )
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def rnaseqc(self):
        sample_file = os.path.join("alignment", "rnaseqc.samples.txt")
        sample_rows = [[sample.name, os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam"), "RNAseq"] for sample in self.samples]
        input_bams = [sample_row[1] for sample_row in sample_rows]
        output_directory = os.path.join("metrics", "rnaseqRep")
        # Use GTF with transcript_id only otherwise RNASeQC fails
        gtf_transcript_id = config.param('rnaseqc', 'gtf_transcript_id', type='filepath')

        job = concat_jobs([
            Job(command="mkdir -p " + output_directory),
            Job(input_bams, [sample_file], command="""\
echo "Sample\tBamFile\tNote
{sample_rows}" \\
  > {sample_file}""".format(sample_rows="\n".join(["\t".join(sample_row) for sample_row in sample_rows]), sample_file=sample_file)),
            metrics.rnaseqc(sample_file, output_directory, self.run_type == "SINGLE_END", gtf_file=gtf_transcript_id),
            Job([], [output_directory + ".zip"], command="zip -r {output_directory}.zip {output_directory}".format(output_directory=output_directory))
        ], name="rnaseqc")

        return [job]

    def wiggle(self):
        jobs = []

        for sample in self.samples:
            bam_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.")
            input_bam = bam_file_prefix + "bam"
            bed_graph_prefix = os.path.join("tracks", sample.name, sample.name)
            big_wig_prefix = os.path.join("tracks", "bigWig", sample.name)

            if config.param('DEFAULT', 'strand_info') != 'fr-unstranded':
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
                ], name="wiggle." + sample.name + ".forward_strandspec")
                # Remove temporary-then-deleted files from job output files, otherwise job is never up to date
                bam_f_job.output_files.remove(input_bam_f1)
                bam_f_job.output_files.remove(input_bam_f2)

                bam_r_job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks", "bigWig")),
                    samtools.view(input_bam, input_bam_r1, "-bh -F 256 -f 97"),
                    samtools.view(input_bam, input_bam_r2, "-bh -F 256 -f 145"),
                    picard.merge_sam_files([input_bam_r1, input_bam_r2], output_bam_r),
                    Job(command="rm " + input_bam_r1 + " " + input_bam_r2)
                ], name="wiggle." + sample.name + ".reverse_strandspec")
                # Remove temporary-then-deleted files from job output files, otherwise job is never up to date
                bam_r_job.output_files.remove(input_bam_r1)
                bam_r_job.output_files.remove(input_bam_r2)

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
                ], name="wiggle." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output)))
                jobs.append(job)

        return jobs

    def raw_counts(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)
            input_bam = alignment_file_prefix + ".QueryNameSorted.bam"

            # Count reads
            output_count = os.path.join("raw_counts", sample.name + ".readcounts.csv")
            stranded = "no" if config.param('DEFAULT', 'strand_info') == "fr-unstranded" else "reverse"
            job = concat_jobs([
                Job(command="mkdir -p raw_counts"),
                htseq.htseq_count(
                    input_bam,
                    config.param('htseq', 'gtf', type='filepath'),
                    output_count,
                    config.param('htseq', 'options'),
                    stranded
                )
            ], name="htseq_count." + sample.name)
            jobs.append(job)

        return jobs

    def raw_counts_metrics(self):
        jobs = []

        # Create raw count matrix
        output_directory = "DGE"
        read_count_files = [os.path.join("raw_counts", sample.name + ".readcounts.csv") for sample in self.samples]
        output_matrix = os.path.join(output_directory, "rawCountMatrix.csv")

        job = Job(read_count_files, [output_matrix], [['raw_counts_metrics', 'module_mugqic_tools']], name="metrics.matrix")

        job.command = """\
mkdir -p {output_directory} && \\
gtf2tmpMatrix.awk \\
  {reference_gtf} \\
  {output_directory}/tmpMatrix.txt && \\
HEAD='Gene\tSymbol' && \\
for read_count_file in \\
  {read_count_files}
do
  sort -k1,1 $read_count_file > {output_directory}/tmpSort.txt && \\
  join -1 1 -2 1 {output_directory}/tmpMatrix.txt {output_directory}/tmpSort.txt > {output_directory}/tmpMatrix.2.txt && \\
  mv {output_directory}/tmpMatrix.2.txt {output_directory}/tmpMatrix.txt && \\
  na=$(basename $read_count_file | cut -d. -f1) && \\
  HEAD="$HEAD\t$na"
done && \\
echo -e $HEAD | cat - {output_directory}/tmpMatrix.txt | tr ' ' '\t' > {output_matrix} && \\
rm {output_directory}/tmpSort.txt {output_directory}/tmpMatrix.txt""".format(
            reference_gtf=config.param('raw_counts_metrics', 'gtf', type='filepath'),
            output_directory=output_directory,
            read_count_files=" \\\n  ".join(read_count_files),
            output_matrix=output_matrix
        )
        jobs.append(job)

        # Create Wiggle tracks archive
        wiggle_directory = os.path.join("tracks", "bigWig")
        wiggle_archive = "tracks.zip"
        big_wig_prefix = os.path.join("tracks", "bigWig", sample.name)
        if config.param('DEFAULT', 'strand_info') != 'fr-unstranded':
            wiggle_files = []
            for sample in self.samples:
                wiggle_files.extend([os.path.join(wiggle_directory, sample.name) + ".forward.bw", os.path.join(wiggle_directory, sample.name) + ".reverse.bw"])
        else:
            wiggle_files = [os.path.join(wiggle_directory, sample.name + ".bw") for sample in self.samples]
        jobs.append(Job(wiggle_files, [wiggle_archive], name="metrics.wigzip", command="zip -r " + wiggle_archive + " " + wiggle_directory))

        # RPKM and Saturation
        count_file = os.path.join("DGE", "rawCountMatrix.csv")
        gene_size_file = config.param('saturation_rpkm', 'gene_size', type='filepath')
        rpkm_directory = "raw_counts"
        saturation_directory = os.path.join("metrics", "saturation")

        job = concat_jobs([
            Job(command="mkdir -p " + saturation_directory),
            metrics.rpkm_saturation(count_file, gene_size_file, rpkm_directory, saturation_directory)
        ], name="saturation_rpkm")
        jobs.append(job)

        return jobs

    def cufflinks(self):
        jobs = []
        
        gtf = config.param('cufflinks','gtf', type='filepath')

        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
            output_directory = os.path.join("cufflinks", sample.name)

            # De Novo FPKM
            job = cufflinks.cufflinks(input_bam, output_directory, gtf)
            job.name = "cufflinks."+sample.name
            jobs.append(job)

        return jobs
    
    def cuffmerge(self):
        
        output_directory = os.path.join("cufflinks", "AllSample")
        sample_file = os.path.join("cufflinks", "cuffmerge.samples.txt")
        input_gtfs = [os.path.join("cufflinks", sample.name, "transcripts.gtf") for sample in self.samples]
        gtf = config.param('cuffmerge','gtf', type='filepath')
        
        
        job = concat_jobs([
            Job(command="mkdir -p " + output_directory),
            Job(input_gtfs, [sample_file], command="""\
`cat > {sample_file} << END
{sample_rows}
END
  
`""".format(sample_rows="\n".join(input_gtfs), sample_file=sample_file)),
            cufflinks.cuffmerge(sample_file, output_directory, gtf_file=gtf)],
            name="cuffmerge")
        
        return [job]
        
    def cuffquant(self):
        jobs = []
        
        gtf = os.path.join("cufflinks", "AllSample","merged.gtf")
        
        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
            output_directory = os.path.join("cufflinks", sample.name)

            #Quantification
            job = cufflinks.cuffquant(input_bam, output_directory, gtf)
            job.name = "cuffquant."+sample.name
            jobs.append(job)

        return jobs
    
    def cuffdiff(self):
        jobs = []

        fpkm_directory = "cufflinks"
        gtf = os.path.join(fpkm_directory, "AllSample","merged.gtf")


        # Perform cuffdiff on each design contrast
        for contrast in self.contrasts:
            job = cufflinks.cuffdiff(
                # Cuffdiff input is a list of lists of replicate bams per control and per treatment
                [[os.path.join(fpkm_directory, sample.name, "abundances.cxb") for sample in group] for group in contrast.controls, contrast.treatments],
                gtf,
                os.path.join("cuffdiff", contrast.name)
            )
            job.name = "cuffdiff." + contrast.name
            jobs.append(job)

        return jobs
    
    def cuffnorm(self):
        jobs = []

        fpkm_directory = "cufflinks"
        gtf = os.path.join(fpkm_directory, "AllSample","merged.gtf")
        sample_labels = ",".join([sample.name for sample in self.samples])


        # Perform cuffnorm using every samples
        job = cufflinks.cuffnorm([os.path.join(fpkm_directory, sample.name, "abundances.cxb") for sample in self.samples],
             gtf,
             "cuffnorm",sample_labels)
        job.name = "cuffnorm" 
        jobs.append(job)
        
        return jobs

    def differential_expression(self):
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)
        output_directory = "DGE"
        count_matrix = os.path.join(output_directory, "rawCountMatrix.csv")

        edger_job = differential_expression.edger(design_file, count_matrix, output_directory)
        edger_job.output_files = [os.path.join(output_directory, contrast.name, "edger_results.csv") for contrast in self.contrasts]

        deseq_job = differential_expression.deseq(design_file, count_matrix, output_directory)
        deseq_job.output_files = [os.path.join(output_directory, contrast.name, "dge_results.csv") for contrast in self.contrasts]

        return [concat_jobs([
            Job(command="mkdir -p " + output_directory),
            edger_job,
            deseq_job
        ], name="differential_expression")]

    def differential_expression_goseq(self):
        jobs = []

        for contrast in self.contrasts:
            # goseq for differential gene expression results
            job = differential_expression.goseq(
                os.path.join("DGE", contrast.name, "dge_results.csv"),
                config.param("differential_expression_goseq", "dge_input_columns"),
                os.path.join("DGE", contrast.name, "gene_ontology_results.csv")
            )
            job.name = "differential_expression_goseq.dge." + contrast.name
            jobs.append(job)

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq(self):
        sample_fpkm_readcounts = [[
            sample.name,
            os.path.join("cufflinks", sample.name, "isoforms.fpkm_tracking"),
            os.path.join("raw_counts", sample.name + ".readcounts.csv")
        ] for sample in self.samples]

        input_file = os.path.join("exploratory", "exploratory.samples.tsv")

        return [concat_jobs([
            Job(command="mkdir -p exploratory"),
            gq_seq_utils.exploratory_analysis_rnaseq(
                os.path.join("DGE", "rawCountMatrix.csv"),
                "cuffnorm",
                config.param('gq_seq_utils_exploratory_analysis_rnaseq', 'genes', type='filepath'),
                "exploratory"
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq")]

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.star,
            self.picard_merge_sam_files,
            self.picard_sort_sam,
            self.picard_mark_duplicates,
            self.rnaseqc,
            self.wiggle,
            self.raw_counts,
            self.raw_counts_metrics,
            self.cufflinks,
            self.cuffmerge,
            self.cuffquant,
            self.cuffdiff,
            self.cuffnorm,
            self.differential_expression,
            self.differential_expression_goseq,
            self.gq_seq_utils_exploratory_analysis_rnaseq
        ]

if __name__ == '__main__':
    RnaSeq()
