#!/usr/bin/env python

from core.job import *
from core.config import *
from bfx import bvatools
from bfx import bwa
from bfx import metrics
from bfx import picard
from bfx import star
from bfx import tools

class RunProcessingAligner(object):
    def __init__(self, output_dir):
        self._output_dir = output_dir

    @property
    def output_dir(self):
        return self._output_dir
    def get_reference_index(self, genome_folder):
        raise NotImplementedError("Please Implement this method")

    def get_alignment_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_metrics_jobs(self, readset):
        raise NotImplementedError("Please Implement this method")

    def get_annotation_files (self, genome_folder):
        raise NotImplementedError("Please Implement this method")


class BwaRunProcessingAligner(RunProcessingAligner):
    @property
    def created_interval_lists(self):
        if not hasattr(self, "_created_interval_lists"):
            self._created_interval_lists= []
        return self._created_interval_lists

    @property
    def downloaded_bed_files(self):
        if not hasattr(self, "_downloaded_bed_files"):
            self._downloaded_bed_files = []
        return self._downloaded_bed_files

    def get_reference_index(self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        return os.path.join(genome_folder,
                            "genome",
                            "bwa_index",
                            folder_name + ".fa")

    def get_annotation_files (self, genome_folder):
        return [];

    def get_alignment_jobs(self, readset):
        jobs = []
        output = readset.bam + ".bam"
        job = concat_jobs([
            Job(command="mkdir -p " + os.path.dirname(output)),
            pipe_jobs([
                bwa.mem(
                    readset.fastq1,
                    readset.fastq2,
                    read_group="'@RG" + \
                        "\tID:" + readset.library + "_" + readset.run + "_" + readset.lane + \
                        "\tSM:" + readset.sample.name + \
                        ("\tLB:" + readset.library if readset.library else "") + \
                        ("\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                        ("\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                        "\tPL:Illumina" + \
                        "'",
                    ref=readset.aligner_reference_index
                ),
                picard.sort_sam(
                    "/dev/stdin",
                    output,
                    "coordinate"
                )
            ])
        ], name="bwa_mem_picard_sort_sam." + readset.name)

        jobs.append(job)
        return jobs

    def get_metrics_jobs(self, readset):
        jobs = []

        input_file_prefix = readset.bam + '.'
        input = input_file_prefix + "bam"

        job = picard.collect_multiple_metrics(input, input_file_prefix + "metrics", reference_sequence=readset.reference_file)
        job.name = "picard_collect_multiple_metrics." + readset.name + ".met"
        jobs.append(job)

        coverage_bed = bvatools.resolve_readset_coverage_bed(readset)
        full_coverage_bed = (self.output_dir + os.sep + coverage_bed) if coverage_bed else None

        if coverage_bed:
            if (not os.path.exists(full_coverage_bed)) and (coverage_bed not in self.downloaded_bed_files):
                # Download the bed file
                command = config.param('DEFAULT', 'fetch_bed_file_command').format(
                        output_directory = self.output_dir,
                        filename = coverage_bed
                )
                job = Job([], [full_coverage_bed], command=command, name="bed_download." + coverage_bed)
                self.downloaded_bed_files.append(coverage_bed)
                jobs.append(job)

            interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

            if not interval_list in self.created_interval_lists:
                # Create one job to generate the interval list from the bed file
                ref_dict = os.path.splitext(readset.reference_file)[0] + '.dict'
                job = tools.bed2interval_list(ref_dict, full_coverage_bed, interval_list)
                job.name = "interval_list." + coverage_bed
                self.created_interval_lists.append(interval_list)
                jobs.append(job)

            job = picard.calculate_hs_metrics(input_file_prefix + "bam", input_file_prefix + "metrics.onTarget.txt", interval_list, reference_sequence=readset.reference_file)
            job.name = "picard_calculate_hs_metrics." + readset.name + ".hs"
            jobs.append(job)

        job = bvatools.depth_of_coverage(
                input, 
                input_file_prefix + "metrics.targetCoverage.txt", 
                full_coverage_bed, 
                other_options = config.param('bvatools_depth_of_coverage', 'other_options', required=False),
                reference_genome = readset.reference_file
        )
        job.name = "bvatools_depth_of_coverage." + readset.name + ".doc"
        jobs.append(job)

        return jobs


class StarRunProcessingAligner(RunProcessingAligner):
    def __init__(self, output_dir, nb_cycles):
        super(StarRunProcessingAligner, self).__init__(output_dir)
        self._nb_cycles = nb_cycles

    @property
    def nb_cycles(self):
        return self._nb_cycles

    def get_reference_index (self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        ini_file = os.path.join(genome_folder + os.sep + folder_name + ".ini")
        genome_config = ConfigParser.SafeConfigParser()
        genome_config.read(ini_file)

        source = genome_config.get("DEFAULT", "source")
        version = genome_config.get("DEFAULT", "version")

        return os.path.join(genome_folder,
                            "genome",
                            "star_index",
                            source + version + ".sjdbOverhang" + str(self.nb_cycles-1))

    def get_annotation_files (self, genome_folder):
        folder_name = os.path.basename(genome_folder)
        ini_file = os.path.join(genome_folder + os.sep + folder_name + ".ini")
        genome_config = ConfigParser.SafeConfigParser()
        genome_config.read(ini_file)

        source = genome_config.get("DEFAULT", "source")
        version = genome_config.get("DEFAULT", "version")

        return [
            os.path.join(genome_folder,
                        "annotations",
                        folder_name + '.' + source + version + ".transcript_id.gtf"),

            os.path.join(genome_folder,
                        "annotations",
                        "ncrna_bwa_index",
                        folder_name + '.' + source + version + ".ncrna.fa")
        ]

    def get_alignment_jobs(self, readset):
        jobs = []
        output = readset.bam + ".bam"

        rg_center = config.param('star_align', 'sequencing_center', required=False)
        star_bam_name = "Aligned.sortedByCoord.out.bam"

        job = concat_jobs([
            star.align(
                reads1=readset.fastq1,
                reads2=readset.fastq2,
                output_directory=os.path.dirname(output),
                sort_bam=True,
                genome_index_folder=readset.aligner_reference_index,
                rg_id=readset.library + "_" + readset.run + "_" + readset.lane,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform="Illumina",
                rg_center=rg_center if rg_center else ""
            ),
            Job(output_files=[output], command="mv " + os.path.dirname(output) + os.sep + star_bam_name + " " + output),
            Job(command="ln -s " + output + " " + os.path.dirname(output) + os.sep + star_bam_name),
            picard.build_bam_index(output, output[::-1].replace(".bam"[::-1], ".bai"[::-1], 1)[::-1])
        ])
        job.name = "star_align." + readset.name
        jobs.append(job)
        return jobs

    def get_metrics_jobs(self, readset):
        jobs = []
        input_bam = readset.bam + ".dup.bam"
        input_bam_directory = os.path.dirname(input_bam)
        sample_file = input_bam + ".sample_file"
        sample_row = readset.sample.name + "\t" + input_bam + "\tRNAseq"
        output_directory = os.path.join(input_bam_directory, "rnaseqc_" + readset.sample.name + "." + readset.library)


        if len(readset.annotation_files) > 1 and os.path.isfile(readset.annotation_files[0]) and os.path.isfile(readset.annotation_files[1]):
            gtf_transcript_id = readset.annotation_files[0]
            ribosomal_fasta = readset.annotation_files[1]
            reference= readset.reference_file
            job = concat_jobs([
                Job(command="mkdir -p " + output_directory),
                Job([input_bam], [sample_file], command="""\
    echo "Sample\tBamFile\tNote
    {sample_row}" \\
    > {sample_file}""".format(sample_row=sample_row, sample_file=sample_file)),
                metrics.rnaseqc(sample_file, output_directory, readset.fastq2 is not None, gtf_file=gtf_transcript_id, ribosomal_fasta=ribosomal_fasta, reference=reference),
                Job(command="cp " + os.path.join(output_directory, "metrics.tsv") + " " +  os.path.join(input_bam_directory, readset.sample.name + "." + readset.library + ".rnaseqc.sorted.dup.metrics.tsv")),
            ], name="rnaseqc" + readset.name + ".rnaseqc")
            jobs.append(job)

        return jobs
