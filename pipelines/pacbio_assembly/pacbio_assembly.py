#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from bio.readset import *

from bio import pacbio_tools
from bio import smrtanalysis
from pipelines import common

log = logging.getLogger(__name__)

class PacBioAssembly(common.MUGQICPipeline):

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_pacbio_readset_file(self.args.readsets.name)
        return self._readsets

    """
    Filtering. This step uses smrtpipe.py (From the SmrtAnalysis package) and will filter reads and subreads based on their length and QVs.
    1- fofnToSmrtpipeInput.py. 
    2- modify RS_Filtering.xml files according to reads filtering values entered in .ini file.
    3- smrtpipe.py with filtering protocol
    4- prinseq-lite.pl to write fasta file based on fastq file.
    Informative run metrics such as loading efficiency, readlengths and base quality are generated in this step as well.
    """
    def smrtanalysis_filtering(self):
        jobs = []

        for sample in self.samples:
            fofn = os.path.join("fofns", sample.name + ".fofn")
            bax_files = [bax_file for readset in sample.readsets for bax_file in readset.bax_files]
            filtering_directory = os.path.join(sample.name, "filtering")

            jobs.append(concat_jobs([
                Job([], [config.param('smrtanalysis_filtering', 'filtering_settings')], command="cp -a -f " + os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "protocols") + " ."),
                Job(command="mkdir -p fofns"),
                Job(bax_files, [fofn], command="""\
`cat > {fofn} << END
{bax_files}
END
`""".format(bax_files="\n".join(bax_files), fofn=fofn)),
                Job(command="mkdir -p " + filtering_directory),
                smrtanalysis.filtering(
                    fofn,
                    os.path.join(filtering_directory, "input.xml"),
                    os.path.join(sample.name, "filtering.xml"),
                    filtering_directory,
                    os.path.join(filtering_directory, "smrtpipe.log")
                )
            ], name="smrtanalysis_filtering." + sample.name))

        return jobs

    """
    Cutoff value for splitting long reads from short reads is done here using 
    estimated coverage and estimated genome size.
    
    You should estimate the overall coverage and length distribution for putting in
    the correct options in the configuration file. You will need to decide a
    length cutoff for the seeding reads. The optimum cutoff length will depend on
    the distribution of the sequencing read lengths, the genome size and the
    overall yield. Here, you provide a percentage value that corresponds to the 
    fraction of coverage you want to use as seeding reads.
    
    First, loop through fasta sequences. put the length of each sequence in an array, 
    sort it, loop through it again and compute the cummulative length coveredby each 
    sequence as we loop though the array. Once that length is > (coverage * genome 
    size) * $percentageCutoff (e.g. 0.10), we have our threshold. The idea is to 
    consider all reads above that threshold to be seeding reads to which will be 
    align lower shorter subreads.
    """
    def pacbio_tools_get_cutoff(self):
        jobs = []

        for sample in self.samples:
            log.info("Sample: " + sample.name)
            sample_nb_base_pairs = sum([readset.nb_base_pairs for readset in sample.readsets])
            log.info("nb_base_pairs: " + str(sample_nb_base_pairs))
            estimated_genome_size = sample.readsets[0].estimated_genome_size
            log.info("estimated_genome_size: " + str(estimated_genome_size))
            estimated_coverage = sample_nb_base_pairs / estimated_genome_size
            log.info("estimated_coverage: " + str(estimated_coverage))

            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                suffix = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, suffix)

                log.info("COVERAGE_CUTOFF: " + coverage_cutoff + "_X_coverage")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + os.path.join(coverage_directory, "preassembly")),
                    pacbio_tools.get_cutoff(
                        os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                        estimated_coverage,
                        estimated_genome_size,
                        coverage_cutoff,
                        os.path.join(coverage_directory, "preassemblyMinReadSize.txt")
                    )
                ], name="pacbio_tools_get_cutoff." + sample.name + ".coverage_cutoff" + suffix))

        return jobs

    """
    Having in hand a cutoff value, filtered reads are splitted between short and long reads. Short reads
    are aligned against long reads and consensus (e.g. corrected reads) are generated from these alignments.
    1- split reads between long and short.
    2- blasr (Aligner for PacBio reads)
    3- m4topre (Converts .m4 blasr output in .pre format.)
    4- pbdagcon (generates corrected reads from alignments)
    """
    def preassembly(self):
        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                suffix = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, suffix)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")
                job_name_suffix = sample.name + ".coverage_cutoff" + suffix
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + preassembly_directory),
                    pacbio_tools.split_reads(
                        os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                        os.path.join(coverage_directory, "preassemblyMinReadSize.txt"),
                        os.path.join(preassembly_directory, "filtered_shortreads.fa"),
                        os.path.join(preassembly_directory, "filtered_longreads.fa")
                    )
                ], name="pacbio_tools_split_reads." + job_name_suffix))

                job = smrtanalysis.blasr(
                    os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                    os.path.join(preassembly_directory, "filtered_longreads.fa"),
                    os.path.join(preassembly_directory, "seeds.m4"),
                    os.path.join(preassembly_directory, "seeds.m4.fofn")
                )
                job.name = "smrtanalysis_blasr." + job_name_suffix
                jobs.append(job)

                job = smrtanalysis.m4topre(
                    os.path.join(preassembly_directory, "seeds.m4.filtered"),
                    os.path.join(preassembly_directory, "seeds.m4.fofn"),
                    os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                    os.path.join(preassembly_directory, "aln.pre")
                )
                job.name = "smrtanalysis_m4topre." + job_name_suffix
                jobs.append(job)

                job = smrtanalysis.pbdagcon(
                    os.path.join(preassembly_directory, "aln.pre"),
                    os.path.join(preassembly_directory, "corrected.fasta"),
                    os.path.join(preassembly_directory, "corrected.fastq")
                )
                job.name = "smrtanalysis_pbdagcon." + job_name_suffix
                jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.smrtanalysis_filtering,
            self.pacbio_tools_get_cutoff,
            self.preassembly
        ]

if __name__ == '__main__': 
    PacBioAssembly().submit_jobs()
