################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job
from ..utils import utils

log = logging.getLogger(__name__)

def align(
    reads1,
    reads2,
    output_directory,
    genome_index_folder,
    rg_id="",
    rg_sample="",
    rg_library="",
    rg_platform_unit="",
    rg_platform="",
    rg_center="",
    sort_bam=False,
    create_wiggle_track=False,
    search_chimeres=False,
    cuff_follow=False,
    two_pass=False,
    allsjdbFiles=None,
    ini_section='star_align'
    ):

    if not genome_index_folder:
        genome_index_folder = global_conf.global_get(ini_section, 'genome_index_folder', required=True).format(
            star_version=global_conf.global_get(ini_section, 'module_star').split('/')[-1]
        )
        if not os.path.exists(os.path.expandvars(genome_index_folder)):
            genome_index_folder = global_conf.global_get(ini_section, 'genome_index_folder', required=True).format(
            star_version=''
        )

    inputs = []
    inputs.extend([reads1,reads2])

    if two_pass:
        inputs.extend(allsjdbFiles)

    bam_name = "Aligned.sortedByCoord.out.bam" if sort_bam else "Aligned.out.bam"
    log_file = os.path.join(output_directory, "Log.final.out")

    job = Job(
        inputs,
        [os.path.join(output_directory, bam_name), os.path.join(output_directory, "SJ.out.tab"), log_file],
        [[ini_section, 'module_star']],
        removable_files=[os.path.join(output_directory, bam_name)]
    )
    ## Get param from config file
    num_threads = global_conf.global_get(ini_section, 'threads', param_type='int')
    ram_limit = global_conf.global_get(ini_section, 'ram')
    max_ram = int(utils.number_symbol_converter(ram_limit))
    io_limit = global_conf.global_get(ini_section, 'io_buffer')
    io_max = int(utils.number_symbol_converter(io_limit))
    stranded = global_conf.global_get(ini_section, 'strand_info')
    wig_prefix = global_conf.global_get(ini_section, 'wig_prefix')
    chimere_segment_min = global_conf.global_get(ini_section,'chimere_segment_min', param_type='int', required=False)
    ## Wiggle information
    if create_wiggle_track:
        wig_cmd = "--outWigType wiggle read1_5p"
        if stranded.lower() == "stranded":
            wig_cmd = wig_cmd + " --outWigStrand Stranded"
        elif stranded.lower() == "unstranded":
            wig_cmd = wig_cmd + " --outWigStrand Unstranded"
        else:
            raise Exception("Stand info\"" + stranded + "\" unrecognized")
        if wig_prefix != "":
            wig_cmd = wig_cmd + " --outWigReferencesPrefix " + str(wig_prefix)
    else:
        wig_cmd = ""

    ## Chimeric information
    if search_chimeres and chimere_segment_min != "":
        chim_cmd = "--chimSegmentMin " + str(chimere_segment_min) + " --chimJunctionOverhangMin " + str(chimere_segment_min) + " --chimSegmentReadGapMax 3"
    else:
        chim_cmd = ""

    ## Check strandness if cufflinks will be run afterwards
    if cuff_follow:
        if stranded.lower() == "stranded":
           cuff_cmd = ""
        elif stranded.lower() == "unstranded":
            cuff_cmd = "--outSAMstrandField intronMotif"
        else:
            raise Exception("Stand info\"" + stranded + "\" unrecognized")
    else:
        cuff_cmd = ""

    if two_pass:
        two_pass_cmd = "  --sjdbFileChrStartEnd" + " ".join(" " + sjdbFile for sjdbFile in allsjdbFiles),
    else:
        two_pass_cmd = ""

    other_options = global_conf.global_get(ini_section, 'other_options', required=False)

    job.command = """\
mkdir -p {output_directory} && \\
STAR --runMode alignReads \\
  --genomeDir {genome_index_folder} \\
  --readFilesIn \\
    {reads1}{reads2} \\
  --runThreadN {num_threads} \\
  --readFilesCommand zcat \\
  --outStd Log \\
  --outSAMunmapped Within \\
  --outSAMtype BAM {sort_value} \\
  --outFileNamePrefix {output_directory}/ \\
  --outSAMattrRGline {rg_id}\t{rg_platform}\t{rg_platform_unit}\t{rg_library}\t{rg_sample}\t{rg_center} \\
  --outTmpDir {tmp_dir}/$(mktemp -u star_XXXXXXXX) \\
  --limitGenomeGenerateRAM {ram}{sort_ram}{io_limit_size}{wig_param}{chim_param}{cuff_cmd}{other_options}""".format(
        output_directory=output_directory,
        genome_index_folder=genome_index_folder,
        reads1=reads1,
        reads2=" \\\n    " + reads2 if reads2 else "",
        num_threads=num_threads if str(num_threads) != "" and isinstance(num_threads, int) and num_threads > 0 else 1,
        sort_value="SortedByCoordinate" if sort_bam else "Unsorted",
        rg_id="ID:\"" + rg_id + "\" " if rg_id != "" else "",
        rg_platform="PL:\"" + rg_platform + "\" " if rg_platform != "" else "",
        rg_platform_unit=" PU:\"" + rg_platform_unit + "\" " if rg_platform_unit != "" else "",
        rg_library=" LB:\"" + rg_library + "\" " if rg_library != "" else "",
        rg_sample="SM:\"" + rg_sample + "\" " if rg_sample != "" else "",
        rg_center="CN:\"" + rg_center + "\" " if rg_center != "" else "",
        ram=int(max_ram),
        sort_ram=" \\\n  --limitBAMsortRAM " + str(int(max_ram)) if sort_bam else "",
        io_limit_size=" \\\n  --limitIObufferSize " + str(io_max) if io_max else "",
        wig_param=" \\\n  " + wig_cmd if wig_cmd else "",
        chim_param=" \\\n  " + chim_cmd if chim_cmd else "",
        cuff_cmd=" \\\n  " + cuff_cmd if cuff_cmd else "",
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir', required=True),
	    two_pass_cmd=" \\\n" + " ".join(two_pass_cmd) if two_pass_cmd else "",
        other_options=" \\\n  " + other_options if other_options else ""
    )

    return job


def index(
    genome_index_folder,
    junction_file,
    genome_length,
    gtf=global_conf.global_get('star_align', 'gtf', param_type='filepath', required=False),
    ini_section='star_index'
    ):

    job = Job(
        [junction_file],
        [os.path.join(genome_index_folder, "SAindex")],
        [[ini_section, 'module_star']],
        removable_files=[genome_index_folder]
    )

    ## get param from config filepath
    reference_fasta = global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath')
    num_threads = global_conf.global_get(ini_section, 'threads', param_type='int')
    ram_limit = global_conf.global_get(ini_section, 'ram')
    max_ram = int(utils.number_symbol_converter(ram_limit))
    io_limit = global_conf.global_get(ini_section, 'io_buffer')
    io_max = int(utils.number_symbol_converter(io_limit))
    read_size = global_conf.global_get(ini_section, 'star_cycle_number', param_type='posint')
    other_options = global_conf.global_get(ini_section, 'other_options', required=False)

    job.command = """\
mkdir -p {genome_index_folder} && \\
STAR --runMode genomeGenerate \\
  --genomeDir {genome_index_folder} \\
  --genomeFastaFiles {reference_fasta} \\
  --genomeSAindexNbases {genome_length} \\
  --runThreadN {num_threads} \\
  --limitGenomeGenerateRAM {ram} \\
  --outTmpDir {tmp_dir}/$(mktemp -u star_XXXXXXXX) \\
  --sjdbFileChrStartEnd {junction_file}{gtf}{io_limit_size}{sjdbOverhang}{other_options}""".format(
        genome_index_folder=genome_index_folder,
        genome_length=genome_length,
        reference_fasta=reference_fasta,
        num_threads=num_threads if str(num_threads) != "" and isinstance(num_threads, int) and  num_threads > 0 else 1,
        ram=max_ram,
        junction_file=junction_file,
        gtf=" \\\n  --sjdbGTFfile " + gtf if gtf else "",
        io_limit_size=" \\\n  --limitIObufferSize " + str(io_max) if io_max else "",
        sjdbOverhang=" \\\n  --sjdbOverhang " + str(read_size) if read_size else "",
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir', required=True),
        other_options=" \\\n  " + other_options if other_options else ""
    )

    return job

def concatenate_junction(
    input_junction_files_list,
    output_junction_file,
    ini_section='star_junction'
    ):

    job = Job(
        input_junction_files_list,
        [output_junction_file],
        [[ini_section, 'module_star']]
    )

    job.command = """\
cat \\
  {file_list} | \\
awk 'BEGIN {{OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"}} {{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' | sort -k1,1h -k2,2n > {output_junction_file}""".format(
    file_list=" \\\n  ".join(input_junction_files_list),
    output_junction_file=output_junction_file
    )

    return job
