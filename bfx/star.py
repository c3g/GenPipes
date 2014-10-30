#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *
from utils import utils 

def align(
    reads1,
    reads2,
    output_directory,
    genome_index_folder=config.param('star_align', 'genome_index_folder', required=True, type='dirpath'),
    rg_id="",
    rg_sample="",
    rg_library="",
    rg_platform_unit="",
    rg_platform="",
    rg_center="",
    create_wiggle_track=False,
    search_chimeres=False,
    cuff_follow=False
    ):

    job = Job(
        [reads1, reads2],
        [os.path.join(output_directory, "Aligned.out.bam"),os.path.join(output_directory, "SJ.out.tab")],
        [['star_align', 'module_star']]
    )
    
    ## get param form config filepat
    num_threads=config.param('star_align', 'threads')
    ram_limit=config.param('star_align', 'ram')
    max_ram=int(utils.number_symbol_converter(ram_limit))
    io_limit=config.param('star_align', 'io_buffer')
    io_max=int(utils.number_symbol_converter(io_limit))
    stranded=config.param('star_align','strand_info')
    wig_prefix=config.param('star_align','wig_prefix')
    chimere_segment_min=config.param('star_align','chimere_segment_min',type='int')
    
    ##wiggle information
    if create_wiggle_track:
        wig_cmd="--outWigType wiggle read1 5p"
        if stranded.lower() == "stranded":
            wig_cmd=wig_cmd+" --outWigStrand Stranded"
        elif stranded.lower() == "unstranded":
            wig_cmd=wig_cmd+" --outWigStrand Unstranded"
        else:
            raise Exception("Stand info\"" + stranded + "\" unrecongized")
        if wig_prefix != "" :
            wig_cmd=wig_cmd+" --outWigReferencesPrefix "+str(wig_prefix)
    else :
        wig_cmd=""
    
    ##chimeric information
    if search_chimeres and chimere_segment_min !+ "" :
        chim_cmd="--chimSegmentMin "+str(chimere_segment_min)
    else :
        chim_cmd=""
    
    ##check strandness if cufflinks will be run afterwards 
    if cuff_follow :
        if stranded.lower() == "stranded":
           cuff_cmd=""
        elif stranded.lower() == "unstranded":
            cuff_cmd="--outSAMstrandField intronMotif"
        else:
            raise Exception("Stand info\"" + stranded + "\" unrecongized")
    else :
        cuff_cmd=""
    
    job.command = """\
mkdir -p {output_directory} && \\
STAR --runMode alignReads \\
  --genomeDir {genome_index_folder} \\
  --readFilesIn {reads1}{reads2}\\
  --runThreadN {num_threads} \\
  --readFilesCommand zcat \\
  --outStd Log \\
  --outSAMunmapped Within \\
  --outSAMtype BAM Unsorted \\
  --outFileNamePrefix {output_directory} \\ 
  --limitBAMsortRAM {ram} \\
  --limitGenomeGenerateRAM {ram} \\
  {io_limit_size} \\
  {wig_param} \\
  {chim_param} \\
  --outSAMattrRGline \"ID:{rg_id} PL:{rg_platform} PU:{rg_platform_unit} LB:{rg_library} SM:{rg_sample} CN:{rg_center}\" \\
  {other_options}""".format(
        genome_index_folder=genome_index_folder
        other_options=config.param('star_align', 'other_options', required=False),
        rg_id=rg_id,
        rg_sample=rg_sample,
        rg_library=rg_library,
        rg_platform_unit=rg_platform_unit,
        rg_platform=rg_platform,
        rg_center=rg_center,
        output_directory=output_directory,
        num_threads=num_threads if str(num_threads) != "" and isinstance(num_threads, int) and  num_threads > 0 else 1,  
        ram=int(max_ram/2),
        reads1=reads1,
        reads2=" \\\n  " + reads2 if reads2 else "",
        io_limit_size="\\\n  --limitIObufferSize " + io_max if io_max else "",
        wig_param=wig_cmd,
        chim_param=chim_cmd
    )

    return job


def index(
    genome_index_folder,
    junction_file,
    gtf=config.param('star_align', 'gtf', required=False, type='filepath'),
    ):
    #STAR --runMode genomeGenerate --genomeDir $odir --genomeFastaFiles $genome --runThreadN $runThreadN --limitGenomeGenerateRAM $limitGenomeGenerateRAM --sjdbOverhang $sjdbOverhang  --sjdbFileChrStartEnd "

     job = Job(
        [junction_file],
        [os.path.join(genome_index_folder, "SAindex")],
        [['star_index', 'module_star']]
    )
     
    ## get param form config filepat
    reference_fasta=config.param('star_index', 'threads', required=True)
    num_threads=config.param('star_index', 'threads')
    ram_limit=config.param('star_index', 'ram')
    max_ram=int(utils.number_symbol_converter(ram_limit))
    io_limit=config.param('star_index', 'io_buffer')
    io_max=int(utils.number_symbol_converter(io_limit))
    read_size=config.param('star_index', 'cycle_number',type='int'))
     
    job.command = """\
mkdir -p {genome_index_folder} && \\
STAR --runMode genomeGenerate \\
  --genomeDir {genome_index_folder} \\
  --genomeFastaFiles {reference_fasta} \\
  --runThreadN {num_threads} \\
  --limitGenomeGenerateRAM  {ram} \\
  --sjdbFileChrStartEnd {junction_file} \\
  {gtf} \\
  {io_limit_size} \\
  {sjdbOverhang} \\
  {other_options}""".format(
        genome_index_folder=genome_index_folder
        other_options=config.param('star_index', 'other_options', required=False),
        num_threads=num_threads if str(num_threads) != "" and isinstance(num_threads, int) and  num_threads > 0 else 1,  
        ram=int(max_ram/2),
        reference_fasta=reference_fasta,
        junction_file=junction_file,
        gtf=" \\\n  --sjdbGTFfile " + gtf if gtf else "",
        sjdbOverhang="\\\n  --sjdbOverhang " + (read_size - 1) if str(read_size) != "" and isinstance(read_size, int) and  read_size > 0 else ""
        io_limit_size="\\\n  --limitIObufferSize " + io_max if io_max else ""
    )
    
    return job

def concatenate_junction(
    input_junction_files_list,
    output_junction_file
    ):
    
    
    job = Job(
        input_junction_files_list,
        [output_junction_file],
        [['star_junction', 'module_star']]
    )
    
    
    job.command = """\
cat {file_list} | awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' |sort -k1,1h -k2,2n > {output_junction_file}""".format(
    file_list=" ".join(input_junction_files_list),
    output_junction_file=output_junction_file
    )
    
    return job
