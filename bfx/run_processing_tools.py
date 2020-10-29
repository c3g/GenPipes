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
import sys
import os
import re

# MUGQIC Modules
from core.config import config
from core.job import Job

def index(
    input,
    output,
    basecalls_dir,
    mismatches,
    lane,
    mask
    ):

    barcode_file = config.param('index', 'barcode_file', type='filepath', required='false')
    if not (barcode_file and os.path.isfile(barcode_file)):
        barcode_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), 'resources', 'barcode.mergedup.txt')

    return Job(
        [input],
        [output],
        [
            ["index", "module_java"],
            ["index", "module_mugqic_tools"]
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} \\
  {java_other_options} \\
  -Xmx{ram} \\
  -jar $JAVA_TOOLS/{jar} \\
  MAX_MISMATCHES={mismatches} \\
  NUM_PROCESSORS={threads} \\
  BARCODE_FILE={barcode_file} \\
  BASECALLS_DIR={basecalls_dir} \\
  LANE={lane_number} \\
  READ_STRUCTURE={read_structure} \\
  METRICS_FILE={output} \\
  TMP_DIR={tmp_dir}""".format(
            tmp_dir=config.param('index', 'tmp_dir'),
            java_other_options=config.param('index', 'java_other_options'),
            ram=config.param('index', 'ram'),
            jar=config.param('index', 'jar'),
            mismatches=mismatches,
            threads=config.param('index', 'threads'),
            barcode_file=barcode_file,
            basecalls_dir=basecalls_dir,
            lane_number=lane,
            read_structure=mask,
            output=output
        )
    )

def bcl2fastq(
    input,
    fastq_outputs,
    output_dir,
    sample_sheet,
    run,
    lane,
    extra_option,
    demultiplex=False,
    mismatches=None,
    mask=None,
    ):

    if demultiplex:
        demultiplex_parameters = """\
  --barcode-mismatches {number_of_mismatches} \\
  --use-bases-mask {mask}""".format(
            number_of_mismatches=mismatches,
            mask=mask
        )

    return Job(
        [input],
        fastq_outputs,
        [
            ['fastq', 'module_bcl_to_fastq']
        ],
        command="""\
bcl2fastq \\
  --runfolder-dir {run_dir} \\
  --output-dir {output_dir} \\
  --tiles {tiles} \\
  --sample-sheet {sample_sheet} \\
  --create-fastq-for-index-reads \\
  {demultiplex_parameters} \\
  {other_options} \\
  {extra_option}""".format(
            run_dir=run,
            output_dir=output_dir,
            tiles="s_" + str(lane),
            sample_sheet=sample_sheet,
            demultiplex_parameters=demultiplex_parameters if demultiplex_parameters else "",
            other_options=config.param('fastq', 'other_options'),
            extra_option=extra_option
        )
    )

def fastqmultx(
    barcode_file,
    mismatches,
    outputs,
    R1_fastq,
    I1_fastq,
    R2_fastq,
    I2_fastq=None
    ):

    metrics_file = outputs[0]
    output_dir = os.path.dirname(metrics_file)
    return Job(
        [
            R1_fastq,
            I1_fastq,
            R2_fastq,
            I2_fastq if I2_fastq else None
        ],
        outputs,
        [
            ['fastq', 'module_fastq_multx']
        ],
        command="""\
fastq-multx \\
  -B {barcode} \\
  -m {mismatches} \\
  {input_I1} \\
  {input_I2} \\
  {input_R1} \\
  {input_R2} \\
  {output_I1_fake} \\
  {output_I2_fake} \\
  {output_R1} \\
  {output_R2} \\
  > {stdout_metrics}""".format(
            barcode=barcode_file,
            mismatches=mismatches,
            input_I1=I1_fastq,
            input_I2=I2_fastq if I2_fastq else "",
            input_R1=R1_fastq,
            input_R2=R2_fastq,
            output_I1_fake="-o " + os.path.join(output_dir, "%_I1.fastq"),
            output_I2_fake="-o " + os.path.join(output_dir, "%_I2.fastq") if I2_fastq else "",
            output_R1="-o " + os.path.join(output_dir, "%_R1.fastq"),
            output_R2="-o " + os.path.join(output_dir, "%_R2.fastq"),
            stdout_metrics=metrics_file
        ),
        removable_files=outputs+[output_dir]
    )

def demux_fastqs(
    sample_sheet,
    mismatches,
    mask,
    outputs,
    metrics_file,
    R1_fastq,
    R2_fastq,
    ini_section='fastq'
    ):

    return Job(
        [
            R1_fastq,
            R2_fastq
        ],
        outputs + [metrics_file],
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_fgbio']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} \\
  {java_other_options} \\
  -Xmx{ram} \\
  -jar $FGBIO_JAR DemuxFastqs \\
  --threads {threads} \\
  --max-mismatches {mismatches} \\
  --metrics {metrics_file} \\
  --inputs {inputs} \\
  --read-structures {read_structure} \\
  --metadata {sample_sheet} \\
  --output {output_dir} \\
  --output-type fastq \\
  --include-all-bases-in-fastqs true""".format(
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            ram=config.param(ini_section, 'ram'),
            threads=config.param(ini_section, 'threads'),
            mismatches=mismatches,
            metrics_file=metrics_file,
            inputs=" ".join([R1_fastq, R2_fastq]),
            read_structure=mask,
            sample_sheet=sample_sheet,
            output_dir=os.path.dirname(metrics_file)
        ),
        removable_files=[os.path.dirname(metrics_file)]
    )

def demux_fastqs_by_chunk(
    sample_sheet,
    mismatches,
    mask,
    outputs,
    metrics_file,
    R1_fastq,
    R2_fastq,
    ini_section='fastq'
    ):

    metrics_folder = os.path.dirname(metrics_file)
    return Job(
        [
            R1_fastq,
            R2_fastq
        ],
        outputs + [metrics_file],
        [
            [ini_section, 'module_parallel'],
            [ini_section, 'module_java'],
            [ini_section, 'module_fgbio'],
            [ini_section, 'module_mugqic_tools']
        ],
        command="""\
FASTQ_1_IN={fastq1}
FASTQ_2_IN={fastq2}
READS_PER_CHUNK={reads_per_chunk}
GZIP_PARALLEL={gzip_thread}   # Rembering that we have $GZIP_PARALLEL processes for fwd and the same again for reverse
FGBIO_PARALLEL={fgbio_thread}  # Number of parallel instances of FGBIO DemuxFastqs
FGBIO_THREADS={threads}   # Number of threads per instance of FGBIO DemuxFastqs
let "BLOCKSIZE = 300 * $READS_PER_CHUNK"

rm -f {tmp_dir} && mkdir -p {tmp_dir}
rm -f {tmp_dir}/chunks_R[12] joblog.*.txt && mkfifo {tmp_dir}/chunks_R1 {tmp_dir}/chunks_R2
rm -rf {tmp_dir}/out && mkdir -p {tmp_dir}/out
rm -rf {tmp_dir}/split && mkdir -p {tmp_dir}/split

# Start chunking up the huge fastqgz into smaller pieces.
# Chunk filenames+paths are written into the named pipes chunks_R1 and chunks_R2.
echo "Chunking of the raw fastq files..."
zcat $FASTQ_1_IN | \\
  parallel \\
    --pipe \\
    --keep-order \\
    -j $GZIP_PARALLEL \\
    --joblog joblog.1.txt \\
    --blocksize $BLOCKSIZE \\
    -N $READS_PER_CHUNK \\
    -L 4 \\
    "gzip -c > {tmp_dir}/split/in.1.child_{{#}}.gz; echo {tmp_dir}/split/in.1.child_{{#}}.gz" > {tmp_dir}/chunks_R1 & \\
zcat $FASTQ_2_IN  \\
| parallel \\
    --pipe \\
    --keep-order \\
    -j $GZIP_PARALLEL \\
    --joblog joblog.2.txt \\
    --blocksize $BLOCKSIZE \\
    -N $READS_PER_CHUNK \\
    -L 4 \\
    "gzip -c > {tmp_dir}/split/in.2.child_{{#}}.gz; echo {tmp_dir}/split/in.2.child_{{#}}.gz" > {tmp_dir}/chunks_R2 & \\

# Take those names pipes and process each pair of chunked fastqgzs in parallel.
echo "Demultinplexing the fastq chunks..."
mkdir -p {metrics_folder}/chunk && \\
parallel \\
  -j $FGBIO_PARALLEL \\
  --joblog joblog.fgbio.txt \\
  "echo \"Demuxing chunks {{1}} {{2}}\" && \\
   java -Djava.io.tmpdir={tmp_dir} \\
    {java_other_options} \\
    -Xmx{ram} \\
    -jar $FGBIO_JAR DemuxFastqs \\
    --threads $FGBIO_THREADS \\
    --max-mismatches {mismatches} \\
    --metrics {metrics_folder}/chunk/metrics.{{#}}.txt \\
    --inputs {{1}} {{2}} \\
    --read-structures {read_structure} \\
    --metadata {sample_sheet} \\
    --output {tmp_dir}/out/out_{{#}} \\
    --output-type fastq \\
    --include-all-bases-in-fastqs true && \\
   echo \"Removing chunks {{1}} {{2}}\" && \\
   rm {{1}} {{2}}" \\
  :::: {tmp_dir}/chunks_R1 \\
  ::::+ {tmp_dir}/chunks_R2 && \\

# Combine fastqs (iterate though samples and combine fastqgzs for each sample)
echo "Combining fastq chunks" && \\
for fastqgz in `ls {tmp_dir}/out/out_1/*.gz`; do
  basename=$(basename $fastqgz)
  zcat $(find {tmp_dir}/out -name $basename) | gzip -c > {output_dir}/$basename &
done
wait

# Combine metrics using mugqic_tools
python $PYTHON_TOOLS/combineDemuxFastqsMetrics.py \\
  -i $(ls {metrics_folder}/chunk/metrics.*.txt)\\
  -o {metrics_file}""".format(
            tmp_dir=config.param(ini_section, 'tmp_dir'),
            gzip_thread=config.param(ini_section, 'gzip_parallel_thread'),
            fgbio_thread=config.param(ini_section, 'fgbio_parallel_thread'),
            reads_per_chunk=config.param(ini_section, 'reads_per_chunk'),
            java_other_options=config.param(ini_section, 'java_other_options'),
            ram=config.param(ini_section, 'ram'),
            threads=config.param(ini_section, 'threads'),
            mismatches=mismatches,
            metrics_folder=metrics_folder,
            metrics_file=metrics_file,
            fastq1=R1_fastq,
            fastq2=R2_fastq,	
            read_structure=mask,
            sample_sheet=sample_sheet,
            output_dir=os.path.dirname(metrics_file)
        ),
        removable_files=[os.path.dirname(metrics_file)]
    )

def bcl2fastq_for_index(
    run_dir,
    output_dir,
    sample_sheet,
    flowcell,
    lane,
    demultiplex=False,
    mismatches=None,
    mask=None,
    ):

    if demultiplex:
        demultiplex_parameters = """\
  --barcode-mismatches {number_of_mismatches} \\
  --use-bases-mask {mask}""".format(
            number_of_mismatches=mismatches,
            mask=mask
        )

    return Job(
        [run_dir],
        [
            os.path.join(output_dir, "Reports/html", flowcell, "all/all/all/lane.html"),
            os.path.join(output_dir, "Stats/Stats.json")
        ],
        [
            ['bcl2fastq_index', 'module_bcl_to_fastq']
        ],
        command="""\
bcl2fastq \\
  --runfolder-dir {run_dir} \\
  --output-dir {output_dir} \\
  --tiles {tiles} \\
  --sample-sheet {sample_sheet} \\
  --create-fastq-for-index-reads \\
  {demultiplex_parameters} \\
  {other_options}""".format(
            run_dir=run_dir,
            output_dir=output_dir,
            tiles="s_" + str(lane),
            sample_sheet=sample_sheet,
            demultiplex_parameters=demultiplex_parameters if demultiplex_parameters else "",
            other_options=config.param('bcl2fastq_index', 'other_options'),
        )
    )

def mgi_t7_basecall(
    input,
    flowcell_id,
    fastq_outputs,
    output_dir,
    json_flag_file,
    lane_config_file,
    ini_section='basecall'
    ):

    return Job(
        [input, json_flag_file],
        fastq_outputs,
        [
            [ini_section, 'module_basecall_t7'],
            [ini_section, 'module_python']
        ],
        command="""\
# Create environment variables to mathose in the config template file
export MGI_MEMORYSIZEMB={mem}
export MGI_THREADCOUNT={thread}
export MGI_WRITEFQMEMORY={mem}
export MGI_SAVEIMAGE_DIR=""
export MGI_WORKSPACE_DIR={workspace}
export MGI_GENFASTQ_DIR={fastq}

# Replace the environment variables by their values and create the config file to use
envsubst < {config_template} > {config_file}

# Create a sym link to the report scripts
BINARY=$(which processor)
ln -s -f ${{BINARY%/*}}/report {target_dir}/

cpulimit -i -l {thread}00 processor {config_file} & \\
sleep 10 && \\
client_linux \\
    ignored 1 1 1 1 \\
    -F \\
    -N {fcid} \\
    --writefq_config {json_flag} \\
    -I {config_file}""".format(
            mem=config.param(ini_section, 'mem'),
            thread=config.param(ini_section, 'thread'),
            workspace=config.param(ini_section, 'mgi_t7_workspace'),
            fastq=output_dir,
            config_template=os.path.expandvars(config.param(ini_section, 'mgi_t7_config_template')),
            config_file=lane_config_file,
            target_dir=os.path.dirname(lane_config_file),
            fcid=flowcell_id,
            json_flag=json_flag_file
        )
    )
