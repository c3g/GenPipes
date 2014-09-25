#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def filtering(
    fofn,
    input_xml,
    params_xml,
    output_dir,
    log
    ):

    ref_params_xml=config.param('smrtanalysis_filtering', 'filtering_settings')
    output = os.path.join(output_dir, "data", "filtered_subreads.fastq")

    job = Job(
        [fofn, ref_params_xml],
        [output],
        [
            ['smrtanalysis_filtering', 'module_perl'],
            ['smrtanalysis_filtering', 'module_python'],
            ['smrtanalysis_filtering', 'module_memtime'],
            ['smrtanalysis_filtering', 'module_prinseq'],
            ['smrtanalysis_filtering', 'module_smrtanalysis']
        ]
    )

    job.command = """\
set +u && \\
source \$SEYMOUR_HOME/etc/setup.sh && \\
set -u && \\
memtime fofnToSmrtpipeInput.py {fofn} > {input_xml} && \\
cp {fofn} {output_dir}/input.fofn && \\
memtime sed -e \'s/MINSUBREADLENGTH/{min_subread_length}/g\' -e \'s/MINREADLENGTH/{min_read_length}/g\' -e \'s/MINQUAL/{min_qual}/g\' \\
  < {ref_params_xml} > {params_xml} && \\
memtime smrtpipe.py \\
  -D NPROC={threads} \\
  -D TMP={tmp_dir} \\
  --params={params_xml} \\
  --output={output_dir} \\
  --debug \\
  xml:{input_xml} \\
  > {log} && \\
memtime prinseq-lite.pl \\
  -verbose \\
  -fastq {output} \\
  -out_format 1 \\
  -out_good {output_dir}/data/filtered_subreads
""".format(
        fofn=fofn,
        input_xml=input_xml,
        min_subread_length=config.param('smrtanalysis_filtering', 'min_subread_length'),
        min_read_length=config.param('smrtanalysis_filtering', 'min_read_length'),
        min_qual=config.param('smrtanalysis_filtering', 'min_qual'),
        ref_params_xml=ref_params_xml,
        params_xml=params_xml,
        threads=config.param('smrtanalysis_filtering', 'threads'),
        tmp_dir=config.param('smrtanalysis_filtering', 'tmp_dir'),
        output_dir=output_dir,
        log=log,
        output=output
    )

    return job
