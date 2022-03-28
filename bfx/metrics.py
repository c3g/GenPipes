################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def dna_sample_metrics(input_directory, output, experiment_type="unknown"):
    return Job(
        [input_directory],
        [output],
        [
            ['dna_sample_metrics', 'module_R'],
            ['dna_sample_metrics', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/DNAsampleMetrics.R \\
  {input_directory} \\
  {output} \\
  {experiment_type}""".format(
            input_directory=input_directory,
            output=output,
            experiment_type=experiment_type
        )
    )

def rnaseqc(sample_file, output_directory, is_single_end=False, gtf_file=None, reference=None, ribosomal_interval_file=None):
    return Job(
        [sample_file],
        [os.path.join(output_directory, "index.html"), os.path.join(output_directory, "metrics.tsv"), os.path.join(output_directory, "corrMatrixSpearman.txt")],
        [
            ['rnaseqc', 'module_java'],
            ['rnaseqc', 'module_bwa'],
            ['rnaseqc', 'module_rnaseqc']
        ],
        command="""\
touch dummy_rRNA.fa && \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $RNASEQC_JAR \\
  -n {number_top_transcripts} \\
  -o {output_directory} \\
  -r {reference_genome_fasta} \\
  -s {sample_file} \\
  -t {gtf_file}{other_options}{single_end}{ribosomal_interval_file}""".format(
            tmp_dir=config.param('rnaseqc', 'tmp_dir'),
            java_other_options=config.param('rnaseqc', 'java_other_options'),
            ram=config.param('rnaseqc', 'ram'),
            number_top_transcripts=config.param('rnaseqc', 'number_top_transcript', param_type='int'),
            output_directory=output_directory,
            reference_genome_fasta=reference if reference else config.param('rnaseqc', 'genome_fasta', param_type='filepath'),
            sample_file=sample_file,
            gtf_file=gtf_file if gtf_file else config.param('rnaseqc', 'gtf', param_type='filepath'),
            other_options=(" \\\n  " + config.param('rnaseqc', 'other_options', required=False)) if config.param('rnaseqc', 'other_options', required=False) else "",
            single_end=" \\\n  -singleEnd" if is_single_end else "",
            ribosomal_interval_file=" \\\n  -rRNA " + ribosomal_interval_file if ribosomal_interval_file else "\\\n  -BWArRNA dummy_rRNA.fa"
        ),
        removable_files=["dummy_rRNA.fa"]
    )

def rpkm_saturation(count_file, gene_size_file, rpkm_directory, saturation_directory):
    return Job(
        [count_file],
        [saturation_directory + ".zip"],
        [
            ['rpkm_saturation', 'module_R'],
            ['rpkm_saturation', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/rpkmSaturation.R \\
  {count_file} \\
  {gene_size_file} \\
  {rpkm_directory} \\
  {saturation_directory} \\
  {threads} \\
  {other_options} && \\
zip -r {saturation_directory}.zip {saturation_directory}""".format(
            count_file=count_file,
            gene_size_file=gene_size_file,
            rpkm_directory=rpkm_directory,
            saturation_directory=saturation_directory,
            threads=config.param('rpkm_saturation', 'threads', param_type='posint'),
            other_options=config.param('rpkm_saturation', 'other_options', required=False)
        ),
        removable_files=[saturation_directory]
    )

def snv_graph_metrics(list, output_basename):
    return Job(
        [list],
        [output_basename + ".snvGraphMetrics_listFiles.txt"],
        [
            ['snv_graph_metrics', 'module_R'],
            ['snv_graph_metrics', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/snvGraphMetrics.R \\
  {list} \\
  {output_basename}""".format(
            list=list,
            output_basename=output_basename
        )
    )

def vcf_stats(input, output, list):
    return Job(
        [input],
        [output, list],
        [
            ['vcf_stats', 'module_python'],
            ['vcf_stats', 'module_mugqic_tools']
        ],
        command="""\
python $PYTHON_TOOLS/vcfStats.py \\
  -v {input} \\
  -d {dictionary} \\
  -o {output} \\
  -f {list}""".format(
            input=input,
            dictionary=config.param('vcf_stats', 'genome_dictionary', param_type='filepath'),
            output=output,
            list=list
        )
    )

def gc_bias(input, output):
    return Job(
        [input],
        [output],
        [
            ['gc_bias', 'module_R'],
            ['gc_bias', 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/GCbias_all.R \\
  {input} > {output}""".format(
            input=input,
            output=output
        )
    )

def ihec_metrics_rnaseq(genome):
  ''' Outputs the ihec metrics file for all samples'''

  ## will parse metrics/rnaseqRep/metrics.tsv to output needed columns only

  command = "python $PYTHON_TOOLS/ihec_metrics_rnaseq.py {genome}".format(genome=genome)

  return Job(input_files=["metrics/rnaseqRep/metrics.tsv", "report/trimAlignmentTable.tsv"],
             output_files=["report/IHEC_metrics_rnaseq_All.txt"],
             module_entries=[["ihec_metrics_rnaseq", "module_mugqic_tools"],
                             ["ihec_metrics_rnaseq", "module_samtools"],
                             ["ihec_metrics_rnaseq", "module_python"]],
             name="ihec_metrics_rnaseq",
             command=command
             )



def target_cpg_profile(
    input,
    output,
    sample,
    coverage_list = [1,10,15,20,25,30]
    ):
    
    
    return Job(
        [input],
        [output],
        [],
        command="""\
echo -e 'Sample\\t1\\t10\\t15\\t20\\t25\\t30' > {output} &&  \\
out={sample} &&  \\
for readcov in {cov_value} ; do cov=$(sed 1d {input} | awk '$11>'$readcov'' |wc -l) ; out="$out\\t$cov" ; done  &&  \\
echo -e "$out" >> {output} """.format(
            input=input,
            output=output,
            sample=sample,
            cov_value=" ".join(str(x) for x in coverage_list)
        )
    )


   

