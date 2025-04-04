################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def batch(
        tumor_bam,
        normal_bam,
        outdir,
        tar_dep=[],
        antitar_dep=[],
        target_bed=None,
        reference=None,
        output_cnn=None
):
    if tumor_bam is not None:
        inputs = [tumor_bam, normal_bam]
    else:
        inputs = [normal_bam]

    return Job(
        inputs,
        [tar_dep, antitar_dep],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py batch {options} \\
  -p {threads} \\
  --fasta {genome} \\
  --access {access} \\
  --annotate {annotate} \\
  {reference} \\
  {target_bed} \\
  {output_cnn} \\
  --output-dir {outdir} \\
  --normal {normal_bam} \\
  {tumor_bam}""".format(
        options=global_conf.global_get('cnvkit_batch', 'batch_options'),
        threads=global_conf.global_get('cnvkit_batch', 'threads'),
        genome=global_conf.global_get('cnvkit_batch', 'genome_fasta', param_type='filepath'),
        access=global_conf.global_get('cnvkit_batch', 'access', param_type='filepath'),
        annotate=global_conf.global_get('cnvkit_batch', 'refFlat', param_type='filepath'),
        reference="--reference " + reference if reference else "",
        target_bed="--targets " + target_bed if target_bed else "",
        output_cnn="--output-reference " + output_cnn if output_cnn else "",
        outdir=outdir,
        normal_bam=normal_bam if normal_bam else "",
        tumor_bam=tumor_bam if tumor_bam else "",
        )
    )

def fix(target_cov,
        antitarget_cov,
        output_cnr,
        reference=None,
        ref_cnn=None):
    return Job(
        [target_cov, antitarget_cov],
        [output_cnr],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py fix {options} \\
  {target_cov} \\
  {antitarget_cov} \\
  {reference} \\
  {ref_cnn} \\
  -o {output_cnr}""".format(
        options=global_conf.global_get('cnvkit_batch', 'fix_options'),
        target_cov=target_cov,
        antitarget_cov=antitarget_cov,
        reference=reference if reference else "",
        ref_cnn=ref_cnn if ref_cnn else "",
        output_cnr=output_cnr,
        )
    )

def segment(input_cnr,
            output_cns,
            vcf=None,
            sample_id=None,
            normal_id=None):

    if vcf:
        inputs = [input_cnr, vcf]
    else:
        inputs = [input_cnr]

    return Job(
        inputs,
        [output_cns],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py segment {options} \\
  {input_cnr} {vcf} {sample_id} \\
  -o {output_cns}""".format(
        options=global_conf.global_get('cnvkit_batch', 'segment_options'),
        input_cnr=input_cnr,
        output_cns=output_cns,
        vcf="--vcf " + vcf if vcf else "",
        sample_id="--sample-id " + sample_id if sample_id else "",
        normal_id="--normal-id " + normal_id if normal_id else "",
        )
    )

def call(input_cns,
         output_cns):
    return Job(
        [input_cns],
        [output_cns],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py call {options} \\
  {input_cns} \\
  -o {output_cns}""".format(
        options=global_conf.global_get('cnvkit_batch', 'call_options'),
        input_cns=input_cns,
        output_cns=output_cns,
        )
    )

def export(tumor_cns,
           output,
           sample_id=None):
    return Job(
        [tumor_cns],
        [output],
        [
            ['cnvkit_batch', 'module_cnvkit'],
        ],
        command="""\
cnvkit.py export {options} \\
  {sample_id} \\
  {output} \\
  {tumor_cns}""".format(
        options=global_conf.global_get('cnvkit_batch', 'export_options'),
        sample_id="-i " + '"' + sample_id + '"' if sample_id else "",
        output="-o " + output if output else "" ,
        tumor_cns=tumor_cns,
        )
    )

def metrics(input_cnr,
            input_cns,
            output):
    return Job(
        [input_cnr, input_cns],
        [output],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py metrics {options} \\
  {input_cnr} \\
  -s {input_cns} \\
  -o {output}""".format(
        options=global_conf.global_get('cnvkit_batch', 'metrics_options'),
        input_cnr=input_cnr,
        input_cns=input_cns,
        output=output,
        )
    )

def segmetrics(input_cnr,
               input_cns,
               output):
    return Job(
        [input_cnr, input_cns],
        [output],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py segmetrics {options} \\
  {input_cnr} \\
  -s {input_cns} \\
  -o {output}""".format(
        options=global_conf.global_get('cnvkit_batch', 'segmetrics_options'),
        input_cnr=input_cnr,
        input_cns=input_cns,
        output=output,
        )
    )

def select_background(input_cnr,
                      input_cns,
                      output):

    return Job(
        [input_cnr, input_cns],
        output,
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
	cnvkit.py metrics {options} \\
	  {input_cnr} \\
	  -s {input_cns} \\
	  -o {output}""".format(
            options=global_conf.global_get('cnvkit_batch', 'metrics_options'),
            input_cnr=input_cnr,
            input_cns=input_cns,
            output=output,
        )
    )

def read_metrics_file(
        in_file
        ):
    with open(in_file) as in_handle:
        header = next(in_handle).strip().split("\t")[1:]
        vals = map(float, next(in_handle).strip().split("\t")[1:])
    return dict(zip(header, vals))


def scatter(input_cnr,
            input_cns,
            output,
            vcf=None,
            normal=None,
            tumor=None):
    return Job(
        [input_cnr, input_cns, vcf],
        [output],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py scatter {options} \\
  {input_cnr} \\
  -s {input_cns} \\
  -o {output} {vcf} {tumor} {normal}""".format(
        options=global_conf.global_get('cnvkit_batch', 'scatter_options'),
        input_cnr=input_cnr,
        input_cns=input_cns,
        output=output,
        vcf="-v " + vcf if vcf else "",
        tumor="-i " + tumor if tumor else "",
        normal="-n " + normal if normal else "",
        )
    )

def diagram(input_cnr,
            input_cns,
            output):
    return Job(
        [input_cnr, input_cns],
        [output],
        [
            ['cnvkit_batch', 'module_cnvkit'],
            ['cnvkit_batch', 'module_R'],
        ],
        command="""\
cnvkit.py diagram {options} \\
  {input_cnr} \\
  -s {input_cns} \\
  -o {output}""".format(
        options=global_conf.global_get('cnvkit_batch', 'diagram_options'),
        input_cnr=input_cnr,
        input_cns=input_cns,
        output=output,
        )
    )

def file_check(
        input,
        output
):
    return Job(
        [input],
        [output],
        [
        ],
        command="""\
`line_count=$(wc -l < {input})
if [ "$line_count" -gt 1 ]; then
    touch {output}
fi`""".format(
            input=input,
            output=output,
        )
    )
