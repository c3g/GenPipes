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

from core.job import *
from core.config import *


def modify_inputinfofile(input_file, sample, histone,  output_inputinfo_file):

    return Job(
        [input_file, output_inputinfo_file],
        [output_inputinfo_file],
        [],
        command="""\
        echo -e "{sample}\\t{histone}\\t{sample}_{histone}.bedgraph.gz" >> {inputinfofile}""".format(
            input_files=input_file,
            sample=sample,
            histone=histone,
            inputinfofile=output_inputinfo_file

        )
    )

def generate_chr_sizes(chr_sizes_file, chr_sizes, chr):
    return Job(
        #to do: remove comment and add output_files. if got errors remove it.
        input_files=[chr_sizes],
        output_files=[chr_sizes_file],
      #  output_files=[chr_sizes_file],
        command="""\
        awk -v OFS="\\t" '{{  if($1=="{chr}") {{print $1,$2}}  }}' {chr_sizes} >> {chr_sizes_file}""".format(
            chr_sizes_file=chr_sizes_file,
            chr_sizes=chr_sizes,
            chr=chr )
    )

def convert_chr_bedgraph(input_file, output_file, chr, outputdir):
    return Job(
        [input_file],
        [output_file],
        [],
        command="""\
          zless {input_file} | awk '{{ if($1=="{chr}") {{print $0}}  }}' | gzip -cf > {output_file}""".format(
            input_files=input_file,
            input_file=input_file,
            output_file=output_file,
            chr=chr
        )

    )

def convert(input_dir, input_file, output_dir, output_files, inputinfofile, histone_mark, sample, chr_sizes_file):
    # input = input_files.extend(inputinfofile)

    return Job(
        [inputinfofile, input_file, chr_sizes_file],
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_convert." + sample + "." + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Convert \\
  -m {histone_mark} \\
  -l {convertsample} \\
  -r {resolution} \\
  {path_to_dataset} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
      java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
      ram=global_config_parser.param('chromimpute', 'ram'),
      histone_mark=histone_mark,
      convertsample=sample,
      resolution=global_config_parser.param('chromimpute', 'resolution'),
      path_to_dataset=input_dir,
      inputinfofile=inputinfofile,
      chrom_sizes=chr_sizes_file,
      output_dir=output_dir
      )
    )

def compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile, histone_mark, chr_sizes_file):

    return Job(
        (input_files),
       # input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_compute_global_dist." + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  ComputeGlobalDist \\
  -m {histone_mark} \\
  -r {resolution} \\
  {converteddir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir}""".format(
      java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
      ram=global_config_parser.param('chromimpute', 'ram'),
      histone_mark=histone_mark,
      resolution=global_config_parser.param('chromimpute', 'resolution'),
      converteddir=converteddir,
      inputinfofile=inputinfofile,
      chrom_sizes=chr_sizes_file,
      output_dir=output_dir
      )
    )

def generate_train_data(input_files, output_dir, output_files, converteddir, distancedir, inputinfofile, histone_mark, chr_sizes_file, chr):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_generate_train_data." + chr + "_" + histone_mark,
        command="""\

java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  GenerateTrainData \\
  -r {resolution} \\
  -c {chr} \\
  {converteddir} \\
  {distancedir} \\
  {inputinfofile} \\
  {chrom_sizes} \\
  {output_dir} \\
  {histone_mark}""".format(
      java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
      ram=global_config_parser.param('chromimpute_generate_train_data', 'ram'),
      resolution=global_config_parser.param('chromimpute', 'resolution'),
      converteddir=converteddir,
      distancedir=distancedir,
      inputinfofile=inputinfofile,
      chrom_sizes=chr_sizes_file,
      output_dir=output_dir,
      histone_mark=histone_mark,
      chr=chr
      )
    )

def temp_inputinfo(input_file, output_file):
    return Job(

        [input_file],
        [output_file],
        [],
        name="chromimpute_preprocess.temp_inputinfo",
        command="""\
sort -n {input_file} | \\
awk -v OFS="\\t" '{{if(NR==1){{sample=$1;rowindex=NR-1; print $0,rowindex}} else {{if(sample==$1){{print $0,rowindex }} else {{sample=$1;rowindex=rowindex+1;print $0,rowindex}} }} }}' > {output_file}""".format(
      input_file=input_file,
      output_file=output_file
      )
    )

def apply(input_dir, output, converteddir, distancedir, predictordir, inputinfofile, output_dir, sample, mark):
    return Job(
        ['chromimpute_metrics_dir', converteddir, distancedir, predictordir],
        [output],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name = "chromimpute_apply."+sample+"_"+mark,
        command = """\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Apply \\
    {chrom} \\
    {resolution} \\
    {converteddir} \\
    {distancedir} \\
    {predictordir} \\
    {inputinfofile} \\
    {chrom_sizes} \\
    {output_dir} \\
    {sample} \\
    {mark}""".format(
        java_other_options = global_config_parser.param('DEFAULT', 'java_other_options'),
        ram = global_config_parser.param('chromimpute', 'ram'),
        chrom = global_config_parser.param('chromimpute', 'chrom'),
        resolution = global_config_parser.param('chromimpute', 'resolution'),
        converteddir = converteddir,
        distancedir = distancedir,
        predictordir = predictordir,
        inputinfofile = inputinfofile,
        chrom_sizes = global_config_parser.param('chromimpute', 'chromsizes'),
        output_dir = output_dir,
        sample = sample,
        mark = mark
        )
    )


def train(input_files, output_dir, output_files, traindatadir, inputinfofile, sample, histone_mark):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_train." + sample + "_" + histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
  Train \\
  {traindatadir} \\
  {inputinfofile} \\
  {predictordir} \\
  {sample} \\
  {histone_mark}""".format(
      java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
      ram=global_config_parser.param('chromimpute', 'ram'),
      traindatadir=traindatadir,
      inputinfofile=inputinfofile,
      predictordir=output_dir,
      sample=sample,
      histone_mark=histone_mark
      )
    )

def apply(input_files, output_dir, converteddir, distancedir, predictordir, inputinfofile, sample, mark, chr, chr_sizes_file, output_files):
    return Job(
        input_files,
        output_files,
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_apply."+sample+"_"+mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Apply \\
    -r {resolution} \\
    -c {chrom} \\
    {converteddir} \\
    {distancedir} \\
    {predictordir} \\
    {inputinfofile} \\
    {chrom_sizes} \\
    {output_dir} \\
    {sample} \\
    {mark}""".format(
        java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
        ram=global_config_parser.param('chromimpute_apply', 'ram'),
        chrom=chr,
        resolution=global_config_parser.param('chromimpute', 'resolution'),
        converteddir=converteddir,
        distancedir=distancedir,
        predictordir=predictordir,
        inputinfofile=inputinfofile,
        chrom_sizes=chr_sizes_file,
        output_dir=output_dir,
        sample=sample,
        mark=mark
        )
    )


def eval(input_files, imputed_file, converted_file, output_file, converteddir, imputeddir, percent1, percent2, chr_sizes_file, sample, histone_mark):
    return Job(
        input_files,
        [output_file],
        [['java', 'module_java'], ['chromimpute', 'module_chromimpute']],
        name="chromimpute_eval."+sample+"_"+histone_mark,
        command="""\
java -Djava.io.tmpdir=$TMPDIR {java_other_options} -Xmx{ram} -jar $CHROMIMPUTE_JAR \\
    Eval \\
    -o {output_file} \\
    -p {percent1} {percent2} \\
    {converteddir} \\
    {converted_file} \\
    {imputeddir} \\
    {imputed_file} \\
    {chrom_sizes}""".format(
        java_other_options=global_config_parser.param('DEFAULT', 'java_other_options'),
        ram=global_config_parser.param('chromimpute_eval', 'ram'),
        percent1=percent1,
        percent2=percent2,
        converteddir=converteddir,
        converted_file=converted_file,
        imputeddir=imputeddir,
        imputed_file=imputed_file,
        chrom_sizes=chr_sizes_file,
        output_file=output_file
        )
    )
