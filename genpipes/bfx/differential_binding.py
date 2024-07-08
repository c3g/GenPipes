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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

#This function is used to render R file and create a html output using knitr and spin
#This is a new feature introduced to Genpipes in 2021
def diffbind(
    input_files,
    comparison,
    design,
    readset,
    output_dir,
    alignment_dir,
    alignment_ext,
    peak_dir,
    minOverlap,
    minMembers,
    method
    ):

    output_file =  "".join((output_dir, "_".join(("/diffbind",comparison, method, "dba.txt"))))
    html_output = "".join((output_dir, "_".join(("/diffbind", comparison, method, "dba.html"))))
    R_filename = "".join((output_dir, "_".join(("/diffbind", comparison, method, "dba.R"))))

    th = global_conf.global_get('differential_binding', 'th')
    thT = global_conf.global_get('differential_binding', 'thT')
    bUsePval = global_conf.global_get('differential_binding', 'bUsePval')
    contrastnb = global_conf.global_get('differential_binding', 'contrastnb', required=False) if global_conf.global_get('differential_binding', 'contrastnb', required=False) else 1

    return Job(
        input_files,
        [output_file, html_output],
        [
            ['differential_binding', 'module_mugqic_tools'],
            ['differential_binding', 'module_R'],
            ['differential_binding', 'module_pandoc']
        ],
        command="""\
        mkdir -p {output_dir} &&
cp $R_TOOLS/DiffBind.R {R_filename} &&
Rscript -e 'cur_dir=getwd();library(knitr);rmarkdown::render("{R_filename}",params=list(cur_wd=cur_dir,d="{design}",r="{readset}",c="{comparison}",o="{output_file}",b="{alignment_dir}",ext="{alignment_ext}",p="{peak_dir}",dir="{output_dir}",minOverlap={minOverlap},minMembers={minMembers},method="{method}",th={th},thT={thT},bUsePval={bUsePval},contrastnb={contrastnb}),output_file=file.path(cur_dir,"{html_output}"));' &&
rm {R_filename}""".format(
            design=design,
            comparison=comparison,
            output_file=output_file,
            output_dir=output_dir,
            readset=readset,
            alignment_dir=alignment_dir,
            alignment_ext=alignment_ext,
            peak_dir=peak_dir,
            minOverlap=minOverlap,
            minMembers=minMembers,
            html_output=html_output,
            R_filename=R_filename,
            method=method,
            th=th,
            thT=thT,
            bUsePval=bUsePval,
            contrastnb=contrastnb
        )
    )

###The below function is not currently used. But if you want to use the old way to call Rscript passing paramters, use and modify this method.
#this function does not produce a HTML report but basic TSV result table.
def diffbind_R( input_files, comparison, design, readset, output_dir, alignment_dir, peak_dir, minOverlap, minMembers, method):

    output_file = "".join((output_dir, "_".join(("/diffbind", comparison, method, "dba.txt"))))

    th = global_conf.global_get('differential_binding', 'th')
    bUsePval = global_conf.global_get('differential_binding', 'bUsePval')

    return Job(
        input_files,
        [output_file],
        [
            ['differential_binding', 'module_mugqic_tools'],
            ['differential_binding', 'module_R']
        ],
        command="""\
        mkdir -p {output_dir} &&
Rscript $R_TOOLS/DiffBind.R \\
  -d {design} \\
  -r {readset} \\
  -c {comparison} \\
  -o {output_file} \\
  -b {alignment_dir} \\
  -p {peak_dir} \\
  -dir {output_dir} \\
  -minOverlap {minOverlap} \\
  -minMembers {minMembers} \\
  -method {method} \\
  -th {th} \\
  -bUsePval {bUsePval}""".format(
        design=design,
        comparison=comparison,
        output_file=output_file,
        output_dir=output_dir,
        readset=readset,
        alignment_dir=alignment_dir,
        peak_dir=peak_dir,
        minOverlap=minOverlap,
        minMembers=minMembers,
        method=method,
        th=th,
        bUsePval=bUsePval
    ))
