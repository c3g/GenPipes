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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def diffbind2(input_files, comparison, design, output_file):

    #merge all the tsvs and create final .csv file
    if not isinstance(input_files, list):
       input_files = [input_files]

    return Job(
        input_files,
        [output_file],
        [
            ['macs2_callpeak', 'module_macs2']
        ],
        command="""\
                    mkdir -p {output_dir}/""".format(
            output_dir=output_file
        ))

#This function is used to render R file and create a html output using knitr and spin
#This is a new feature introduced to Genpipes in 2021
def diffbind( input_files, comparison, design, readset, output_dir, alignment_dir, peak_dir, minOverlap, minMembers):

    output_file =  "".join((output_dir, "_".join(("/diffbind",comparison,"dba.txt"))))
    html_output = "".join((output_dir, "_".join(("/diffbind", comparison, "dba.html"))))
    R_filename = "".join((output_dir, "_".join(("/diffbind", comparison, "dba.R"))))

    return Job(
        input_files,
        [output_file],
        [
            ['differential_binding', 'module_mugqic_tools'],
            ['differential_binding', 'module_R']
        ],
        command="""\
        mkdir -p {output_dir} &&
        cp /home/pubudu/projects/rrg-bourqueg-ad/pubudu/chipseq_diff/analysis.R {R_filename} &&
#Rscript $R_TOOLS/diffbind.R \\
Rscript -e 'cur_dir=getwd();library(knitr);rmarkdown::render("{R_filename}",params=list(cur_wd=cur_dir,d="{design}",r="{readset}",c="{comparison}",o="{output_file}",b="{alignment_dir}",p="{peak_dir}",dir="{output_dir}",minOverlap="{minOverlap}",minMembers="{minMembers}"),output_file=file.path(cur_dir,"{html_output}"));' &&
rm {R_filename}""".format(
        design=design,
        comparison=comparison,
        output_file=output_file,
        output_dir=output_dir,
        readset=readset,
        alignment_dir=alignment_dir,
        peak_dir=peak_dir,
        minOverlap=minOverlap,
        minMembers=minMembers,
        html_output=html_output,
        R_filename=R_filename
    ))

# Rscript -e 'library(knitr);cur_dir=getwd();design=paste0(cur_dir,"/","{design}");
# readset=paste0(cur_dir,"/","{readset}");output_file=paste0(cur_dir,"/","{output_file}");
# alignment_dir=paste0(cur_dir,"/","{alignment_dir}");peak_dir=paste0(cur_dir,"/","{peak_dir}");
# output_dir=paste0(cur_dir,"/","{output_dir}");html_output=paste0(cur_dir,"/","{html_output}");
# rmarkdown::render("/home/pubudu/projects/rrg-bourqueg-ad/pubudu/chipseq_diff/analysis.R",params=list(d=design,r=readset,c="{comparison}",o=output_file,b=alignment_dir,p=peak_dir,dir=output_dir,minOverlap="{minOverlap}",minMembers="{minMembers}"),output_file=html_output);'""".format(
###The function is not currently used. But if you want to use the old way to call Rscript passing paramters use and modify this method.
def diffbind_R( input_files, comparison, design, readset, output_dir, alignment_dir, peak_dir, minOverlap, minMembers):

    output_file =  "".join((output_dir, "_".join(("/diffbind",comparison,"dba.txt"))))

    return Job(
        input_files,
        [output_file],
        [
            ['differential_binding', 'module_mugqic_tools'],
            ['differential_binding', 'module_R']
        ],
        command="""\
        mkdir -p {output_dir} &&
#Rscript $R_TOOLS/diffbind.R \\
Rscript /home/pubudu/projects/rrg-bourqueg-ad/pubudu/chipseq_diff/analysis.R \\
  -d {design} \\
  -r {readset} \\
  -c {comparison} \\
  -o {output_file} \\
  -b {alignment_dir} \\
  -p {peak_dir} \\
  -dir {output_dir} \\
  -minOverlap {minOverlap} \\
  -minMembers {minMembers}""".format(
        design=design,
        comparison=comparison,
        output_file=output_file,
        output_dir=output_dir,
        readset=readset,
        alignment_dir=alignment_dir,
        peak_dir=peak_dir,
        minOverlap=minOverlap,
        minMembers=minMembers
    ))
