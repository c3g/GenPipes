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

def create_fend_object(chromosome_lengths_file, output_file,output_dir, res_chr):

    return Job(
        [chromosome_lengths_file],
        [output_file],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            mkdir -p {output_dir} &&
            hifive fends \\
            -L {chromosome_lengths} \\
            --binned={res_chr} \\
            {output_file}""".format(
            chromosome_lengths=chromosome_lengths_file,
            output_file=output_file,
            res_chr=res_chr,
            output_dir=output_dir
        ))


def restructure_matrix(input_file_path, output_file_path,  output_dir):

    return Job(

        [input_file_path],
        [output_file_path],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            mkdir -p {output_dir}/temp &&
            awk -v OFS="\\t" '{{print $1":"$2"-"$3,$0}}' {input_file} | \\
            cut -f2,3,4 --complement | \\
            sed 's/ *//g' | \\
            Rscript -e 'library(data.table); a <- fread("file:///dev/stdin", data.table=F,sep="\\t", na.strings=c("",NA,"NULL")); a[is.na(a)] <- 0; rownames(a)<-a[,1]; a<- a[,-1]; if(!all(as.matrix(a)%%1==0)){{a<-a[,] *10;}};  a<-trunc(a[,]); write.table(a,"file:///dev/stdout",row.names=TRUE,col.names=FALSE,sep="\\t",quote=FALSE)'  > \\
            {output_file}""".format(
            input_file=input_file_path,
            output_file=output_file_path,
            output_dir=output_dir
        ))


def quality_analysis(quasr_temp_files, input_files, output_file_path, output_dir, fend_file, quasar_res, quasar_coverage):
    output_data_file = "_".join((
                    output_file_path, "hic.data"))

    output_project_file = "_".join((
                    output_file_path, "hic.project"))

    report_file="_".join((
                    output_file_path, "report.txt"))
    quasar_file="_".join((
                    output_file_path, "hic.quasar"))
    return Job(
        quasr_temp_files,
        [report_file],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            hifive hic-data -X "{input_files}" {fend_file} {output_data_file} &&
            hifive hic-project {output_data_file} {output_project_file} &&
            hifive quasar -p {output_project_file} -r {res} -d {coverage} -o {report_file} {quasar_file} --seed 12345
            """.format(
            input_files=input_files,
            output_data_file=output_data_file,
            output_project_file=output_project_file,
            output_dir=output_dir,
            fend_file=fend_file,
            report_file=report_file,
            quasar_file=quasar_file,
            res=quasar_res,
            coverage=quasar_coverage
        ))

