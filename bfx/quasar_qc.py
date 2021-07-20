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

def create_fend_object(chromosome_lengths_file, output_file,output_dir, res_chr, temp_dir):

    return Job(
        [chromosome_lengths_file],
        [output_file],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            mkdir -p {output_dir}/{temp_dir} &&
            hifive fends \\
            -L {chromosome_lengths} \\
            --binned={res_chr} \\
            {output_file}""".format(
            chromosome_lengths=chromosome_lengths_file,
            output_file=output_file,
            res_chr=res_chr,
            output_dir=output_dir,
            temp_dir=temp_dir
        ))


def restructure_matrix(input_file_path, output_file_path,  output_dir, bin, temp_dir):

    return Job(

     [input_file_path],
        [output_file_path],
        [
            ['quality_scores', 'module_R']
        ],
        command="""\
            mkdir -p {output_dir}/{temp_dir} &&
            awk -v OFS="\\t" '{{if(NR!=1) {{split($1,a,"-"); print a[1]":"a[2]"-"{bin}+a[2],$0 }} else {{print "V1",$0}} }}' {input_file} | \\
            cut -f2 --complement | \\
            Rscript -e 'library(data.table); a <- fread("file:///dev/stdin", data.table=F,sep="\\t", na.strings=c("",NA,"NULL")); a[is.na(a)] <- 0; rownames(a)<-a[,1]; a<- a[,-1]; if(!all(as.matrix(a)%%1==0)){{a<-a[,] *10;}};  write.table(a,"file:///dev/stdout",row.names=TRUE,col.names=FALSE,sep="\\t",quote=FALSE)' | \\
            awk '{{ printf $1"\\t"; for( i = 2; i <= NF; i++) {{ printf "%.0f\\t", $i;  }} printf "\\n"}}' > \\
            {output_file}""".format(
            input_file=input_file_path,
            output_file=output_file_path,
            output_dir=output_dir,
            bin=bin,
            temp_dir=temp_dir
        ))


def quality_analysis(quasr_temp_files, input_files, output_file_prefix, output_dir, fend_file, quasar_res, quasar_coverage, report_file):
    output_data_file = "_".join((
                    output_file_prefix, "hic.data"))

    output_project_file = "_".join((
                    output_file_prefix, "hic.project"))

    quasar_file="_".join((
                    output_file_prefix, "hic.quasar"))
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


def merge_all_reports(out_dir, input_files, output_file, res, enzyme, quasr_prefix,temp_dir, q_res ):

    return Job(
        input_files,
        [output_file],
        [
            ['quality_scores', 'module_python']
        ],
        command="""\
            if test -f "{output_file}"; then
            rm {output_file} 
            fi &&
            i=1 &&
            for f in {out_dir}/{temp_dir}/*_{quasr_prefix}_{res}_{enzyme}_report.txt ;do awk -v OFS="\\t" '{{if($1=="Resolution") {{print $0}} else {{print FILENAME,$0 }} }}' $f | \\
            tr -d " \\r" | \\
            awk -v OFS="\\t" '{{if(NF>3){{print $0}} }}' | \\
            awk -v i="$i" -v OFS="\\t" '{{if(NR==1 && i==1){{print "Sample","Sample",$0}} else {{if($2=="{quasar_resolution}") {{split($1,fname,"_{quasr_prefix}_"); split(fname[1],sname,"{out_dir}/{temp_dir}/"); print sname[2], $0}} }} }}' | \\
            cut -f2 --complement >> \\
            {output_file}; i=$(($i+1)); done
            """.format(
            input_files=input_files,
            output_file=output_file,
            res=res,
            enzyme=enzyme,
            out_dir=out_dir,
            quasr_prefix=quasr_prefix,
            temp_dir=temp_dir,
            quasar_resolution=q_res

        ))

