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

def exploratory_analysis_rnaseq(htseq_count_file, cuffnorm_dir, genes_file, output_dir):

    return Job(
        [htseq_count_file, os.path.join(cuffnorm_dir, "isoforms.fpkm_table"), os.path.join(cuffnorm_dir, "isoforms.attr_table")],
        [os.path.join(output_dir, "index.tsv"), os.path.join(output_dir, "top_sd_heatmap_cufflinks_logFPKMs.pdf")],
        [
            ['gq_seq_utils_exploratory_analysis_rnaseq', 'module_R'],
            ['gq_seq_utils_exploratory_analysis_rnaseq', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryAnalysisRNAseq(htseq.counts.path="{htseq_count_file}", cuffnorm.fpkms.dir="{cuffnorm_dir}", genes.path="{genes_file}", output.dir="{output_dir}")
desc = readRDS(file.path("{output_dir}","index.RData"))
write.table(desc,file=file.path("{output_dir}","index.tsv"),sep='\t',quote=F,col.names=T,row.names=F)
print("done.")

EOF""".format(
        htseq_count_file=htseq_count_file,
        cuffnorm_dir=cuffnorm_dir,
        genes_file=genes_file,
        output_dir=output_dir
    ))

def exploratory_analysis_rnaseq_denovo(count_file, genes_file, output_dir):

    return Job(
        [count_file, genes_file ],
        [os.path.join(output_dir, "index.tsv"), os.path.join(output_dir, "top_sd_heatmap_log2CPM.pdf")],
        [
            ['gq_seq_utils_exploratory_analysis_rnaseq_denovo', 'module_R'],
            ['gq_seq_utils_exploratory_analysis_rnaseq_denovo', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))
exploratoryAnalysisRNAseqdenovo(read.counts.path="{count_file}", genes.path="{genes_file}", output.dir="{output_dir}")
desc = readRDS(file.path("{output_dir}","index.RData"))
write.table(desc,file=file.path("{output_dir}","index.tsv"),sep='\t',quote=F,col.names=T,row.names=F)
print("done.")

EOF""".format(
        count_file=count_file,
        genes_file=genes_file,
        output_dir=output_dir
    ))

def exploratory_analysis_rnaseq_light(count_file, genes_file, output_dir):

    return Job(
        [count_file, genes_file ],
        [os.path.join(output_dir, "index.tsv"), os.path.join(output_dir, "top_sd_heatmap_log2CPM.pdf")],
        [
            ['gq_seq_utils_exploratory_analysis_rnaseq_light', 'module_R'],
            ['gq_seq_utils_exploratory_analysis_rnaseq_light', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))
exploratoryAnalysisRNAseqLightKallisto(read.counts.path="{count_file}", genes.path="{genes_file}", output.dir="{output_dir}")
desc = readRDS(file.path("{output_dir}","index.RData"))
write.table(desc,file=file.path("{output_dir}","index.tsv"),sep='\t',quote=F,col.names=T,row.names=F)
print("done.")

EOF""".format(
        count_file=count_file,
        genes_file=genes_file,
        output_dir=output_dir
    ))

def report(ini_filepaths, project_path, pipeline_type, output_directory):

    title = global_config_parser.param('gq_seq_utils_report', 'report_title', required=False)
    path = os.path.join(output_directory, global_config_parser.param('gq_seq_utils_report', 'report_dir'))
    author = global_config_parser.param('gq_seq_utils_report', 'report_author', required=False)
    contact = global_config_parser.param('gq_seq_utils_report', 'report_contact', required=False)

    return Job(
        # Job input files must be set in pipeline class since they are different for each pipeline
        [],
        [os.path.join(path, "index.html")],
        [
            ['gq_seq_utils_report', 'module_R'],
            ['gq_seq_utils_report', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save -e 'library(gqSeqUtils); mugqicPipelineReport(pipeline="{pipeline}"{title}{path}{author}{contact}, ini.file.path=c({ini_filepaths}), project.path="{project_path}")'""".format(
        pipeline=pipeline_type,
        title=", report.title=\"" + title + "\"" if title else "",
        path=", report.path=\"" + path + "\"" if path else "",
        author=", report.author=\"" + author + "\"" if author else "",
        contact=", report.contact=\"" + contact + "\"" if contact else "",
        ini_filepaths=",".join(['"' + os.path.abspath(ini_filepath) + '"' for ini_filepath in ini_filepaths]),
        project_path=project_path
    ))
