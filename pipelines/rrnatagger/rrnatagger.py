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
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

from bfx import rrna_amplicons
#from bfx import microbial_ecology

from pipelines import common

# Global scope variables
log = logging.getLogger(__name__)

class RRNATagger(common.Illumina):
    #outdir             = config.param('default', 'currentDir')
    #clustering_method  = config.param('clustering', 'clusteringMethod')
    #lib_type           = config.param('default', 'libraryType')
    #project_name       = config.param('default', 'projectName')
    #if not project_name:
    #    project_name       = "RRNATagger-project" 
    #organism_type      = config.param('default', 'organism_type')
    #barcodes           = os.path.abspath(self.args.barcodes)
    
    #outdir 
    #clustering_method  
    #lib_type
    #project_name
    #organism_type
    #barcodes

    def merge_barcodes(self):
        # Merge all demultiplexed fastq files in one file. One file for reads1 and one file for reads2 of Illumina paired-end.
        # If library is Illumina single end or 454 or PacBio, only one file will be generated.
        #lib_type = config.param('default', 'libraryType', type='string')
        #raw_reads_dir = os.path.join("raw_reads", sample.name)
        outdir = os.path.curdir + "/raw_reads"
        jobs = []
        reads1 = []
        reads2 = []
        
        if self.lib_type == "ex":
            for readset in self.readsets: 
                reads1.append(os.path.join("raw_reads", readset.sample.name, "run" + readset.run + "_" + readset.lane, readset.sample.name + "." + readset.library + "." + str(readset.quality_offset) + ".pair1.fastq.gz"))
                reads2.append(os.path.join("raw_reads", readset.sample.name, "run" + readset.run + "_" + readset.lane, readset.sample.name + "." + readset.library + "." + str(readset.quality_offset) + ".pair2.fastq.gz"))
                # exemple of path raw_reads/11_a_25/runM00833_0215_1/11_a_25.1000023354_HSP-A11.33.pair1.fastq.gz
            
            job = rrna_amplicons.merge_barcodes(
                reads1,
                reads2,
                outdir
            )
            job.name = "merge_barcodes"
            jobs.append(job)

        elif self.lib_type == "nc1":
            for readset in self.readsets: 
                reads1.append(os.path.join("raw_reads", readset.sample.name, "run" + readset.run + "_" + readset.lane, readset.sample.name + "." + readset.library + "." + str(readset.quality_offset) + ".pair1.fastq.gz"))
            
            job = rrna_amplicons.merge_barcodes_single_end_reads(
                reads1,
                outdir,
                "1"
            )
            job.name = "merge_barcodes_reads1"
            jobs.append(job)
        
        elif self.lib_type == "nc2":
            for readset in self.readsets: 
                reads2.append(os.path.join("raw_reads", readset.sample.name, "run" + readset.run + "_" + readset.lane, readset.sample.name + "." + readset.library + "." + str(readset.quality_offset) + ".pair2.fastq.gz"))
            
            job = rrna_amplicons.merge_barcodes_single_end_reads(
                reads2,
                outdir,
                "2"
            )
            job.name = "merge_barcodes_reads2"
            jobs.append(job)
        
        else:
            raise Exception("Error: run type \"" + readset.run_type +
            "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
                
        return jobs


    def remove_contam(self):
        duk_dir = outdir + "/duk"

        job = rrna_amplicons.duk_wrapper(
            infile_fastq,
            duk_dir + "/contam.fastq",
            duk_dir + "/ncontam.fastq",
            duk_dir + "/duk_contam_log.txt",
            config.param('DB', 'contaminants')
        );
        job.name = "duk_wrapper"
        jobs.append(job)

        sys.stderr.write('not implemented yet\n')

    def split_barcodes(self):
        sys.stderr.write('not implemented yet\n')

    def qscores(self):
        sys.stderr.write('not implemented yet\n')

    def qc(self):
        sys.stderr.write('not implemented yet\n')

    def qscores_qced(self):
        sys.stderr.write('not implemented yet\n')

    def generate_clusters(self):
        sys.stderr.write('not implemented yet\n')

    def classify(self): # Add PacBio's blast here!
        sys.stderr.write('not implemented yet\n')
    
    def generate_otu_table(self):
        sys.stderr.write('not implemented yet\n')

    def filter_otu_table(self):
        sys.stderr.write('not implemented yet\n')

    def summarize_taxonomy(self):
        sys.stderr.write('not implemented yet\n')

    def plot_taxonomy(self):
        sys.stderr.write('not implemented yet\n')

    def deliverables(self):
        sys.stderr.write('not implemented yet\n')

    def cleanup(self):
        sys.stderr.write('not implemented yet\n')
    
    # Override illumina.py readsets to make sure we are parsing a nanuq sample sheet
    # and not a readset sheet.
    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = parse_nanuq_readset_file(self.args.readsets.name)
        return self._readsets

    @property
    def steps(self):
        return [
            self.merge_barcodes
        ]

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-b", "--barcodes", help="barcodes in fasta format", type=file)
        sys.stderr.write('Running rRNA amplicons pipeline\n')
        
        
        #self.outdir = config.param('default', 'currentDir', type='filePath')
        self.clustering_method = config.param('clustering', 'clusteringMethod', type='int')
        #lib_type = config.param('default', 'libraryType', type='string')
        self.project_name = config.param('default', 'projectName', type='string')
        if not self.project_name:
            self.project_name = "RRNATagger-project" 
        #organism_type = config.param('default', 'organism_type', type='string')
        self.lib_type = config.param('default', 'libraryType', type='string')
        self.barcodes = os.path.abspath(self.args.barcodes)
        
        sys.stderr.write('lib_type:' + self.lib_type + '\n')
        super(RRNATagger, self).__init__()
       
if __name__ == "__main__":
    # RRNATagger pipeline is under development
    raise NotImplementedError
    RRNATagger()

# Steps from old Perl pipeline:

# MERGE
#1 mergeBarcodes

# REMOVE CONTAM
#2 duk_wrapper_contam
#3 duk_wrapper_phix

# SPLIT BY BARCODES
#4 barcodes

# QSCORES
#5 qscore_sheet_R1_raw
#6 qscore_plot_R1_raw

# QC
#7 itags_QC

# QSCORE QCED
#8 qscore_sheet_R1QCed
#9 qscore_plot_R1QCed

# CLUSTER SEQUENCED
#10 clustering_method_3

# CLASSIFY
#11 RDP 

# GENERATE RAW OTU TABLE
#12 add_taxonomy

# FILTER OTU TABLE
#13 filter_obs
#14 split_otu_table
#15 convert_otu_to_biom
#16 convert_otu_to_biom2
#17 convert_otu_to_biom3
#18 rarefy
#19 rarefy
#20 filter_obs_table2X2
#21 filter_obs_table2X2
#22 filter_obs_table1X1
#23 filter_obs_table1X1

# SUMMARIZE TAXONOMY
#24 summarize_taxonomy_absolute_L1
#25 summarize_taxonomy_absolute_L2
#26 summarize_taxonomy_absolute_L3
#27 summarize_taxonomy_absolute_L4
#28 summarize_taxonomy_absolute_L5
#29 summarize_taxonomy_absolute_L6
#30 summarize_taxonomy_absolute_L7
#31 summarize_taxonomy_relative_L1
#32 summarize_taxonomy_relative_L2
#33 summarize_taxonomy_relative_L3
#34 summarize_taxonomy_relative_L4
#35 summarize_taxonomy_relative_L5
#36 summarize_taxonomy_relative_L6
#37 summarize_taxonomy_relative_L7
#38 summarize_taxonomy_absolute_raw_L1
#39 summarize_taxonomy_absolute_raw_L2
#40 summarize_taxonomy_absolute_raw_L3
#41 summarize_taxonomy_absolute_raw_L4
#42 summarize_taxonomy_absolute_raw_L5
#43 summarize_taxonomy_absolute_raw_L6
#44 summarize_taxonomy_absolute_raw_L7
#45 summarize_taxonomy_relative_raw_L1
#46 summarize_taxonomy_relative_raw_L2
#47 summarize_taxonomy_relative_raw_L3
#48 summarize_taxonomy_relative_raw_L4
#49 summarize_taxonomy_relative_raw_L5
#50 summarize_taxonomy_relative_raw_L6
#51 summarize_taxonomy_relative_raw_L7

# PLOT TAXONOMY
#52 plot_taxonomy_absolute
#53 plot_taxonomy_absolute_raw
#54 phylum_barplot_all
#55 phylum_barplot_all_rel
#56 phylum_barplot_bacteria_rarefied

# BLAST (PACBIO ONLY)
#57 blast raw OTUs

# DELIVERABLES
#58 count_report
#59 txtToPdf
#60 merge_pdf
#61 Nozzle deliverables

# CLEANUP
#62 cleanup
