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
import re

# MUGQIC Modules
from core.config import *
from core.job import *

def check_readName_format (input_bam, illuminaPE):
	
	## if bam read names end with /1 and /2 then -illuminaPE flag should be used, otherwise not. Check first read in bam to see is it contains /1  or /2

	## extract first read name in bam
	command = """readID=$(samtools view {input_bam} | head -n 1 | cut -f 1)
	
    if grep -q '^.*/[12]$' <<< $readID; then
        if [ "$illuminaPE" = "False" ]; then
        echo "Error! Your bam readIds contain /1 , /2 endings. Please edit the ini File to illuminaPE=True OR run the pipeline with the fastq_readName_Edit method!" 1>&2
        exit 1
        fi
    else
       if [ "$illuminaPE" = "True" ]; then
        echo "Error! Your bam readIds do not contain /1 , /2 endings. Please edit the ini File to illuminaPE=False OR run the pipeline without the fastq_readName_Edit method!" 1>&2
        exit 1
        fi
    fi""".format(input_bam = input_bam, illuminaPE = illuminaPE)

	return Job(input_files = [input_bam],
            module_entries = [["homer_tag_directory", "module_samtools"]],
            command = command
            )


def makeTagDir_hic (output_dir, input_bam, genome, restriction_site, illuminaPE=False, ini_section='homer_tag_directory'):

	command_tagDir = """rm -rf {output_dir} && \\
		 	makeTagDirectory {output_dir} {input_bam},{input_bam} \\
		 	-genome {genome} \\
		 	-restrictionSite {restriction_site} \\
		 	-checkGC{illuminaPE}{makeDirTag_hic_other_options}""".format(
		 		output_dir = output_dir, 
		 		input_bam = input_bam, 
		 		genome = genome, 
		 		restriction_site = restriction_site,
		 		illuminaPE=" \\\n -illuminaPE" if illuminaPE else "",
				makeDirTag_hic_other_options=config.param(ini_section, 'other_options', required=False))


	return Job(input_files = [input_bam],
            output_files = [output_dir],
            module_entries = [["homer_tag_directory", "module_perl"], ["homer_tag_directory", "module_homer"], ["homer_tag_directory", "module_samtools"]],
            command = command_tagDir
            )


def tagDirQcPlots_hic (name, working_dir, output_dir = "HomerQcPlots"):

	command_QcPlots = "Rscript {script} {name} {working_dir} {output_dir}".format(
		script = os.path.expandvars("${R_TOOLS}/HomerHiCQcPlotGenerator.R"),  
		name = name, 
		working_dir = working_dir, 
		output_dir = output_dir)

	return Job(input_files = [working_dir],
            output_files = [os.path.join(working_dir, output_dir)],
            module_entries = [["homer_tag_directory", "module_R"], ["homer_tag_directory", "module_mugqic_tools"]],
            command = command_QcPlots
            )


def archive_contigs_hic (homerTagDir, output_dir = "archive"):

    command_archive = """cd {homerTagDir} && \\
    mkdir -p {output_dir} && \\
    mv -t {output_dir} *random*.tsv *chrUn*.tsv *hap*.tsv chrM*.tsv MT*.tsv *Y*.tsv *EBV*.tsv *GL*.tsv NT_*.tsv || echo "not all files found"
    cd ../../""".format(
    	homerTagDir = homerTagDir, 
    	output_dir = output_dir)

    return Job(input_files = [homerTagDir],
            output_files = [os.path.join(homerTagDir, output_dir)],
            command = command_archive
        )


def interactionMatrix_chr_hic (name, output_dir, homer_dir, res, chr, fileName, fileNameRN=None, norm="raw", format = True):

# norm can be "raw", "simpleNorm", "norm"

    commandChrMatrix = """mkdir -p {output_dir} && \\
                analyzeHiC {homer_dir} -res {res} -{norm} -chr {chr} > {fileName}""".format(
                output_dir = output_dir, 
                homer_dir = homer_dir, 
                res = res, 
                norm = norm,
                chr = chr, 
                fileName = fileName, 
                fileNameRN = fileNameRN)


    if fileNameRN is None:
        fileNameRN = re.sub("\.txt", "RN.txt", fileName)


    commandFormatMatrix = """cut -f 2- {fileName} > {fileNameRN}""".format(
        fileName = fileName, 
        fileNameRN = fileNameRN)

    if format:
        command = commandChrMatrix + " && " + commandFormatMatrix
    else:
        command = commandChrMatrix
                    
    return Job(input_files = [homer_dir],
        output_files = [fileNameRN, fileName],
        module_entries = [["interaction_matrices_Chr", "module_perl"], ["interaction_matrices_Chr", "module_homer"]],
        command = command,
        name = "interaction_matrices_Chr.plotting." + name + "_" + chr + "_res" + res,
        removable_files = [fileName]
        )

def interactionMatrix_genome_hic(name, output_dir, homer_dir, res, fileName, fileNameRN=None, norm="raw", format = True):

    commandMatrix="""mkdir -p {output_dir} && \\
            analyzeHiC {homer_dir} -res {res} -{norm} > {fileName}""".format(
        output_dir = output_dir, 
        homer_dir = homer_dir, 
        res = res,
        norm = norm, 
        fileName = fileName
        )

    commandFormatMatrix = """cut -f 2- {fileName} > {fileNameRN}""".format(
        fileName = fileName, 
        fileNameRN = fileNameRN)

    if fileNameRN is None:
        fileNameRN = re.sub("\.txt", "RN.txt", fileName)


    if format:
        command = commandMatrix + " && " + commandFormatMatrix
    else:
        command = commandMatrix


    return Job(input_files = [homer_dir],
                    output_files = [fileNameRN, fileName],
                    module_entries = [["interaction_matrices_genome", "module_perl"], ["interaction_matrices_genome", "module_homer"]],
                    name = "interaction_matrices_genome." + name  + "_res" + res,
                    command = command,
                    removable_files = [fileName]
                    )

def compartments_hic (name, output_dir, fileName, homer_dir, res, genome, fileName_PC1, fileName_Comp, cpu):

    command = """mkdir -p {output_dir} && \\
    runHiCpca.pl {fileName} {homer_dir} -res {res} -genome {genome} -cpu {cpu} && \\
    findHiCCompartments.pl {fileName_PC1}  > {fileName_Comp}""".format(
        output_dir = output_dir, 
        fileName = fileName, 
        homer_dir = homer_dir, 
        res = res, 
        genome = genome,
        cpu = cpu, 
        fileName_PC1 = fileName_PC1, 
        fileName_Comp = fileName_Comp)

    return Job(input_files = [homer_dir],
            output_files = [fileName_PC1, fileName_Comp],
            module_entries = [["identify_compartments", "module_perl"], ["identify_compartments", "module_homer"], ["identify_compartments", "module_R"]],
            name = "identify_compartments." + name + "_res" + res,
            command = command,
            removable_files = [fileName]
            )


def peaks_hic(name, output_dir, homer_dir, res, genome, fileName, fileName_anno, cpu=1):

    command = """mkdir -p {output_dir} && \\
    findHiCInteractionsByChr.pl {homer_dir} -res {res} -cpu {cpu} > {fileName} && \\
    annotateInteractions.pl {fileName} {genome} {fileName_anno}""".format(
        output_dir = output_dir, 
        homer_dir = homer_dir, 
        res = res, 
        fileName = fileName, 
        genome = genome, 
        fileName_anno = fileName_anno, 
        cpu = cpu)


    return Job(input_files = [homer_dir],
        output_files = [fileName, fileName_anno],
        module_entries = [["identify_peaks", "module_perl"], ["identify_peaks", "module_homer"]],
        name = "identify_peaks." + name + "_res" + res,
        command = command
        )



