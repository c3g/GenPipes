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
import re

# MUGQIC Modules
from core.config import *
from core.job import *

##### General Homer:

def archive_contigs(homerTagDir, output_dir = "archive"):

    command_archive = """cd {homerTagDir} && \\
mkdir -p {output_dir} && \\
mv -t {output_dir} *random*.tsv *chrUn*.tsv *hap*.tsv chrM*.tsv MT*.tsv *Y*.tsv *EBV*.tsv *GL*.tsv NT_*.tsv || echo "not all files found"
cd ../../""".format(
        homerTagDir=homerTagDir,
        output_dir=output_dir)

    return Job(input_files = [os.path.join(homerTagDir, "tagInfo.txt")],
            output_files = [os.path.join(homerTagDir, output_dir)],
            command = command_archive
        )

def makeTagDir(output_dir, input_bam, genome, restriction_site=None, illuminaPE=False, other_options=None):
    ## if hic experiment then input_bam should be provided twice: {input_bam},{input_bam}
    command_tagDir = """makeTagDirectory {output_dir} \\
    {input_bam} \\
    -genome {genome} \\
    -checkGC{illuminaPE}{restriction_site}{other_options}""".format(
                output_dir=output_dir,
                input_bam=input_bam,
                genome=genome,
                illuminaPE=" \\\n -illuminaPE" if illuminaPE else "",
                restriction_site=" \\\n -restrictionSite " + restriction_site if restriction_site is not None else "",
                other_options=" \\\n " + other_options if other_options is not None else ""
                )

    input_bam_file = input_bam.split(",")[0]

    return Job(input_files=[input_bam_file],
            output_files=[os.path.join(output_dir, "tagInfo.txt"), os.path.join(output_dir, "tagGCcontent.txt"), os.path.join(output_dir, "genomeGCcontent.txt"), os.path.join(output_dir, "tagLengthDistribution.txt")],
            module_entries=[["homer_tag_directory", "module_perl"], ["homer_tag_directory", "module_homer"], ["homer_tag_directory", "module_samtools"]],
            command=command_tagDir,
            name="homer_make_tag_directory"
            )



##### Homer for HiC:



def hic_tagDirQcPlots (name, working_dir, output_dir = "HomerQcPlots"):

    command_QcPlots="Rscript {script} {name} {working_dir} {output_dir}".format(
        script=os.path.expandvars("${R_TOOLS}/HomerHiCQcPlotGenerator.R"),  
        name=name, 
        working_dir=working_dir, 
        output_dir=output_dir)

    return Job(input_files=[os.path.join(working_dir, "tagInfo.txt")],
            output_files=[os.path.join(working_dir, output_dir)],
            module_entries=[["homer_tag_directory", "module_R"], ["homer_tag_directory", "module_mugqic_tools"]],
            command=command_QcPlots
            )



def hic_interactionMatrix_chr (name, output_dir, homer_dir, res, chr, fileName, fileNameRN=None, norm="raw", format = True):

# norm can be "raw", "simpleNorm", "norm"

    commandChrMatrix = """mkdir -p {output_dir} && \\
                analyzeHiC {homer_dir} -res {res} -{norm} -chr {chr} > {fileName}""".format(
                output_dir=output_dir,
                homer_dir=homer_dir,
                res=res,
                norm=norm,
                chr=chr,
                fileName=fileName,
                fileNameRN=fileNameRN)

    if fileNameRN is None:
        fileNameRN = re.sub("\.txt", "RN.txt", fileName)

    commandFormatMatrix = """cut -f 2- {fileName} > {fileNameRN}""".format(
        fileName=fileName,
        fileNameRN=fileNameRN)

    if format:
        command = commandChrMatrix + " && " + commandFormatMatrix
    else:
        command = commandChrMatrix

    return Job(input_files=[os.path.join(homer_dir, "tagInfo.txt")],
        output_files=[fileNameRN, fileName],
        module_entries=[["interaction_matrices_Chr", "module_perl"], ["interaction_matrices_Chr", "module_homer"]],
        command=command,
        name="interaction_matrices_Chr." + name + "_" + chr + "_res" + res,
        removable_files=[fileName]
        )

def hic_interactionMatrix_genome(name, output_dir, homer_dir, res, fileName, fileNameRN=None, norm="raw", format = True):

    commandMatrix="""mkdir -p {output_dir} && \\
analyzeHiC {homer_dir} -res {res} -{norm} > {fileName}""".format(
        output_dir=output_dir,
        homer_dir=homer_dir,
        res=res,
        norm=norm,
        fileName=fileName
        )

    commandFormatMatrix = """cut -f 2- {fileName} > {fileNameRN}""".format(
        fileName=fileName,
        fileNameRN=fileNameRN)

    if fileNameRN is None:
        fileNameRN = re.sub("\.txt", "RN.txt", fileName)


    if format:
        command = commandMatrix + " && " + commandFormatMatrix
    else:
        command = commandMatrix


    return Job(input_files=[os.path.join(homer_dir, "tagInfo.txt")],
                    output_files=[fileNameRN, fileName],
                    module_entries=[["interaction_matrices_genome", "module_perl"], ["interaction_matrices_genome", "module_homer"]],
                    name="interaction_matrices_genome." + name  + "_res" + res,
                    command=command,
                    removable_files=[fileName]
                    )

def hic_compartments (name, output_dir, fileName, homer_dir, res, genome, fileName_PC1, fileName_Comp, cpu):

    command = """mkdir -p {output_dir} && \\
runHiCpca.pl {fileName} {homer_dir} -res {res} -genome {genome} -cpu {cpu} && \\
findHiCCompartments.pl {fileName_PC1}  > {fileName_Comp}""".format(
        output_dir=output_dir,
        fileName=fileName,
        homer_dir=homer_dir,
        res=res,
        genome=genome,
        cpu=cpu,
        fileName_PC1=fileName_PC1,
        fileName_Comp=fileName_Comp)

    return Job(input_files=[os.path.join(homer_dir, "tagInfo.txt")],
            output_files=[fileName_PC1, fileName_Comp],
            module_entries=[["identify_compartments", "module_perl"], ["identify_compartments", "module_homer"], ["identify_compartments", "module_R"]],
            name="identify_compartments." + name + "_res" + res,
            command=command,
            removable_files=[fileName]
            )


def hic_peaks(name, output_dir, homer_dir, res, genome, fileName, fileName_anno, cpu=1):

    command = """mkdir -p {output_dir} && \\
findHiCInteractionsByChr.pl {homer_dir} -res {res} -cpu {cpu} > {fileName} && \\
annotateInteractions.pl {fileName} {genome} {fileName_anno}""".format(
        output_dir=output_dir,
        homer_dir=homer_dir,
        res=res,
        fileName=fileName,
        genome=genome,
        fileName_anno=fileName_anno,
        cpu=cpu)


    return Job(input_files=[os.path.join(homer_dir, "tagInfo.txt")],
        output_files=[fileName, fileName_anno],
        module_entries=[["identify_peaks", "module_perl"], ["identify_peaks", "module_homer"]],
        name="identify_peaks." + name + "_res" + res,
        command=command
        )

##### Homer for ChIP-Seq:


def makeUCSCfile(tag_dir, bedgraph_file):

    bedgraph_file_gz = bedgraph_file + ".gz"

    cmd = """makeUCSCfile \\
{tag_dir} > {bedgraph_file} && \\
gzip -c {bedgraph_file} > {bedgraph_file_gz}""".format(
                    tag_dir=tag_dir,
                    bedgraph_file=bedgraph_file,
                    bedgraph_file_gz=bedgraph_file_gz)

    return Job(input_files=[os.path.join(tag_dir, "tagInfo.txt")],
            output_files=[bedgraph_file, bedgraph_file_gz],
            module_entries=[["homer_make_ucsc_files", "module_perl"], ["homer_make_ucsc_files", "module_homer"]],
            command=cmd,
            name="homer_make_ucsc_file"
            )


def annotatePeaks(peak_file, genome, output_dir, annotation_file, genome_size):
    cmd = """annotatePeaks.pl \\
    {peak_file} \\
    {genome} \\
    -gsize {genome_size} \\
    -cons -CpG \\
    -go {output_dir} \\
    -genomeOntology {output_dir} > {annotation_file}""".format(
                            peak_file=peak_file,
                            genome=genome,
                            genome_size=genome_size,
                            output_dir=output_dir,
                            annotation_file=annotation_file)

    return Job(input_files=[peak_file],
            output_files=[
                annotation_file,
                os.path.join(output_dir, "geneOntology.html"),
                os.path.join(output_dir, "GenomeOntology.html")
                ],
            module_entries=[
                ["homer_annotate_peaks", "module_perl"],
                ["homer_annotate_peaks", "module_homer"]
                ],
            command=cmd,
            name="homer_make_ucsc_file"
            )

def findMotifsGenome(peak_file, genome, output_dir, threads):

    cmd = """findMotifsGenome.pl \\
  {peak_file} \\
  {genome} \\
  {output_dir} \\
  -preparsedDir {output_dir}/preparsed \\
  -p {threads}""".format(
                        peak_file=peak_file,
                        genome=genome,
                        output_dir=output_dir,
                        threads=threads,
                        name="homer_find_motifs_genome"
                )

    return Job(input_files=[peak_file],
            output_files=[
                os.path.join(output_dir, "homerResults.html"),
                os.path.join(output_dir, "knownResults.html")
                ],
            module_entries=[
                ["homer_find_motifs_genome", "module_perl"],
                ["homer_find_motifs_genome", "module_homer"],
                ["homer_find_motifs_genome", "module_weblogo"]
                ],
            command=cmd,
            name="homer_find_motifs_genome"
            )
