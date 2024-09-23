################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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
import re

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job 

##### General Homer:
def archive_contigs(
    homerTagDir,
    output_dir="archive"
    ):
    return Job(
        input_files=[
            os.path.join(homerTagDir, "tagInfo.txt")
        ],
        output_files=[
            os.path.join(homerTagDir, output_dir)
        ],
        command="""\
cd {homerTagDir} && \\
mkdir -p {output_dir} && \\
mv -t {output_dir} *random*.tsv *chrUn*.tsv *hap*.tsv chrM*.tsv MT*.tsv *Y*.tsv *EBV*.tsv *GL*.tsv NT_*.tsv || echo "not all files found"
cd ../../""".format(
            homerTagDir=homerTagDir,
            output_dir=output_dir
        )
    )

def makeTagDir(
    output_dir,
    input_bam,
    genome,
    restriction_site=None,
    illuminaPE=False,
    other_options=None
    ):
    ## if hic experiment then input_bam should be provided twice: {input_bam},{input_bam}
    input_bam_file = input_bam.split(",")[0]
    outputs=[
        os.path.join(output_dir, "tagInfo.txt"),
        os.path.join(output_dir, "tagGCcontent.txt"),
        os.path.join(output_dir, "genomeGCcontent.txt"),
        os.path.join(output_dir, "tagLengthDistribution.txt")
    ]

    return Job(
        input_files=[
            input_bam_file
        ],
        output_files=outputs,
        module_entries=[
            ["homer_tag_directory", "module_perl"],
            ["homer_tag_directory", "module_homer"],
            ["homer_tag_directory", "module_samtools"]
        ],
        command="""\
makeTagDirectory \\
  {output_dir} \\
  {input_bam} \\
  -genome {genome} \\
  -checkGC{illuminaPE}{restriction_site}{other_options}""".format(
            output_dir=output_dir,
            input_bam=input_bam,
            genome=genome,
            illuminaPE=" \\\n -illuminaPE" if illuminaPE else "",
            restriction_site=" \\\n -restrictionSite " + restriction_site if restriction_site else "",
            other_options=" \\\n " + other_options if other_options else ""
        ),
        name="homer_make_tag_directory",
        report_files=outputs
    )

##### Homer for HiC:

def hic_tagDirQcPlots(
    name,
    working_dir,
    chrlist,
    output_dir="HomerQcPlots"
    ):

    return Job(
        input_files=[
            os.path.join(working_dir, "tagInfo.txt")
        ],
        output_files=[
            os.path.join(working_dir, output_dir)
        ],
        module_entries=[
            ["homer_tag_directory", "module_R"],
            ["homer_tag_directory", "module_mugqic_tools"]
        ],
        command="""\
Rscript $R_TOOLS/HomerHiCQcPlotGenerator.R \\
  {name} \\
  {working_dir} \\
  {chrlist} \\
  {output_dir}""".format(
            name=name,
            working_dir=working_dir,
            chrlist=chrlist,
            output_dir=output_dir
        )
    )

def hic_interactionMatrix_chr(
    homer_dir,
    res,
    chr,
    fileName,
    fileNameRN=None,
    norm="raw",
    format=True
    ):

    # norm can be "raw", "simpleNorm", "norm"

    if fileNameRN is None:
        fileNameRN = re.sub(r"\.txt", "RN.txt", fileName)

    commandFormatMatrix = None
    if format:
        commandFormatMatrix = """\
cut -f 2- \\
  {fileName} \\
  > {fileNameRN}""".format(
            fileName=fileName,
            fileNameRN=fileNameRN
        )

    return Job(
        input_files=[
            os.path.join(homer_dir, "tagInfo.txt")
        ],
        output_files=[
            fileNameRN,
            fileName
        ],
        module_entries=[
            ["interaction_matrices_Chr", "module_perl"],
            ["interaction_matrices_Chr", "module_homer"]
        ],
        command="""\
analyzeHiC {homer_dir} \\
  -res {res} \\
  -{norm} \\
  -chr {chr} \\
  > {fileName}{extra_cmd}""".format(
            homer_dir=homer_dir,
            res=res,
            norm=norm,
            chr=chr,
            fileName=fileName,
            fileNameRN=fileNameRN,
            extra_cmd=f" && {commandFormatMatrix}" if commandFormatMatrix else ""
        ),
        removable_files=[fileName]
    )

def hic_interactionMatrix_genome(
    homer_dir,
    res,
    fileName,
    fileNameRN=None,
    norm="raw",
    format = True
    ):
    # norm can be "raw", "simpleNorm", "norm"

    if fileNameRN is None:
        fileNameRN = re.sub(r"\.txt", "RN.txt", fileName)
 
    commandFormatMatrix = None
    if format:
        commandFormatMatrix = """\
cut -f 2- {fileName} \\
  > {fileNameRN}""".format(
            fileName=fileName,
            fileNameRN=fileNameRN
        )

    return Job(
        input_files=[
            os.path.join(homer_dir, "tagInfo.txt")
        ],
        output_files=[
            fileNameRN,
            fileName
        ],
        module_entries=[
            ["interaction_matrices_genome", "module_perl"],
            ["interaction_matrices_genome", "module_homer"]
        ],
        command="""\
analyzeHiC \\
  {homer_dir} \\
  -res {res} \\
  -{norm} \\
  > {fileName}{format_cmd}""".format(
            homer_dir=homer_dir,
            res=res,
            norm=norm,
            fileName=fileName,
            format_cmd=f"&& {commandFormatMatrix}" if commandFormatMatrix else ""
        ),
        removable_files=[fileName]
    )

def hic_compartments(
    fileName,
    homer_dir,
    res,
    genome,
    fileName_PC1,
    fileName_Comp
    ):
    return Job(
        input_files=[
            os.path.join(homer_dir, "tagInfo.txt")
        ],
        output_files=[
            fileName_PC1,
            fileName_Comp
        ],
        module_entries=[
            ["identify_compartments", "module_perl"],
            ["identify_compartments", "module_homer"],
            ["identify_compartments", "module_R"]
        ],
        command="""\
runHiCpca.pl \
  {fileName} \\
  {homer_dir} \\
  -res {res} \\
  -genome {genome} \\
  -cpu {cpu} && \\
findHiCCompartments.pl \\
  {fileName_PC1} \\
  > {fileName_Comp}""".format(
            fileName=fileName,
            homer_dir=homer_dir,
            res=res,
            genome=genome,
            cpu=global_conf.global_get("identify_compartments", "threads", param_type='posint'),
            fileName_PC1=fileName_PC1,
            fileName_Comp=fileName_Comp
        ),
        removable_files=[fileName]
    )

def hic_peaks(
    homer_dir,
    res,
    genome,
    fileName,
    fileName_anno
    ):
    return Job(
        input_files=[
            os.path.join(homer_dir, "tagInfo.txt")
        ],
        output_files=[
            fileName,
            fileName_anno
        ],
        module_entries=[
            ["identify_peaks", "module_perl"],
            ["identify_peaks", "module_homer"]
        ],
        command="""\
findHiCInteractionsByChr.pl \\
  {homer_dir} \\
  -res {res} \\
  -cpu {cpu} \\
  > {fileName} && \\
annotateInteractions.pl \\
  {fileName} \\
  {genome} \\
  {fileName_anno}""".format(
            homer_dir=homer_dir,
            res=res,
            fileName=fileName,
            genome=genome,
            fileName_anno=fileName_anno,
            cpu=global_conf.global_get("identify_peaks", "threads", param_type='posint')
        )
    )

##### Homer for ChIP-Seq:
def makeUCSCfile(
    tag_dir,
    bedgraph_file
    ):
    bedgraph_file_gz = bedgraph_file + ".gz"
    return Job(
        input_files=[
            os.path.join(tag_dir, "tagInfo.txt")
        ],
        output_files=[
            bedgraph_file,
            bedgraph_file_gz
        ],
        module_entries=[
            ["homer_make_ucsc_files", "module_perl"],
            ["homer_make_ucsc_files", "module_homer"]
        ],
        command="""\
makeUCSCfile \\
  {tag_dir} \\
  > {bedgraph_file} && \\
gzip -c {bedgraph_file} \\
  > {bedgraph_file_gz}""".format(
            tag_dir=tag_dir,
            bedgraph_file=bedgraph_file,
            bedgraph_file_gz=bedgraph_file_gz
        ),
        name="homer_make_ucsc_file"
    )

def annotatePeaks(
    peak_file,
    genome,
    output_dir,
    annotation_file,
    genome_size
    ):
    return Job(
        input_files=[
            peak_file
        ],
        output_files=[
            annotation_file,
            os.path.join(output_dir, "geneOntology.html"),
            os.path.join(output_dir, "GenomeOntology.html")
        ],
        module_entries=[
            ["homer_annotate_peaks", "module_perl"],
            ["homer_annotate_peaks", "module_homer"]
        ],
        command="""\
annotatePeaks.pl \\
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
            annotation_file=annotation_file
        ),
        name="homer_make_ucsc_file"
    )

def findMotifsGenome(
    peak_file,
    genome,
    output_dir
    ):
    return Job(
        input_files=[
            peak_file
        ],
        output_files=[
            os.path.join(output_dir, "homerResults.html"),
            os.path.join(output_dir, "knownResults.html")
        ],
        module_entries=[
            ["homer_find_motifs_genome", "module_perl"],
            ["homer_find_motifs_genome", "module_homer"],
            ["homer_find_motifs_genome", "module_weblogo"]
        ],
        command="""\
findMotifsGenome.pl \\
  {peak_file} \\
  {genome} \\
  {output_dir} \\
  -preparsedDir {output_dir}/preparsed \\
  -p {threads}""".format(
            peak_file=peak_file,
            genome=genome,
            output_dir=output_dir,
            threads=global_conf.global_get('homer_find_motifs_genome', 'threads', param_type='posint')
        ),
        name="homer_find_motifs_genome"
    )
