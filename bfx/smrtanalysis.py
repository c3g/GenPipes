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
from core.config import global_config_parser
from core.job import Job

def blasr(
    infile,
    infile_long,
    outfile,
    outfile_fofn,
    sam=False
    ):

    outfile_filtered = outfile + ".filtered"

    return Job(
        [infile, infile_long],
        [outfile_filtered, outfile_fofn],
        [
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
blasr \\
  {infile} \\
  {infile_long} \\
  -out {outfile} \\
  -m {m} \\
  -nproc {threads} \\
  -bestn {bestn} \\
  -nCandidates {n_candidates} \\
  -noSplitSubreads \\
  -minReadLength {min_read_length} \\
  -maxScore {max_score} \\
  -maxLCPLength {max_lcp_length}{sam} && \\
echo {outfile} > {outfile_fofn} &&
filterm4.py {outfile} > {outfile_filtered} 2> {outfile_filtered}.log'""".format(
            infile=infile,
            infile_long=infile_long,
            outfile=outfile,
            m=global_config_parser.param('smrtanalysis_blasr', 'm', param_type='int'),
            threads=global_config_parser.param('smrtanalysis_blasr', 'threads', param_type='posint'),
            bestn=global_config_parser.param('smrtanalysis_blasr', 'bestn', param_type='int'),
            n_candidates=global_config_parser.param('smrtanalysis_blasr', 'n_candidates', param_type='int'),
            min_read_length=global_config_parser.param('smrtanalysis_blasr', 'min_read_length', param_type='int'),
            max_score=global_config_parser.param('smrtanalysis_blasr', 'max_score', param_type='int'),
            max_lcp_length=global_config_parser.param('smrtanalysis_blasr', 'max_lcp_length', param_type='int'),
            sam=" \\\n  -sam" if sam else "",
            outfile_fofn=outfile_fofn,
            outfile_filtered=outfile_filtered
    ))

def cmph5tools_sort(
    cmph5,
    cmph5_out
    ):

    return Job(
        [cmph5],
        [cmph5_out],
        [
            ['smrtanalysis_cmph5tools_sort', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
cmph5tools.py -vv sort --deep --inPlace \\
  --outFile {cmph5_out} \\
  {cmph5} \\
  > /dev/null'""".format(
            cmph5=cmph5,
            cmph5_out=cmph5_out
    ))

def fastq_to_ca(
    libraryname,
    reads,
    outfile
    ):

    return Job(
        [reads],
        [outfile],
        [
            ['smrtanalysis_fastq_to_ca', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
fastqToCA \\
  -technology sanger \\
  -type sanger \\
  -libraryname {libraryname} \\
  -reads {reads} \\
  > {outfile}'""".format(
            libraryname=libraryname,
            reads=reads,
            outfile=outfile)
    )


def filtering(fofn, input_xml, params_xml, output_dir, log):

    ref_params_xml = global_config_parser.param('smrtanalysis_filtering', 'filtering_settings')

    output_prefix = os.path.join(output_dir, "data", "filtered_subreads.")
    input_fofn = os.path.join(output_dir, "input.fofn")
    output_fastq = output_prefix + "fastq"

    white_path = global_config_parser.param('smrtanalysis_filtering', 'whitelist_path', required=False)

    temp_dir = global_config_parser.param('smrtanalysis_filtering', 'tmp_dir')
    if (global_config_parser.param('smrtanalysis_filtering', 'tmp_dir').startswith("$")):
        temp_dir = temp_dir.replace('$', '$SMRT_ORIGUSERENV_')
        temp_dir = temp_dir.replace('{', '')
        temp_dir = temp_dir.replace('}', '')

    return Job(
        [fofn, ref_params_xml],
        [input_fofn, output_fastq, output_prefix + "fasta", os.path.join(output_dir, "data", "filtered_summary.csv")],
        [
            ['smrtanalysis_filtering', 'module_prinseq'],
            ['smrtanalysis_filtering', 'module_smrtanalysis']
        ],
        command="""\
bash -ic 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
fofnToSmrtpipeInput.py {fofn} > {input_xml} && \\
cp {fofn} {input_fofn} && \\
sed -e "s|MINSUBREADLENGTH|{min_subread_length}|g" -e "s|MINREADLENGTH|{min_read_length}|g" -e "s|MINQUAL|{min_qual}|g"{whitelist_param} \\
  < {ref_params_xml} > {params_xml} && \\
smrtpipe.py \\
  -D NPROC={threads} \\
  -D TMP={tmp_dir} \\
  --params={params_xml} \\
  --output={output_dir} \\
  --debug \\
  xml:{input_xml} \\
  > {log}' && \\
prinseq-lite.pl \\
  -verbose \\
  -fastq {output_fastq} \\
  -out_format 1 \\
  -out_good {output_dir}/data/filtered_subreads""".format(
            fofn=fofn,
            input_fofn=input_fofn,
            input_xml=input_xml,
            min_subread_length=global_config_parser.param('smrtanalysis_filtering', 'min_subread_length'),
            min_read_length=global_config_parser.param('smrtanalysis_filtering', 'min_read_length'),
            min_qual=global_config_parser.param('smrtanalysis_filtering', 'min_qual'),
            whitelist_param=' -e "s|<\!-- WHITELISTCOM||g" -e "s|WHITELISTCOM -->||g" -e "s|WHITELISTFILEPATH|'+white_path+'|g"' if white_path != "" else '',
            ref_params_xml=ref_params_xml,
            params_xml=params_xml,
            threads=global_config_parser.param('smrtanalysis_filtering', 'threads'),
            tmp_dir=temp_dir,
            output_dir=output_dir,
            log=log,
            output_fastq=output_fastq
    ))

def load_chemistry(
    cmph5,
    input_fofn,
    cmph5_output
    ):

    return Job(
        [input_fofn, cmph5],
        [cmph5_output],
        [
            ['smrtanalysis_load_chemistry', 'module_smrtanalysis']
        ],
        # Copy cmph5 file since loadChemistry.py modifies the input cmph5 directly
        command="""\
cp {cmph5} {cmph5_output} && \\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
loadChemistry.py \\
  {input_fofn} \\
  {cmph5_output}'""".format(
            input_fofn=input_fofn,
            cmph5=cmph5,
            cmph5_output=cmph5_output
    ))

def load_pulses(
    cmph5,
    input_fofn,
    cmph5_output
    ):

    return Job(
        [input_fofn, cmph5],
        [cmph5_output],
        [
            ['smrtanalysis_load_pulses', 'module_smrtanalysis']
        ],
        # Copy cmph5 file since loadPulses modifies the input cmph5 directly
        command="""\
cp {cmph5} {cmph5_output} && \\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
loadPulses \\
  {input_fofn} \\
  {cmph5_output} \\
  -metrics DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag -byread'""".format(
            input_fofn=input_fofn,
            cmph5=cmph5,
            cmph5_output=cmph5_output
    ))

def m4topre(
    infile,
    allm4,
    subreads,
    outfile
    ):

    return Job(
        [infile],
        [outfile],
        [
            ['smrtanalysis_m4topre', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
m4topre.py \\
  {infile} \\
  {allm4} \\
  {subreads} \\
  {bestn} \\
  > {outfile}'""".format(
            infile=infile,
            allm4=allm4,
            subreads=subreads,
            bestn=global_config_parser.param('smrtanalysis_m4topre', 'bestn', param_type='int'),
            outfile=outfile
    ))

def pbalign(
    polishing_round_directory,
    sample_name,
    sample_cutoff_mer_size_polishing_round
    ):

    input_fofn = os.path.join(sample_name, "filtering", "input.fofn")
    ref_upload = os.path.join(polishing_round_directory, sample_cutoff_mer_size_polishing_round, "sequence", sample_cutoff_mer_size_polishing_round + ".fasta")

    cmph5 = os.path.join(polishing_round_directory, "data", "aligned_reads.cmp.h5")

    control_regions_fofn = os.path.join(sample_name, "filtering", "data", "filtered_regions.fofn")

    temp_dir = os.path.join(global_config_parser.param('smrtanalysis_filtering', 'tmp_dir'), sample_cutoff_mer_size_polishing_round)
    if (global_config_parser.param('smrtanalysis_pbutgcns', 'tmp_dir').startswith("$")):
        temp_dir = temp_dir.replace('$', '$SMRT_ORIGUSERENV_')
        temp_dir = temp_dir.replace('{', '')
        temp_dir = temp_dir.replace('}', '')

    return Job(
        [input_fofn, ref_upload],
        [cmph5],
        [
            ['smrtanalysis_pbalign', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
pbalign \\
  {input_fofn} \\
  {ref_upload} \\
  {cmph5} \\
   --seed=1 --minAccuracy=0.75 --minLength=50 --algorithmOptions="-useQuality" --algorithmOptions=" -minMatch 12 -bestn 10 -minPctIdentity 70.0" --hitPolicy=randombest \\
  --tmpDir={tmp_dir} \\
  -vv \\
  --nproc={threads} \\
  --regionTable={control_regions_fofn}'""".format(
            input_fofn=input_fofn,
            ref_upload=ref_upload,
            cmph5=cmph5,
            tmp_dir=temp_dir,
            threads=global_config_parser.param('smrtanalysis_pbalign', 'threads', param_type='posint'),
            control_regions_fofn=control_regions_fofn
    ))

def pbdagcon(
    infile,
    outfile,
    outfile_fastq
    ):

    return Job(
        [infile],
        [outfile, outfile_fastq],
        [
            ['smrtanalysis_pbdagcon', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
pbdagcon -a -j {threads} \\
  {infile} \\
  > {outfile}' && \\
awk '{{if ($0~/>/) {{sub(/>/,"@",$0);print;}} else {{l=length($0);q=""; while (l--) {{q=q "9"}}printf("%s\\n+\\n%s\\n",$0,q)}}}}' {outfile} \\
  > {outfile_fastq}""".format(
            threads=global_config_parser.param('smrtanalysis_pbdagcon', 'threads', param_type='posint'),
            infile=infile,
            outfile=outfile,
            outfile_fastq=outfile_fastq
    ))


def pbutgcns(assembly_directory, sample_cutoff_mer_size, mer_size_directory):

    gpk_store = os.path.join(assembly_directory, sample_cutoff_mer_size + ".gkpStore")
    tig_store = os.path.join(assembly_directory, sample_cutoff_mer_size + ".tigStore")

    outfile = os.path.join(assembly_directory, "9-terminator", sample_cutoff_mer_size + ".ctg.fasta")

    unitigs_list = os.path.join(mer_size_directory, "unitigs.lst")
    prefix = os.path.join(assembly_directory, sample_cutoff_mer_size)
    outdir = os.path.join(assembly_directory, "9-terminator")

    temp_dir = os.path.join(global_config_parser.param('smrtanalysis_filtering', 'tmp_dir'), sample_cutoff_mer_size)
    if (global_config_parser.param('smrtanalysis_pbutgcns', 'tmp_dir').startswith("$")):
        temp_dir = temp_dir.replace('$', '$SMRT_ORIGUSERENV_')
        temp_dir = temp_dir.replace('{', '')
        temp_dir = temp_dir.replace('}', '')

    return Job(
        [gpk_store, tig_store],
        [outfile],
        [
            ['smrtanalysis_pbutgcns', 'module_smrtanalysis']
        ],
        command="""\
tigStore \\
  -g {gpk_store} \\
  -t {tig_store} 1 \\
  -d properties -U | \\
awk 'BEGIN{{t=0}}$1=="numFrags"{{if ($2 > 1) {{print t, $2}} t++}}' | sort -nrk2,2 \\
  > {unitigs_list} && \\
mkdir -p {outdir} && \\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && mkdir -p {tmp_dir} && \\
tmp={tmp_dir} \\
cap=$PWD/{prefix} \\
utg=$PWD/{unitigs_list} \\
nprocs={threads} \\
cns=$PWD/{outfile} \\
pbutgcns_wf.sh'""".format(
            gpk_store=gpk_store,
            tig_store=tig_store,
            unitigs_list=unitigs_list,
            outdir=outdir,
            tmp_dir=temp_dir,
            prefix=prefix,
            threads=global_config_parser.param('smrtanalysis_pbutgcns', 'threads', param_type='posint'),
            outfile=outfile
    ))

def reference_uploader(
    prefix,
    sample_name,
    fasta
    ):

    return Job(
        [fasta],
        [os.path.join(prefix, sample_name, "sequence", sample_name + ".fasta")],
        [
            ['smrtanalysis_reference_uploader', 'module_smrtanalysis']
        ],
        # Preload assembled contigs as reference
        command="""\
$(set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
referenceUploader \\
  --skipIndexUpdate \\
  --create \\
  --refRepos {prefix} \\
  --name {sample_name} \\
  --fastaFile {fasta} \\
  --saw="sawriter -blt 8 -welter" --jobId="Anonymous" \\
  --samIdx="samtools faidx" --jobId="Anonymous" --verbose)""".format(
            prefix=prefix,
            sample_name=sample_name,
            fasta=fasta
    ))

def run_ca(
    infile,
    ini,
    prefix,
    outdir
    ):

    return Job(
        [infile, ini],
        [
            os.path.join(outdir, prefix + ".ovlStore.list"),
            os.path.join(outdir, prefix + ".tigStore"),
            os.path.join(outdir, prefix + ".gkpStore")
        ],
        [
            ['smrtanalysis_run_ca', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
runCA \\
  -s {ini} \\
  -p {prefix} \\
  -d {outdir} \\
  {infile}'""".format(
            infile=infile,
            ini=ini,
            prefix=prefix,
            outdir=outdir
    ))

def summarize_polishing(
    sample_name,
    reference,
    aligned_reads_cmph5,
    alignment_summary,
    coverage_bed,
    input_fofn,
    aligned_reads_sam,
    variants_gff,
    variants_bed,
    variants_vcf
    ):

    return Job(
        [aligned_reads_cmph5, input_fofn, os.path.join(os.path.dirname(reference), "data", "consensus.fasta")],
        [alignment_summary, coverage_bed, aligned_reads_sam, variants_bed, variants_vcf],
        [
            ['smrtanalysis_run_ca', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
summarize_coverage.py \\
  --reference {reference} \\
  --numRegions=500 \\
  {aligned_reads_cmph5} \\
  {alignment_summary} && \\
gffToBed.py \\
  --name=meanCoverage \\
  --description="Mean coverage of genome in fixed interval regions" \\
  coverage {alignment_summary} \\
  > {coverage_bed} && \\
h5repack -f GZIP=1 \\
  {aligned_reads_cmph5} \\
  {aligned_reads_cmph5}.repacked && \\
pbsamtools --bam \\
  --outfile {aligned_reads_sam} \\
  --refrepos {reference} \\
  --readGroup movie {aligned_reads_cmph5}.repacked && \\
summarizeConsensus.py \\
  --variantsGff {variants_gff} \\
  {alignment_summary} \\
  --output {alignment_summary}.tmp && \\
mv {alignment_summary}.tmp {alignment_summary} && \\
gffToBed.py --name=variants \\
  --description="PacBio: snps, insertions, and deletions derived from consensus calls against reference" \\
  variants {variants_gff} \\
  > {variants_bed} && \\
gffToVcf.py \\
  --globalReference={sample_name} \\
  {variants_gff} \\
  > {variants_vcf}'""".format(
            reference=reference,
            aligned_reads_cmph5=aligned_reads_cmph5,
            alignment_summary=alignment_summary,
            coverage_bed=coverage_bed,
            input_fofn=input_fofn,
            aligned_reads_sam=aligned_reads_sam,
            variants_gff=variants_gff,
            variants_bed=variants_bed,
            sample_name=sample_name,
            variants_vcf=variants_vcf
    ))

def variant_caller(
    cmph5,
    ref_fasta,
    outfile_variants,
    outfile_fasta,
    outfile_fastq
    ):

    outfile_fasta_uncompressed = re.sub("\.(gz|gzip)$", "", outfile_fasta)

    return Job(
        [cmph5, ref_fasta],
        [outfile_variants, outfile_fasta, outfile_fastq, re.sub("\.gz$", "", outfile_fasta_uncompressed)],
        [
            ['smrtanalysis_variant_caller', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
variantCaller.py \\
  -P {protocol} \\
  -v \\
  -j {threads} \\
  --algorithm={algorithm} \\
  {cmph5} \\
  -r {ref_fasta} \\
  -o {outfile_variants} \\
  -o {outfile_fasta} \\
  -o {outfile_fastq} \\
  > /dev/null' && \\
gunzip -c \\
  {outfile_fasta} \\
  > {outfile_fasta_uncompressed}""".format(
            protocol=global_config_parser.param('smrtanalysis_variant_caller', 'protocol', param_type='dirpath'),
            threads=global_config_parser.param('smrtanalysis_variant_caller', 'threads', param_type='posint'),
            algorithm=global_config_parser.param('smrtanalysis_variant_caller', 'algorithm'),
            cmph5=cmph5,
            ref_fasta=ref_fasta,
            outfile_variants=outfile_variants,
            outfile_fasta=outfile_fasta,
            outfile_fastq=outfile_fastq,
            outfile_fasta_uncompressed=outfile_fasta_uncompressed
    ))


def basemodification(
    fasta_consensus,
    aligned_reads,
    basemodification_file_prefix,
    mer_size_directory,
    polishing_rounds
    ):

    return Job(
        [fasta_consensus, aligned_reads],
        [basemodification_file_prefix],
        [
            ['basemodification', 'module_smrtanalysis']
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
ipdSummary.py {aligned_reads} \\
  --reference {fasta_consensus} \\
  --identify m6A,m4C \\
  --methylFraction \\
  --paramsPath /cvmfs/soft.mugqic/CentOS6/software/smrtanalysis/smrtanalysis_2.3.0.140936.p5/analysis/etc/algorithm_parameters/2015-11/kineticsTools/ \\
  --numWorkers 12 \\
  --outfile {output_gff}'""".format(
            fasta_consensus=os.path.join(mer_size_directory, "polishing" + str(polishing_rounds - 1), "data", "consensus.fasta"),
            aligned_reads=os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "aligned_reads.sorted.cmp.h5"),
            output_gff=basemodification_file_prefix
    ))

def motifMaker(
    fasta_consensus,
    basemodification_file_prefix,
    mer_size_directory,
    polishing_rounds,
    motifMaker_file
    ):

    return Job(
        [fasta_consensus, basemodification_file_prefix],
        [motifMaker_file],
        [
                  ['motifMaker', 'module_smrtanalysis'],
                  ['motifMaker', 'module_java'],
        ],
        command="""\
bash -c 'set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
motifMaker.sh find \\
  -f {fasta_consensus} \\
  -g {output_gff} \\
  -o {output}'""".format(
            fasta_consensus=os.path.join(mer_size_directory, "polishing" + str(polishing_rounds - 1), "data", "consensus.fasta"),
            output_gff=basemodification_file_prefix + ".gff",
            output=motifMaker_file
    ))
