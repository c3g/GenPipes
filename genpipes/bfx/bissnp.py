################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def bisulfite_genotyper(input, cpg_output, snp_output):

    return Job(
        [input],
        [cpg_output, snp_output],
        [
            ['bissnp', 'module_java'],
            ['bissnp', 'module_bissnp']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $BISSNP_JAR \\
  --analysis_type BisulfiteGenotyper \\
  --reference_sequence {reference_genome} \\
  --input_file {input} \\
  --dbsnp {dbsnp} \\
  --vcf_file_name_1 {cpg_output} \\
  --vcf_file_name_2 {snp_output} \\
  --num_threads {threads}""".format(
            tmp_dir=global_conf.global_get('bissnp', 'tmp_dir'),
            java_other_options=global_conf.global_get('bissnp', 'java_other_options'),
            ram=global_conf.global_get('bissnp', 'ram'),
            reference_genome=global_conf.global_get('bissnp', 'genome_fasta'),
            input=input,
            dbsnp=global_conf.global_get('bissnp', 'known_variants'),
            cpg_output=cpg_output,
            snp_output=snp_output,
            threads=global_conf.global_get('bissnp', 'threads'),
        )
    )
