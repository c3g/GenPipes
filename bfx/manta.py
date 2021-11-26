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


from core.config import *
from core.job import *

def manta_config(input_normal, input_tumor, output_dir, callRegion=None):
    output = [os.path.join(output_dir, "results", "variants", "diploidSV.vcf.gz"), os.path.join(output_dir, "results", "variants", "candidateSmallIndels.vcf.gz")]
    if input_tumor is not None:
        inputs = [input_normal, input_tumor, callRegion]
    else:
        inputs = [input_normal, callRegion]
        
    return Job(
        inputs,
        [output_dir],
        [
            ['manta_sv', 'module_python'],
            ['manta_sv', 'module_manta']
        ],
        command="""\
 python $MANTA_HOME/bin/configManta.py \\
        --normalBam {normal} \\
        {tumor} \\
        --referenceFasta {genome} \\
        {experiment_type} {callRegion} \\
        --runDir {output}""".format(
            normal=input_normal,
            tumor="--tumorBam " + input_tumor if input_tumor else "",
            genome=config.param('manta_sv','genome_fasta',type='filepath'),
            experiment_type=config.param('manta_sv','experiment_type_option') if config.param('manta_sv','experiment_type_option') else "",
            callRegion="\\\n        --callRegions " + callRegion if callRegion else "",
            output=output_dir
        )
    )

def manta_run(input_dir, output_dep):
    return Job(
        [input_dir],
        output_dep,
        [
            ['manta_sv', 'module_python'],
            ['manta_sv', 'module_manta']
        ],    
        command="""\
python {input_dir}/runWorkflow.py \\
        -m {mode}  \\
        -j {nodes} \\
        -g {ram} \\
        --quiet""".format(
            input_dir=input_dir,
            mode=config.param('manta_sv','option_mode'),
            nodes=config.param('manta_sv','option_nodes'),
            ram=config.param('manta_sv','ram')
        )
    )
