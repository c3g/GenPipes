################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

#MUGQIC Modules
from core.config import *
from core.job import *

def run_ichorCNA(
        input_tumor,
        input_normal,
        input_id,
        output_dir
        ):
    return Job(
            [input_tumor, input_normal],
            [output_dir],
            [
                ['run_ichorCNA', 'module_R'],
                ['run_ichorCNA', 'module_ichorCNA']
            ],

            command="""\
Rscript $ICHORCNA_HOME/runIchorCNA.R {other_ichorCNA_options} \\
    --id {input_id} \\
    --WIG {input_tumor} \\
    --NORMWIG {input_normal} \\
    --ploidy {ploidy} \\
    --normal {normal_contam} \\
    --maxCN {CN_states} \\
    --gcWig {gcWig} \\
    --mapWig {mapWig} \\
    --centromere {centromere} \\
    --includeHOMD {includeHOMD} \\
    --chrs {chromosomes} \\
    --chrTrain {chr_train} \\
    --estimateNormal {estimateNormal} \\
    --estimatePloidy {estimatePloidy} \\
    --estimateScPrevalence {estimateScPrevalence} \\
    --scStates {scStates} \\
    --txnE {txnE} \\
    --txnStrength {txnStrength} \\
    --outDir {output_dir}""".format(
        other_ichorCNA_options=config.param('run_ichorCNA', 'other_ichorCNA_options', required=False),
        input_id=input_id,
        input_tumor=input_tumor,
        input_normal=input_normal,
        ploidy=config.param('run_ichorCNA', 'ploidy'),
        normal_contam=config.param('run_ichorCNA', 'normal_contam'),
        CN_states=config.param('run_ichorCNA', 'CN_states'),
        gcWig=config.param('run_ichorCNA', 'gcWig', param_type='filepath'),
        mapWig=config.param('run_ichorCNA', 'mapWig', param_type='filepath'),
        centromere=config.param('run_ichorCNA', 'centromere', param_type='filepath'),
        includeHOMD=config.param('run_ichorCNA', 'includeHOMD'),
        chromosomes=config.param('run_ichorCNA', 'chromosomes'),
        chr_train=config.param('run_ichorCNA', 'chr_train'),
        estimateNormal=config.param('run_ichorCNA', 'estimateNormal'),
        estimatePloidy=config.param('run_ichorCNA', 'estimatePloidy'),
        estimateScPrevalence=config.param('run_ichorCNA', 'estimateScPrevalence'),
        scStates=config.param('run_ichorCNA', 'scStates'),
        txnE=config.param('run_ichorCNA', 'txnE'),
        txnStrength=config.param('run_ichorCNA', 'txnStrength'),
        output_dir=output_dir
        )
    )
