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

# MUGQIC Modules
from core.config import *
from core.job import *

def paired_somatic(input, normal_name, tumor_name, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
zless {input} | \\
        vawk --header \\
        '(S${tumor_name}$GT!="0/0" && S${tumor_name}$GT!="./." \\
        && S${tumor_name}$GT!=S${normal_name}$GT) \\
        && (S${normal_name}$GT=="0/0" || S${normal_name}$GT=="./.")' \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            tumor_name=tumor_name,
            output="> " + output if output else "",
        )
    )

def paired_germline(input, normal_name, tumor_name, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
zless {input} | \\
        vawk --header \\
        '(S${normal_name}$GT!="0/0" && S${normal_name}$GT!="./." && S${tumor_name}$GT!="./.")' \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            tumor_name=tumor_name,
            output="> " + output if output else "",
        )
    )

def single_germline(input, normal_name, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
zless {input} | \\
        vawk --header \\
        '(S${normal_name}$GT!="0/0" && S${normal_name}$GT!="./.")' \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            output="> " + output if output else "",
        )
    )

def sv(input, normal_name, tumor_name, caller, output):
    return Job(
        [os.path.abspath(input)],
        [os.path.abspath(output)],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
vawk -v CALLER={caller} -v SNAME={tumor_name} \\
        '{{if (($7=="PASS" || $7 == ".") && (S${tumor_name}$GT!="0/0" && S${tumor_name}$GT!="./." && S${normal_name}$GT!="./.")) \\
        print CALLER,I$SV_HIGHEST_TIER,SNAME,$1,$2,I$END,I$SVTYPE=="BND" ? I$SVTYPE":"$3":"I$MATEID":"$5 : I$SVTYPE, I$LOF, I$SIMPLE_ANN, I$PE, I$SR, \\
        S${tumor_name}$SR, S${normal_name}$SR, S${tumor_name}$PE, S${normal_name}$PE, S${tumor_name}$PR, S${normal_name}$PR, \\
        S${tumor_name}$RS, S${tumor_name}$AS, S${normal_name}$RS, S${normal_name}$AS, S${tumor_name}$DV, S${normal_name}$DV }}' \\
        {input} \\
        {output}""".format(
            input=input,
            normal_name=normal_name,
            tumor_name=tumor_name,
            caller=caller,
            output="> " + output if output else "",
        )
    )

def format_cna(input, output):
    return Job(
        [input],
        [output],
        [
            ['vawk', 'module_vawk'],
        ],
        command="""\
cat <(echo "Chromosome;Start;End;Segment_Mean" | tr ';' '\t') \\
<(zless {input} | \\
        vawk \\
        '{{print $1, $2, I$END, I$LOG_FOLD_CHANGE}}') \\
        {output}""".format(
            input=input,
            output="> " + output if output else "",
        )
    )
