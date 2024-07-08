################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def create(readset, output, fastqc=False):
    if fastqc == True:
        if readset.run_type =="PAIRED_END" and readset.adapter2:  # Paired end reads
            return Job(
                command="""\
`cat > {adapter_file} << END
>Adapter1\t{sequence1}
>Adapter2\t{sequence2}
END`""".format(
            adapter_file=output,
            sequence1=readset.adapter1,
            sequence2=readset.adapter2,
            )
        )
        elif readset.run_type =="SINGLE_END" and readset.adapter1:  # Single end reads
            return Job(
                command="""\
`cat > {adapter_file} << END
>Adapter1\t{sequence1}
END`""".format(
            adapter_file=output,
            sequence1=readset.adapter1,
            )
        )
    else:
        if readset.run_type =="PAIRED_END" and readset.adapter2:  # Paired end reads
            return Job(
                command="""\
`cat > {adapter_file} << END
>Adapter1\n{sequence1}
>Adapter2\n{sequence2}
END`""".format(
            adapter_file=output,
            sequence1=readset.adapter1,
            sequence2=readset.adapter2,
            )
        )
        elif readset.run_type =="SINGLE_END" and readset.adapter1:  # Single end reads
            return Job(
                command="""\
`cat > {adapter_file} << END
>Adapter1\n{sequence1}
END`""".format(
            adapter_file=output,
            sequence1=readset.adapter1,
            )
        )
