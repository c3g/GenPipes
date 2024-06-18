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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def quant(
    inputs,
    output_dir,
    transcriptome,
    parameters
    ):
    """
    Call to kalliso quant
    """

    return Job(
        inputs,
        [
            os.path.join(output_dir, "transcripts.tsv"),
            os.path.join(output_dir, "kallisto_quant.log"),
            os.path.join(output_dir, "abundance.h5")
            ],
        [
            ['kallisto', 'module_kallisto']
        ],
        command="""\
kallisto quant \\
  {parameters} \\
  -t {threads} \\
  -i {transcriptome} \\
  -o {output_dir} \\
  {infiles} \\
  2> {output_dir}/kallisto_quant.log""".format(
            parameters=parameters,
            threads=global_conf.global_get('kallisto', 'threads'),
            transcriptome=transcriptome,
            output_dir=output_dir,
            infiles=" \\\n  ".join(inputs),
        )
    )

def parse_mean_insert_size_metrics_pt(input_file):
    """
    Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
module load mugqic/python/3.10.4 &&
export mean_insert_size=`python <<EOF
import h5py
import numpy as np
with h5py.File('{input_file}', 'r') as h5file:
    frequencies = np.asarray(h5file['aux']['fld'])
    values = np.arange(0, len(frequencies), 1)
    print(np.around(np.average(values, weights=frequencies), decimals=1))
EOF`"""
        )

def parse_median_insert_size_metrics_pt(input_file):
    """
    Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
module load mugqic/python/3.10.4 &&
export median_insert_size=`python <<EOF
import h5py
import numpy as np
with h5py.File('{input_file}', 'r') as h5file:
    frequencies = np.asarray(h5file['aux']['fld'])
    values = np.arange(0, len(frequencies), 1)
    ord = np.argsort(values)
    cdf = np.cumsum(frequencies[ord])
    print(values[ord][np.searchsorted(cdf, cdf[-1] // 2)].astype('float64'))
EOF`"""
        )
