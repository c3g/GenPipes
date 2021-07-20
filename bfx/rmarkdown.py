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

# MUGQIC Modules
from core.config import *
from core.job import *

# This is general to all jobs 
def render(job_input, job_name, input_rmarkdown_file, samples, render_output_dir, module_section="DEFAULT", prerun_r=""):
    # input_rmarkdown_file='/root/blu/awe.some.Rmd';render_output_dir='report'
    output_markdown_file = os.path.join(render_output_dir, os.path.splitext( os.path.basename(input_rmarkdown_file) )[0] + '.md')

    if not isinstance(job_input, list):
        job_input = [job_input]

    return Job(
        job_input, [output_markdown_file],
        [
            [module_section, 'module_R'],
            [module_section, 'module_pandoc']
        ],
        command="""\
R --no-save --no-restore <<-'EOF'
{prerun_r}
input_rmarkdown_file = '{input_rmarkdown_file}'
render_output_dir    = '{render_output_dir}'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF""".format(input_rmarkdown_file=input_rmarkdown_file, render_output_dir=render_output_dir, prerun_r=prerun_r ),
        name=job_name,
        samples=samples,
        report_files=[output_markdown_file],
    )
