#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def gemini_annotations(variants, gemini_output, tmp_dir):

    return Job(
        [variants],
        [gemini_output],
        [
            ['gemini_annotations', 'module_gemini'],
            ['gemini_annotations', 'module_htslib']
        ],
        command="""\
gemini load -v {variants} \\
  {options} \\
  --tempdir {temp} \\
  {output}""".format(
        options=config.param('gemini_annotations', 'options'),
        variants=variants,
        output=gemini_output,
        temp=tmp_dir
        )
    )
