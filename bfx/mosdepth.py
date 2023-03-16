#####
# TBD
#####

# MUGQIC Modules
from core.config import config
from core.job import Job

def mosdepth(input, output_prefix, per_base=False, regions=None):
    outputs = [
            output_prefix + ".mosdepth.global.dist.txt",
            output_prefix + ".mosdepth.summary.txt",
            output_prefix + ".mosdepth.region.dist.txt"
            ]
    return Job(
            [input],
            outputs,
            [
                ['mosdepth', 'module_python'],
                ['mosdepth', 'module_mosdepth']
            ],
            command="""\
mosdepth {other_options} \\
{per_base} \\
{regions} \\
{chrom} \\
{output_prefix} \\
{input}""".format(
    other_options=config.param('mosdepth', 'other_options'),
    per_base="--no-per-base " if not per_base else "",
    regions="--by " + regions if regions else "",
    chrom="--chrom " + config.param('mosdepth', 'chrom') if config.param('mosdepth', 'chrom', required=False) else "",
    output_prefix=output_prefix,
    input=input
        )
    )
