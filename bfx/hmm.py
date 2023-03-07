#####
# copyright TBD
#####

#MUGQIC Modules
from core.config import *
from core.job import *

def readCounter(
        input,
        output
        ):
    return Job(
            [input],
            [output],
            [
#                ['hmm_readCounter', 'module_hmm']
            ],

            command="""\
$HMM_HOME/readCounter \\
    --window {window_size} \\
    --quality {threshold} \\
    --chromosome {chr_list} \\
    {input} \\
    > {output}""".format(
        window_size=config.param('hmm_readCounter', 'window_size'),
        threshold=config.param('hmm_readCounter', 'threshold'),
        chr_list=config.param('hmm_readCounter', 'chr_list'),
        input=input,
        output=output
        )
    )
