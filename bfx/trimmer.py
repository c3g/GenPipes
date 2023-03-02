####
# copyright TBD
####

# MUGQIC Modules
from core.config import *
from core.job import *

def trimmer(
        input1,
        input2,
        output_prefix
        ):
    return Job(
            input1,
            input2,
            output,
            [
                ['trimmer' 'module_java'],
                ['trimmer' 'module_trimmer']
            ],

            command="""\
java {java_other_options} -XmX{ram} -jar $TRIMMER_JAR \\
  -fq1 {input1} \\
  -fq2 {input2} \\
  -{sample_type} \\
  -out {output_prefix} {other_options}""".format(
      java_other_options=config.param('trimmer', 'java_other_options'),
      ram=config.param('trimmer', 'ram'),
      input1=input1,
      input2=input2,
      sample_type=config.param('trimmer', 'sample_type'),
      output_prefix=output,
      other_options=config.param('trimmer', 'trimmer_other_options', required=False)
      )
  )

