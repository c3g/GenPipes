#####
# copyright
######

#MUGQIC Modules
from core.config import *
from core.job import *

def dedup(input,
            output,
            bed=None,
            mode=None):
    return Job(
            input,
            output,
            [
                ['locatit' 'module_java'],
                ['locatit' 'module_locatit']
            ],

            command="""\
java {java_other_options} -XmX{ram} -jar $LOCATIT_JAR \\
  {mode} \\
  {other_options} \\
  {bed} \\
  -o {output} \\
  {input}""".format(
      java_other_options=config.param('locatit', 'java_other_options'),
      ram=config.param('locatit', 'java_other_options'),
      mode="-"+mode if mode else "",
      other_options=config.param('locatit', 'locatit_other_options')
      bed="-l " + bed if bed else "",
      output=output,
      input=input
    ) 
)

