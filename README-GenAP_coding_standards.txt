# Here are the coding standards used by the GenAP/C3G developpers
# These should be followed by anyone willing to participate in the GenAP/C3G applications development

# General
Code should be well commented, clear and readable.

# Indentation :
Use 4 space indentations instead of tabs.

# Wrappers vs. command-lines
Useful methods should be written in common script files (e.g. as a python script within the bfx/ folder, or as a new method within pipeline/common.py) : no command-line should be hard-coded within a pipeline python wrapper (e.g. within dnaseq.py).
When a command has to be implemented, the developper has to wrap it as a method in a python script (e.g. pipelines/common.py or bfx/picard.py).
Methods should be made as general as possible with the possibility of adding command line options, with the uses of the configuration .ini files, in order to make them usable across pipelines.

# Command-line format
When coding commands, a single command should be divided across different lines (\) with every parameter on a line to enhance code readability. Also, each line should be preceded by 2 space characters, except for the first line of the command.
e.g.
motifMaker.sh find \\
  -f {fasta_consensus} \\
  -g {output_gff} \\
  -o {output}

# Use of modules
All software calls should be done using the module system. Modules should be specified within the .ini file(s). Preferably use the mugqic modules which hare distributed among the Compute Canada HPC clusters, this ensures stability of the analysis and reproductibily of the results. If using a non-mugqic module, ensure that the module is available on all the cluster. Avoid system modules which are not be transferable across different platforms.

# Pipeline development
JIRA tickets are usually created for each pipeline development, where benchmarking parameters and other decisions made can be logged for later reference.
Code is always developed on a branch, i.e. try as much as possible to not edit the master branch directly. Once code is ready to be merged, a pull request is created. After being reviewed by at least two developers on the team (if possible, developpers have to be uninvolved in the reviewed project), and after all the reviewer comments and enquiries have been adressed and accepted, the development branch can then be merged to master.

# Transferability
Pipelines should be transferable to a host of other clusters; so keep that in mind as you write your code.
