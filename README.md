Repo Policy
------------------

This repo holds perl libs, wrappers and scripts of several bioinformatics pipelines.

The main repo dir organization is:

mugqic_pipeline  -  root dir

root/lib       - dir containing all libraries of all pipelines. 

root/tool_shed - dir containing all tools used by the pipelines (usually written in languages such as bash, awk, python, R etc...)

root/pipelines - dir containing the pipelines itself. An addtional dir should be created in this dir with the pipeline name (root/pipeline/Pipeline_name)



Documentation:

Perl documentation on *.pm and *.pl files should (as much as possible) be created using POD (Pod::Usage). 

