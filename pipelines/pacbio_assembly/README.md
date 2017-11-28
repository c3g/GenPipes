usage: pacbio_assembly.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                          [-o OUTPUT_DIR] [-j {pbs,batch,daemon}] [-f]
                          [--report] [--clean]
                          [-l {debug,info,warning,error,critical}]
                          [-r READSETS] [-v]

Version: 3.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch,daemon}, --job-scheduler {pbs,batch,daemon}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- smrtanalysis_filtering
2- pacbio_tools_get_cutoff
3- preassembly
4- assembly
5- polishing
6- pacbio_tools_assembly_stats
7- blast
8- mummer
9- compile
10- circlator
