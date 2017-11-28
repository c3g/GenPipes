usage: illumina_run_processing.py [-h] [--help] [-c CONFIG [CONFIG ...]]
                                  [-s STEPS] [-o OUTPUT_DIR]
                                  [-j {pbs,batch,daemon}] [-f] [--report]
                                  [--clean]
                                  [-l {debug,info,warning,error,critical}]
                                  [-d RUN_DIR] [--lane LANE_NUMBER]
                                  [-r READSETS] [-i CASAVA_SHEET_FILE]
                                  [-x FIRST_INDEX] [-y LAST_INDEX]
                                  [-m NUMBER_OF_MISMATCHES] [-w] [-v]

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
  -d RUN_DIR, --run RUN_DIR
                        run directory
  --lane LANE_NUMBER    lane number
  -r READSETS, --readsets READSETS
                        nanuq readset file. The default file is
                        'run.nanuq.csv' in the output folder. Will be
                        automatically downloaded if not present.
  -i CASAVA_SHEET_FILE  illumina casava sheet. The default file is
                        'SampleSheet.nanuq.csv' in the output folder. Will be
                        automatically downloaded if not present
  -x FIRST_INDEX        first index base to use for demultiplexing
                        (inclusive). The index from the sample sheet will be
                        adjusted according to that value.
  -y LAST_INDEX         last index base to use for demultiplexing (inclusive)
  -m NUMBER_OF_MISMATCHES
                        number of index mistmaches allowed for demultiplexing
                        (default 1). Barcode collisions are always checked.
  -w, --force-download  force the download of the samples sheets (default:
                        false)
  -v, --version         show the version information and exit

Steps:
------
1- index
2- fastq
3- align
4- picard_mark_duplicates
5- metrics
6- blast
7- qc_graphs
8- md5
9- copy
10- end_copy_notification
