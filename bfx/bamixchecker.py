
    return Job(
        inputfiles,
        outputfiles,
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_java'],
            [ini_section, 'module_bedtools'],
            [ini_section, 'module_R'],
            [ini_section, 'module_pandoc'],
            [ini_section, 'module_gatk'],
            [ini_section, 'module_bamixchecker']
        ],
        command="""\
python $BAMixChecker_PATH/BAMixChecker.py \\
  -l {inputlist} \\
  -r {reference} \\
  -o {out_dir} \\
  {options}""".format(
            inputlist=inputlist,
            out_dir=out_dir,
            options=options,
            reference=reference
        )
    )