#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys

class ValidationError(Exception):
    pass

def load_file(file_path):
    '''
    Loads file and returns a list of dictionaries, also checks for extra tabs
    '''
    with open(file_path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        data = []
        for row_num, row in enumerate(reader, start=1):
            if len(row) > len(header):
                print(f"Warning: Extra tab(s) detected at the end of line {row_num + 1}. Please remove any extra tabs.")
            data.append(dict(zip(header, row)))
        return data


def unicode_check(string_array):
    '''
    Checks if there are non-ascii characters in file
    '''
    for row in string_array:
        try:
            [s.encode(encoding='utf-8').decode('ascii') for s in row.values() if s is not None]
        except UnicodeDecodeError:
            raise ValidationError(f"File contains non-English characters in line:\n{row}")
    return True

def row_size_check(string_array):
    '''
    Checks that all rows in the file have the same size
    '''
    columns = [len(row) for row in string_array]
    if not all(columns[0] == n for n in columns):
        raise ValidationError("Rows are not all the same size! Please make all the rows of equal size!")
    return True

def trailing_spaces_check(string_array):
    '''
    Checks that words don't have spaces around them
    '''
    for row in string_array:
        for s in row.values():
            if s is not None and re.search(r"[ \t]", s):
                raise ValidationError(f"There are trailing spaces around the following word: {s} ; Please fix it and try again!")
    return True

def dos_newline_check(string_array):
    '''
    Checks if the file has carriage returns and DOS newline characters
    '''
    for row in string_array:
        for s in row.values():
            if s is not None and re.search(r"[\r]", s):
                raise ValidationError(f"There is a carriage return after the following word: {s} ; To fix it, you can run 'dos2unix myfile'. Please fix it and try again!")
    return True

def readset_header_check(header, pipeline):
    '''
    Checks that the readset file header contains the required columns based on the pipeline
    '''
    valid_headers = {
        'ampliconseq': ['format5'],
        'chipseq': ['format2'],
        'covseq': ['format1'],
        'dnaseq': ['format1'],
        'methylseq': ['format1'],
        'longread_dnaseq': ['format4'],
        'nanopore_covseq': ['format3'],
        'rnaseq': ['format1'],
        'rnaseq_denovo_assembly': ['format1'],
        'rnaseq_light': ['format1'],
    }

    formats = {
        'format1': {
            'required': ['Sample', 'Readset', 'Library', 'RunType', 'Adapter1', 'Adapter2'],
            'optional': ['Run', 'Lane', 'QualityOffset', 'BED', 'FASTQ1', 'FASTQ2', 'BAM']
        },
        'format2': {
            'required': ['Sample', 'Readset', 'MarkName', 'MarkType', 'Library', 'Adapter1', 'Adapter2'],
            'optional': ['Run', 'Lane', 'QualityOffset', 'BED', 'FASTQ1', 'FASTQ2', 'BAM']
        },
        'format3': {
            'required': ['Sample', 'Readset', 'Flowcell', 'Library', 'FAST5'],
            'optional': ['Run', 'Summary', 'FASTQ', 'Barcode', 'AnalysisName']
        },
        'format4': {
            'required': ['Sample', 'Readset', 'Flowcell', 'Library'],
            'optional': ['Run', 'Summary', 'FASTQ', 'FAST5', 'BAM']
        },
        'format5': {
            'required': ['Sample', 'Readset', 'RunType', 'Adapter1', 'Adapter2', 'primer1', 'primer2'],
            'optional': ['FASTQ1', 'FASTQ2']
        }
    }

    pipeline_format = valid_headers.get(pipeline)
    if not pipeline_format:
        raise ValidationError(f"Pipeline {pipeline} is not recognized. Please revise and try again!")

    expected_header = formats[pipeline_format[0]]['required']
    optional_header = formats[pipeline_format[0]]['optional']

    missing_columns = [column for column in expected_header if column not in header]
    if missing_columns:
        raise ValidationError(f"File header is not correct for {pipeline}. It should contain the following required columns:\n\n{'\t'.join(expected_header)}\n\nMissing columns: {', '.join(missing_columns)}\n\nPlease fix this and try again!")

    return expected_header + optional_header

def check_column_dependencies(row, pipeline, row_num):
    '''
    Checks dependencies between columns based on the pipeline
    '''
    errors = []

    # Check mandatory columns for all pipelines
    if not row.get('Sample'):
        errors.append(f"Row {row_num}: Sample must be provided.")
    if not row.get('Readset'):
        errors.append(f"Row {row_num}: Readset must be provided.")

    # Check pipeline-specific dependencies
    if pipeline in ['covseq', 'dnaseq', 'methylseq', 'rnaseq', 'rnaseq_denovo_assembly', 'rnaseq_light']:
        if not row.get('Run'):
            errors.append(f"Row {row_num}: Run must be provided.")
        if not row.get('Lane'):
            errors.append(f"Row {row_num}: Lane must be provided.")
        if not row.get('FASTQ1') and not row.get('BAM'):
            errors.append(f"Row {row_num}: Either FASTQ1 or BAM must be provided.")
        if row.get('RunType') == 'PAIRED_END' and not (row.get('FASTQ2') or row.get('BAM')):
            errors.append(f"Row {row_num}: FASTQ2 must be provided for PAIRED_END RunType.")
        if row.get('FASTQ1') and row.get('BAM'):
            errors.append(f"Row {row_num}: BAM should be ignored if FASTQ1 is provided.")
    elif pipeline == 'chipseq':
        if not row.get('Run'):
            errors.append(f"Row {row_num}: Run must be provided.")
        if not row.get('Lane'):
            errors.append(f"Row {row_num}: Lane must be provided.")
        if not row.get('FASTQ1') and not row.get('BAM'):
            errors.append(f"Row {row_num}: Either FASTQ1 or BAM must be provided.")
        if row.get('RunType') == 'PAIRED_END' and not (row.get('FASTQ2') or row.get('BAM')):
            errors.append(f"Row {row_num}: FASTQ2 must be provided for PAIRED_END RunType.")
        if row.get('FASTQ1') and row.get('BAM'):
            errors.append(f"Row {row_num}: BAM should be ignored if FASTQ1 is provided.")
        if not row.get('MarkName'):
            errors.append(f"Row {row_num}: MarkName must be provided (either histone mark or input).")
        if not row.get('MarkType') or row.get('MarkType') not in ('B', 'N', 'I'):
            errors.append(f"Row {row_num}: MarkType must be provided and should be 'B' (for broad), 'N' (for narrow) , or 'I' (for input).")
    elif pipeline == 'longread_dnaseq':
        if not row.get('Run'):
            errors.append(f"Row {row_num}: Run must be provided.")
        if not row.get('BAM') and not row.get('FASTQ') and not row.get('FAST5'):
            errors.append(f"Row {row_num}: Either BAM, FASTQ or FAST5 must be provided for longread_dnaseq.")
        if (row.get('FASTQ') or row.get('FAST5')) and not row.get('Summary'):
            errors.append(f"Row {row_num}: Summary must be provided for the nanopore protocol.")
    elif pipeline == 'nanopore_covseq':
        if not row.get('Run'):
            errors.append(f"Row {row_num}: Run must be provided.")
        if not row.get('Summary'):
            errors.append(f"Row {row_num}: Summary must be provided.")
        if not row.get('FASTQ'):
            errors.append(f"Row {row_num}: FASTQ must be provided for nanopore_covseq.")
    elif pipeline == 'ampliconseq':
        if not row.get('FASTQ1'):
            errors.append(f"Row {row_num}: FASTQ1 must be provided for ampliconseq.")
        if row.get('RunType') == 'PAIRED_END' and not row.get('FASTQ2'):
            errors.append(f"Row {row_num}: FASTQ2 must be provided for PAIRED_END RunType.")

    return errors

def readset_run_type_check(string_array):
    '''
    Checks if Run type is "SINGLE_END" or "PAIRED_END"
    '''
    errors = []
    for row_num, row in enumerate(string_array, start=1):
        run_type = row.get('RunType')
        if run_type not in ["SINGLE_END", "PAIRED_END"]:
            errors.append(f"Row {row_num}: RunType must be 'SINGLE_END' or 'PAIRED_END'.")
    return errors

def readset_adapter_check(string_array):
    '''
    Checks if Adapter sequences are a combination of {A,C,T,G}
    '''
    errors = []
    valid_nucleotides = {'A', 'C', 'T', 'G'}
    for row_num, row in enumerate(string_array, start=1):
        adapters = [row.get('Adapter1', '').strip().upper(), row.get('Adapter2', '').strip().upper()]
        for adapter in adapters:
            if adapter and not all(nucleotide in valid_nucleotides for nucleotide in adapter):
                errors.append(f"Row {row_num}: {adapter} contains invalid nucleotides. Allowed: {valid_nucleotides}.")
    return errors

def readset_integer_check(string_array):
    '''
    Checks if QualityOffset is an integer
    '''
    errors = []
    for row_num, row in enumerate(string_array, start=1):
        quality_offset = row.get('QualityOffset')
        if quality_offset:
            try:
                int(quality_offset.strip())
            except ValueError:
                errors.append(f"Row {row_num}: QualityOffset must be an integer.")
    return errors

def readset_uniqueness_check(string_array):
    '''
    Checks if Readset values are unique and reports the row number and non-unique Readset value
    '''
    readset_dict = {}
    errors = []
    for row_num, row in enumerate(string_array, start=1):
        readset = row.get('Readset')
        if readset in readset_dict:
            errors.append(f"Row {row_num}: Readset '{readset}' is not unique. Previously found in row {readset_dict[readset]}.")
        else:
            readset_dict[readset] = row_num
    if errors:
        raise ValidationError("\n".join(errors))
    return True

def design_structure_check(string_array, pipeline):
    '''
    Checks if the entries in a design file are valid based on the pipeline
    '''
    errors = []
    if pipeline == 'chipseq':
        if 'MarkName' not in string_array[0]:
            errors.append("The design file must have a 'MarkName' column for the chipseq pipeline.")
    else:
        if 'Sample' not in string_array[0]:
            errors.append("The design file must have a 'Sample' column.")

    valid_entries = {0, 1, 2}
    for row_num, row in enumerate(string_array, start=1):
        if pipeline == 'chipseq' and not row.get('MarkName'):
            errors.append(f"Row {row_num}: MarkName must not be empty.")
        if not row.get('Sample'):
            errors.append(f"Row {row_num}: Sample must not be empty.")
        try:
            nums = [int(s.strip()) for s in row.values() if s.strip().isdigit()]
            if not all(n in valid_entries for n in nums):
                errors.append(f"Row {row_num}: {nums} contains values that are not permitted ({valid_entries}).")
        except ValueError:
            errors.append(f"Row {row_num}: Entries in design file matrix should be integers.")

    return errors

def output_pass(checks, file_type):
    '''
    Checks if all tests have passed
    '''
    if all(checks.values()):
        print(f'''
---------------------------------------------------------------------------
                Your {file_type} file has passed the check!
---------------------------------------------------------------------------
        ''')
    else:
        for key, value in checks.items():
            if not value:
                print(f"{key} has failed the test. Please fix this before launching the pipelines.")

def readset_validator(readset_file, pipeline):
    '''
    Validator for readset file. activated by -r
    '''
    readset_dict = {}
    readset_list_dict = load_file(readset_file)
    try:
        expected_header = readset_header_check(list(readset_list_dict[0].keys()), pipeline)

        readset_dict['unicodePass'] = unicode_check(readset_list_dict)
        readset_dict['DosNewlinePass'] = dos_newline_check(readset_list_dict)
        readset_dict['headerPass'] = True  # Already checked in readset_header_check
        readset_dict['rowSizePass'] = row_size_check(readset_list_dict)
        readset_dict['trailingSpacesPass'] = trailing_spaces_check(readset_list_dict)
        readset_dict['uniquenessPass'] = readset_uniqueness_check(readset_list_dict)

        all_errors = []
        if 'RunType' in expected_header:
            all_errors.extend(readset_run_type_check(readset_list_dict))
        if 'Adapter1' in expected_header and 'Adapter2' in expected_header:
            all_errors.extend(readset_adapter_check(readset_list_dict))
        if 'QualityOffset' in expected_header:
            all_errors.extend(readset_integer_check(readset_list_dict))

        for row_num, row in enumerate(readset_list_dict, start=1):
            errors = check_column_dependencies(row, pipeline, row_num)
            all_errors.extend(errors)

        if all_errors:
            raise ValidationError("\n".join(all_errors))
    except ValidationError as e:
        print(e)
        sys.exit(2)
    output_pass(readset_dict, "readset")

def design_validator(design_file, pipeline):
    '''
    Validator for design file. activated by -d
    '''
    design_dict = {}
    design_str = load_file(design_file)
    try:
        design_dict['unicodePass'] = unicode_check(design_str)
        design_dict['DosNewlinePass'] = dos_newline_check(design_str)
        design_dict['rowSizePass'] = row_size_check(design_str)
        design_dict['trailingSpacesPass'] = trailing_spaces_check(design_str)

        all_errors = design_structure_check(design_str, pipeline)

        if all_errors:
            raise ValidationError("\n".join(all_errors))
    except ValidationError as e:
        print(e)
        sys.exit(2)
    output_pass(design_dict, "design")

def main():
    script_name = os.path.basename(__file__)
    parser = argparse.ArgumentParser(
        description='Validate the basic structure and integrity of files used by the GenPipes pipelines.',
        epilog=f'''
        {script_name} validates the basic structure and integrity of files used
        by the GenPipes pipelines. User must provide a readset file or a design file.'''
    )
    parser.add_argument('-r', '--readset', type=argparse.FileType('r'), help='Readset file to validate')
    parser.add_argument('-d', '--design', type=argparse.FileType('r'), help='Design file to validate')
    parser.add_argument('-p', '--pipeline', type=str, required=True, help='Pipeline name to validate against', choices=['ampliconseq', 'chipseq', 'covseq', 'dnaseq', 'methylseq', 'longread_dnaseq', 'nanopore_covseq', 'rnaseq', 'rnaseq_denovo_assembly', 'rnaseq_light'])
    args = parser.parse_args()

    if not args.readset and not args.design:
        parser.error("You need to provide a readset file (-r) or a design file (-d) to be validated!")

    if args.readset:
        readset_validator(args.readset.name, args.pipeline)

    if args.design:
        design_validator(args.design.name, args.pipeline)

if __name__ == "__main__":
    main()
