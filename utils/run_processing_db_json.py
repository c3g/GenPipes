#!/usr/bin/env python3

import os
import re
import sys
import json
import logging
import argparse
import inspect

from datetime import datetime

def arg_parse():
    """
    Argument parser.

    Returns argparse parsed object.
    """
    parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__),
        description="Aggregate database information in json format following successful run processing"
        )
    parser.add_argument(
        "runinfofile",
        metavar="RUNINFOFILE",
        help="Run info file supplied by Freezeman as json"
        )
    parser.add_argument(
        "json",
        metavar="JSON",
        nargs="+",
        help="Run validation json file(s)"
        )
    parser.add_argument(
        "-o", "--output_folder",
        metavar="PATH",
        help="Output folder"
        )
    parser.add_argument(
        "-p", "--platform",
        help="Operation platform (Default: abacus)",
        default="abacus",
        )
    parser.add_argument(
        "-v", "--verbose",
        help="Increase verbose level (Default: Warning)",
        action='count',
        default=0,
        )
    parser.add_argument(
        "-q", "--quiet",
        help="Decrease verbose level (Levels: Info, Warning)",
        action='count',
        default=0,
        )
    return parser.parse_args()

def read_jsons(json_list):
    """
    Convert the list of jsons files supplied via arguments into a dictionary
    of json parsed objects.

    Args:
        json_list   list
            Paths of json files to load in the dictionary. These should be
            run_validation json obtained from GenPipes' run_processing.

    Returns:
        Dictionary of the json, with content available to call with keys. 
    """
    read_dict = {}
    for json_file in json_list:
        with open(json_file, mode="r") as readfile:
            read_dict[json_file] = json.load(readfile)
    return read_dict, json_list

def write_json(dictionary, filename=None):
    """
    Write a dictionary into a json file.

    Args:
        dictionary  dict
            Dictionary to be written into a json file.
        filename    str
            Filename to write the json file. If not supplied, the json will be
            written to stdout.
    """
    if filename:
        os.makedirs(os.path.dirname(f"{filename}.json"), exist_ok=True)
        with open(f"{filename}.json", mode="w") as writefile:
            writefile.write(json.dumps(dictionary, indent=4))
    else:
        sys.stdout.write(json.dumps(dictionary, indent=4))

class DatabaseJson():
    """
    Class to generate the database json object from the run_validation jsons
    """

    def __init__(self, args, project):
        self._json_dict, self._json_names = read_jsons(args.json)
        with open(args.runinfofile, mode="r") as runinfofile:
            self._runinfo_dict = json.load(runinfofile)
        self._args = args
        self._project = project

    def input_files(self):
        """
        Return the input files used to generate the database json object.
        """
        return self._json_dict.keys()

    @property
    def operation_platform(self):
        """
        Return the operation platform.
        """
        platform = self._args.platform
        if not re.search(self._args.platform, os.uname()[1]):
            LOGGER.warning(inspect.cleandoc(
                f"Operation platform supplied is \"{platform}\", but running from a different hostname \"{os.uname()[1]}\". The specified platform will be used but it is recommended to run on the same platform."
                ))
        return platform

    @property
    def project_ext_src(self):
        """
        Return the project_ext_src src to be FREEZMAN.
        """
        return "FREEZMAN"

    @property
    def project_ext_id(self):
        """
        Return the project FMS ID.
        """
        return int(self._project[0])

    @property
    def project_name(self):
        """
        Return the project name.
        """
        if self._project[1].startswith("MOH-Q"):
            return "MOH-Q"
        return self._project[1]

    @property
    def run_ext_id(self):
        """
        Return the run FMS ID.
        """
        if self.is_same_value_between_jsons("run_obj_id"):
            return self._json_dict[sorted(self._json_dict)[0]]["run_obj_id"]

    # @property
    # def run_name(self):
    #     """
    #     Return the run name.
    #     """
    #     if self.is_same_value_between_jsons("run"):
    #         return self._runinfo_dict["run_name"]

    @property
    def run_name(self):
        """
        Return the run folder.
        """
        return self._json_names[0].split("/")[-3]
        # if self.is_same_value_between_jsons("run"):
        #     return "_".join([
        #         datetime.strftime(
        #             datetime.strptime(
        #                 self._runinfo_dict["run_start_date"], "%Y-%m-%d"
        #                 ),
        #             "%y%m%d"
        #             ),
        #         self._json_dict[sorted(self._json_dict)[0]]["run"],
        #         self._runinfo_dict["container_barcode"],
        #         f"-{self._json_dict[sorted(self._json_dict)[0]]["seqtype"]}"
        #         ])

    @property
    def run_instrument(self):
        """
        Return the run instrument.
        """
        if self.is_same_value_between_jsons("seqtype"):
            return self._json_dict[sorted(self._json_dict)[0]]["seqtype"]

    @property
    def run_date(self):
        """
        Return the run date.
        """
        return f"{datetime.strptime(self._runinfo_dict["run_start_date"], '%Y-%m-%d')}"

    @property
    def specimen(self):
        """
        Return the specimen object.
        """
        list_of_dict = []

        # Helper function to find existing specimen
        def find_specimen(specimen_name):
            for specimen in list_of_dict:
                if specimen["specimen_name"] == specimen_name:
                    return specimen
            return None

        # Helper function to find existing sample by name
        def find_sample(specimen, sample_name):
            for sample in specimen["sample"]:
                if sample["sample_name"] == sample_name:
                    return sample
            return None

        for json_path, section in self._json_dict.items():
            for rdst in section["readsets"]:
                sample_name = section["readsets"][rdst]["sample_name"]
                # TODO
                # Extract specimen name, cohort and institution from the sample name only for MoHQ samples.
                # To be adapted for other projects
                result = re.search(r"^((MoHQ-(JG|CM|GC|MU|MR|XX|HM)-\w+)-\w+)-\w+-\w+(D|R)(T|N)", sample_name)
                if result:
                    specimen_name = result.group(1)
                    cohort = result.group(2)
                    institution = result.group(3)
                else:
                    continue  # Skip if the specimen name pattern does not match

                # Check if the specimen already exists
                specimen = find_specimen(specimen_name)
                if not specimen:
                    specimen = {
                        "specimen_ext_id": None,
                        "specimen_ext_src": None,
                        "specimen_name": specimen_name,
                        "specimen_cohort": cohort,
                        "specimen_institution": institution,
                        "sample": []
                    }
                    list_of_dict.append(specimen)

                # Check if the sample already exists within the specimen by name
                sample = find_sample(specimen, sample_name)
                if not sample:
                    sample = {
                        "sample_ext_id": int(section["readsets"][rdst]["derived_sample_obj_id"]),
                        "sample_ext_src": self.project_ext_src,
                        "sample_name": sample_name,
                        "sample_tumour": sample_name.endswith("T"),
                        "readset": []
                    }
                    specimen["sample"].append(sample)

                runinfo_sample = self.runinfo_sample(json_path, rdst)
                readset = {
                    "experiment_sequencing_technology": None,
                    "experiment_type": section["readsets"][rdst]["library_type"],
                    "experiment_nucleic_acid_type": "RNA" if section["readsets"][rdst]["library_type"] == "RNASeq" else "DNA",
                    "experiment_library_kit": runinfo_sample["library_kit"],
                    "experiment_kit_expiration_date": None,
                    "readset_name": f"{rdst}.{section['run']}_{section['lane']}",
                    "readset_lane": section["lane"],
                    "readset_index_name": section["readsets"][rdst]["barcodes"][0]["INDEX_NAME"],
                    "readset_index1": section["readsets"][rdst]["barcodes"][0]["INDEX1"],
                    "readset_index2": section["readsets"][rdst]["barcodes"][0]["INDEX2"],
                    "readset_adapter1": section["readsets"][rdst]["barcodes"][0]["ADAPTERi7"],
                    "readset_adapter2": section["readsets"][rdst]["barcodes"][0]["ADAPTERi5"],
                    "readset_barcode": section["readsets"][rdst]["barcodes"][0]["BARCODE_SEQUENCE"],
                    "readset_sequencing_type": section["sequencing_method"],
                    "file": self.return_files(json_path, rdst),
                    "metric": self.return_metrics(json_path, rdst),
                }
                sample["readset"].append(readset)

        return list_of_dict



    def properties(self):
        """
        Return the properties of the DatabaseJson object.
        """
        class_items = self.__class__.__dict__.items()
        return {k: getattr(self, k) for k, v in class_items if isinstance(v, property)}


    def is_same_value_between_jsons(self, run_validation_field):
        """
        Check if the value of the field is the same in all JSONs.
        """
        values = set()

        for section in self._json_dict.values():
            value = section.get(run_validation_field)
            if value is not None:
                values.add(value)

        if len(values) == 1:
            return True
        else:
            LOGGER.error(f"Multiple {run_validation_field} detected")
            raise Exception(
                f"Conflict in input values for \"{run_validation_field}\". Confirm that the input files are from a single run"
            )


    def return_files(self, json_path, rdst):
        """
        Return the files for the readset.
        """
        filetypes = ["fastq_1", "fastq_2", "bam", "bai"]
        files = []

        readsets = self._json_dict.get(json_path, {}).get("readsets", {})
        readset = readsets.get(rdst, {})

        for ftype in filetypes:
            file_info = readset.get(ftype, {})
            final_path = file_info.get("final_path")

            if final_path:
                files.append({
                    "location_uri": f"abacus://{final_path}",
                    "file_name": os.path.basename(final_path),
                    "file_deliverable": False
                })

        return files

    def return_metrics(self, json_path, rdst):
        """
        Return the metrics for the readset.
        """
        metric_matches = {
            "raw_reads_count": ["qc", "nb_reads"],
            "raw_duplication_rate": ["qc", "duplicate_rate"],
            "raw_median_insert_size": ["alignment", "median_aligned_insert_size"],
            "raw_mean_insert_size": ["alignment", "average_aligned_insert_size"],
            "raw_mean_coverage": ["alignment", "mean_coverage"]
        }
        metrics = []
        run_validations = self._json_dict.get(json_path, {}).get("run_validation", [])

        for db_key, run_v_key in metric_matches.items():
            for run_v in run_validations:
                if run_v.get("sample") == rdst:
                    if run_v_key[0] in run_v and run_v_key[1] in run_v[run_v_key[0]]:
                        metric_value = run_v[run_v_key[0]][run_v_key[1]]
                        if metric_value is not None:
                            metrics.append({"metric_name": db_key, "metric_value": metric_value})
                        elif run_v_key[0] is None:
                            metrics.append({"metric_name": db_key, "metric_value": None, "metric_flag": "MISSING"})
                        elif metric_value is None and self.runinfo_sample(json_path, rdst).get("reference_genome") is None:
                            metrics.append({"metric_name": db_key, "metric_value": None, "metric_flag": "NOT_APPLICABLE"})

        return metrics

    def runinfo_sample(self, json_path, sample_name):
        """
        Return the sample object from the runinfofile json.
        """
        samples = self._runinfo_dict["samples"]
        json_samples = self._json_dict[json_path]["readsets"]

        for sample in samples:
            if sample["sample_name"] == json_samples[sample_name]["sample_name"]:
                return sample
        return None  # In case no match is found


def main():
    """
    Main function to be run when script is called.
    """
    args = arg_parse()
    global LOGGER
    logging.basicConfig(level=30-args.verbose*10+args.quiet*10)
    LOGGER = logging.getLogger(__name__)
    json_dict, _ = read_jsons(args.json)
    projects = {}
    # Separate projects into different DatabaseJson objects
    for _, json_content in json_dict.items():
        for _, content in json_content["readsets"].items():
            if content["project_obj_id"] not in projects or content["project_name"] not in projects.values():
                projects[content["project_obj_id"]] = content["project_name"]
                LOGGER.info(
                    f"Project ID:{content["project_obj_id"]} \"{content["project_name"]}\" identified."
                    )
    databases = {}
    for project in projects.items():
        databases[project] = DatabaseJson(args, project)
        if args.output_folder:
            # TODO
            # Special case for MOH-Q, to be adapted for other projects
            if project[1].startswith("MOH-Q"):
                project_name = "MOH-Q"
            else:
                project_name = project[1]
            run_name = databases[project].run_name
            write_json(
                databases[project].properties(),
                filename=os.path.join(args.output_folder, f"{project_name}_{run_name}")
                )
        else:
            write_json(databases[project].properties())


if __name__ == "__main__":
    main()
