################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import csv
import os

import shtab

# GenPipes Modules
from ...core.config import global_conf
from ...core.job import concat_jobs
from ...core.pipeline import Pipeline
from ...bfx import bash_cmd as bash
from ...bfx import biodiscml
from ..common import GenPipesPipeline

# Notebook bundled with the pipeline; override with `report_notebook` in the [biodiscml] ini section.
DEFAULT_REPORT_NOTEBOOK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "biodiscml_report.ipynb")


class BioDiscML(Pipeline):
    """
BioDiscML Pipeline
===================

Runs [BioDiscML](https://github.com/mickaelleclercq/BioDiscML) on an arbitrary omics feature matrix (gene
expression, proteomics, metabolomics, methylation, clinico-pathological variables...) to discover a minimal
biomarker signature and the best-performing classifier for one or more contrasts.

Unlike GenPipes' sequencing pipelines, BioDiscML does not take FASTQ/BAM readsets as input: it takes a
features x samples matrix plus a GenPipes-style design file (one row per sample, one column per contrast,
coded 1 = control / 2 = case / 0 = excluded). It is deliberately standalone rather than a step embedded in
another pipeline, since BioDiscML itself is omics-agnostic: it can be pointed at the output of rnaseq
(rawCountMatrix.tsv), methylseq, or any other feature table, GenPipes-produced or not.

Added as part of CFI Innovation Fund 2025 project #45950, Activity 1.4 ("New ML tools in GenPipes").
    """

    def __init__(self, *args, protocol="biodiscml", feature_matrix=None, design_file=None, contrasts=None, **kwargs):
        self.protocol = protocol
        self._feature_matrix = feature_matrix
        self._design_file = design_file
        self._contrasts_arg = contrasts
        super(BioDiscML, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument(
            "--feature-matrix",
            dest="feature_matrix",
            help="Feature matrix to mine for biomarkers: features x samples, tab-separated, first column "
                 "= feature ID, header row = sample names (e.g. a GenPipes rawCountMatrix.tsv, or any "
                 "other omics feature table).",
            type=argparse.FileType('r'),
            required=True
            ).complete = shtab.FILE
        cls._argparser.add_argument(
            "--design",
            dest="design_file",
            help="GenPipes-style design file: one row per sample, one column per contrast, coded "
                 "1 (control) / 2 (case) / 0 (excluded).",
            type=argparse.FileType('r'),
            required=True
            ).complete = shtab.FILE
        cls._argparser.add_argument(
            "--contrasts",
            dest="contrasts",
            help="Comma-separated contrast column names (from the design file) to run BioDiscML on. "
                 "Default: all columns in the design file."
            )
        return cls._argparser

    @property
    def feature_matrix(self):
        return self._feature_matrix

    @property
    def design_file(self):
        return self._design_file

    @property
    def contrasts(self):
        if self._contrasts_arg:
            return [contrast.strip() for contrast in self._contrasts_arg.split(",") if contrast.strip()]
        with open(self.design_file.name) as handle:
            header = next(csv.reader(handle, delimiter='\t'))
        return header[1:]

    def biomarker_discovery(self):
        """
        For each requested contrast, reformats the feature matrix and design file into a BioDiscML
        training file (samples x features, with a class column built from the design file's 1/2 coding),
        generates a BioDiscML config.conf, trains models across BioDiscML's classifiers, then selects the
        best model(s) based on the sorting metric configured in the [biodiscml] ini section.
        Best model selection runs single-threaded (`cpus=0`, in `config.bestmodel.conf`), to avoid a
        concurrency bug in BioDiscML 1.8.14 that can make it loop forever while retrieving correlated features.
        """
        feature_matrix = os.path.relpath(self.feature_matrix.name, self.output_dir)
        design_file = os.path.relpath(self.design_file.name, self.output_dir)

        jobs = []
        for contrast in self.contrasts:
            contrast_directory = os.path.join(self.output_dir, contrast)
            train_file = os.path.join(contrast_directory, f"{contrast}.biodiscml_input.csv")
            config_file = os.path.join(contrast_directory, "config.conf")
            # BioDiscML (v1.8.14) always inserts ".classification." for
            # doClassification=true runs, regardless of outputDir (see
            # bfx.biodiscml.config's docstring) -- confirmed by running it.
            results_file = os.path.join(contrast_directory, f"{contrast}_c.classification.results.csv")
            sentinel = os.path.join(contrast_directory, f"{contrast}.bestmodel.done")
            report_html = os.path.join(contrast_directory, f"{contrast}.biodiscml_report.html")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(contrast_directory),
                        biodiscml.format_input(
                            feature_matrix,
                            design_file,
                            contrast,
                            train_file
                        ),
                        biodiscml.config(
                            contrast,
                            train_file,
                            contrast_directory,
                            config_file
                        ),
                        biodiscml.train(
                            contrast_directory,
                            config_file,
                            train_file,
                            [results_file]
                        ),
                        biodiscml.best_model(
                            contrast_directory,
                            config_file,
                            results_file,
                            sentinel
                        )
                    ],
                    name=f"biomarker_discovery.{contrast}"
                )
            )

        return jobs

    def report(self):
        """
        For each contrast, renders a self-contained interactive HTML report from the BioDiscML outputs
        (`<contrast>.biodiscml_report.html`): ranked model leaderboard with the model(s) retained by
        `-bestmodel` highlighted, and for each retained model a PCA with confidence ellipses, a
        hierarchically clustered heatmap, the confusion matrix and per-feature boxplots. No external
        dependency: Plotly is embedded, the file opens in any browser.
        """
        notebook = global_conf.global_get('biodiscml', 'report_notebook', required=False) or DEFAULT_REPORT_NOTEBOOK

        jobs = []
        for contrast in self.contrasts:
            contrast_directory = os.path.join(self.output_dir, contrast)
            sentinel = os.path.join(contrast_directory, f"{contrast}.bestmodel.done")
            report_html = os.path.join(contrast_directory, f"{contrast}.biodiscml_report.html")

            jobs.append(
                concat_jobs(
                    [
                        biodiscml.report(
                            contrast_directory,
                            contrast,
                            notebook,
                            sentinel,
                            report_html
                        )
                    ],
                    name=f"report.{contrast}"
                )
            )

        return jobs

    # Not a GenPipesPipeline (no readsets), but its log_report step only relies on the base Pipeline
    log_report = GenPipesPipeline.log_report

    def protocols(self):
        return {
            "biodiscml": [
                self.biomarker_discovery,
                self.report,
                self.log_report
            ]
        }

    @property
    def step_list(self):
        return self.protocols()[self._protocol]


def main(parsed_args):
    """
    The function that will call this pipeline!
    Parameters:
        parsed_args (argparse.Namespace): The parsed arguments from the command line.
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    json_pt = parsed_args.json_pt
    force = parsed_args.force
    force_mem_per_cpu = parsed_args.force_mem_per_cpu
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps

    # Specific pipeline options
    feature_matrix = parsed_args.feature_matrix
    design_file = parsed_args.design_file
    contrasts = parsed_args.contrasts

    pipeline = BioDiscML(
        config_files,
        genpipes_file=genpipes_file,
        steps=steps,
        feature_matrix=feature_matrix,
        design_file=design_file,
        contrasts=contrasts,
        clean=clean,
        force=force,
        force_mem_per_cpu=force_mem_per_cpu,
        job_scheduler=job_scheduler,
        output_dir=output_dir,
        json_pt=json_pt,
        container=container
    )

    pipeline.submit_jobs()
