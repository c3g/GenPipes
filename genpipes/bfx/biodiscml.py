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

"""
Wrappers around BioDiscML (https://github.com/mickaelleclercq/BioDiscML), a
machine-learning-based biomarker discovery tool developed at Universite
Laval. Added to GenPipes as part of CFI Innovation Fund 2025 project #45950
(C3G core facility renewal), Activity 1.4 "New ML tools in GenPipes".

BioDiscML expects a single samples x features CSV file with a class column,
which is not the shape GenPipes count matrices come in (features x samples,
tab-separated, no class information). ``format_input`` bridges the two by
reformatting a GenPipes count matrix and design file into a BioDiscML-ready
training file for a given contrast.
"""

# Python Standard Modules
import os

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job


def format_input(
    count_matrix,
    design_file,
    contrast_name,
    output_csv
    ):
    """
    Reformat a GenPipes count matrix (features x samples, tab-separated)
    and its design file into the samples x features CSV expected by
    BioDiscML, appending a class column built from the design file's 1/2
    coding for the given contrast (1 = control, 2 = case). Samples coded 0,
    or absent from the contrast, are excluded.
    """
    class_column = global_conf.global_get('biodiscml', 'class_column_name')

    return Job(
        [count_matrix, design_file],
        [output_csv],
        [
            ['biodiscml', 'module_python']
        ],
        command="""\
( python3 - <<'EOF'
import csv

class_column = "{class_column}"
contrast_name = "{contrast_name}"

with open("{design_file}") as handle:
    reader = csv.reader(handle, delimiter='\\t')
    header = next(reader)
    contrast_index = header.index(contrast_name)
    labels = {{}}
    for row in reader:
        if not row:
            continue
        code = row[contrast_index]
        if code == "1":
            labels[row[0]] = "control"
        elif code == "2":
            labels[row[0]] = "case"

with open("{count_matrix}") as handle:
    reader = csv.reader(handle, delimiter='\\t')
    header = next(reader)
    sample_names = header[1:]
    kept_indices = [i for i, name in enumerate(sample_names) if name in labels]
    feature_rows = [row for row in reader if row]

with open("{output_csv}", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["sample"] + [row[0] for row in feature_rows] + [class_column])
    for i in kept_indices:
        sample_name = sample_names[i]
        writer.writerow([sample_name] + [row[i + 1] for row in feature_rows] + [labels[sample_name]])
EOF
)""".format(
            class_column=class_column,
            contrast_name=contrast_name,
            design_file=design_file,
            count_matrix=count_matrix,
            output_csv=output_csv
        )
    )


def config(
    project_name,
    train_file,
    work_dir,
    output_config
    ):
    """
    Generate a BioDiscML config.conf file for a training run on
    ``train_file`` (as produced by ``format_input``).

    Note: BioDiscML (checked against v1.8.14) does not actually honour the
    ``outputDir`` config key for its main output files (``*_a.*``, ``*_b.*``,
    ``*_c.*``, ``*_d.*`` -- confirmed by running it for real): it writes
    them to the process' current working directory regardless. So
    ``outputDir`` is deliberately left out of the generated config here;
    ``work_dir`` is instead used by ``train``/``best_model`` to ``cd`` into
    before invoking the jar, which is what actually controls where the
    outputs land.
    """
    class_column = global_conf.global_get('biodiscml', 'class_column_name')
    sampling = global_conf.global_get('biodiscml', 'sampling')
    number_of_best_models = global_conf.global_get('biodiscml', 'number_of_best_models')
    sorting_metric = global_conf.global_get('biodiscml', 'number_of_best_models_sorting_metric')
    cpus = global_conf.global_get('biodiscml', 'cpus')
    other_options = global_conf.global_get('biodiscml', 'other_options', required=False)

    return Job(
        [train_file],
        [output_config],
        command="""\
mkdir -p {work_dir} && \\
( cat > {output_config} << CONF
project={project_name}
trainFile={train_file}
excluded=sample
sampling={sampling}
doClassification=true
classificationClassName={class_column}
numberOfBestModels={number_of_best_models}
numberOfBestModelsSortingMetric={sorting_metric}
cpus={cpus}
{other_options}
CONF
)""".format(
            work_dir=work_dir,
            output_config=output_config,
            project_name=project_name,
            train_file=train_file,
            sampling=sampling,
            class_column=class_column,
            number_of_best_models=number_of_best_models,
            sorting_metric=sorting_metric,
            cpus=cpus,
            other_options=other_options
        )
    )


def train(
    work_dir,
    config_file,
    train_file,
    outputs
    ):
    """
    Train BioDiscML models: feature ranking, then model/feature-subset
    search and evaluation across the classifiers configured in
    ``config_file``. Produces the ``<project>_c.classification.results.csv``
    summary table used by ``best_model`` to extract the final signature(s).

    Runs with ``work_dir`` as the current directory, since that -- not
    ``outputDir`` in the config -- is what determines where BioDiscML
    writes its output files (see ``config``'s docstring).
    """
    return Job(
        [config_file, train_file],
        outputs,
        [
            ['biodiscml', 'module_java'],
            ['biodiscml', 'module_biodiscml']
        ],
        command="""\
cd {work_dir} && \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $BIODISCML_JAR \\
  -config {config_file} \\
  -train""".format(
            work_dir=work_dir,
            tmp_dir=global_conf.global_get('biodiscml', 'tmp_dir'),
            java_other_options=global_conf.global_get('biodiscml', 'java_other_options', required=False),
            ram=global_conf.global_get('biodiscml', 'ram'),
            config_file=config_file
        )
    )


def best_model(
    work_dir,
    config_file,
    results_file,
    sentinel
    ):
    """
    Select the best model(s) from a completed BioDiscML training run,
    according to the sorting metric configured in ``config_file``, and
    generate the corresponding signature/feature reports.

    Runs with ``work_dir`` as the current directory (see ``train``). The
    per-model detail/feature/model files BioDiscML writes are named after
    the winning classifier and its hyperparameters, so their exact names
    cannot be predicted ahead of the run; ``sentinel`` is an empty marker
    file touched on completion so GenPipes has a stable dependency to
    track. Validating this behaviour against real runs, and tightening the
    declared outputs accordingly, is part of the Activity 1.4 workplan.
    """
    return Job(
        [config_file, results_file],
        [sentinel],
        [
            ['biodiscml', 'module_java'],
            ['biodiscml', 'module_biodiscml']
        ],
        command="""\
cd {work_dir} && \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $BIODISCML_JAR \\
  -config {config_file} \\
  -bestmodel && \\
touch {sentinel}""".format(
            work_dir=work_dir,
            tmp_dir=global_conf.global_get('biodiscml', 'tmp_dir'),
            java_other_options=global_conf.global_get('biodiscml', 'java_other_options', required=False),
            ram=global_conf.global_get('biodiscml', 'ram'),
            config_file=config_file,
            sentinel=sentinel
        )
    )


def report(
    work_dir,
    project_name,
    notebook,
    sentinel,
    output_html
    ):
    """
    Render the self-contained interactive HTML report for a completed
    BioDiscML run (model leaderboard, per-model PCA / clustered heatmap /
    confusion matrix / boxplots). Executes the bundled Jupyter notebook
    ``notebook`` with nbconvert, with code cells hidden. The output has no
    external dependency (Plotly embedded, no CDN).

    ``sentinel`` is the marker file touched by ``best_model``: the report
    reads the results table and the per-model ``*_d.*features.csv`` files
    written by ``-bestmodel``. The sorting metric is the same one used for
    ``-bestmodel``, so the retained model(s) are highlighted consistently.
    """
    sorting_metric = global_conf.global_get('biodiscml', 'number_of_best_models_sorting_metric')

    return Job(
        [sentinel],
        [output_html],
        [
            # The report module declares its Python as a prereq, so load it first.
            ['biodiscml', 'module_python_report_base'],
            ['biodiscml', 'module_python_report']
        ],
        command="""\
BIODISCML_OUTPUT_DIR=$(cd {work_dir} && pwd) \\
PROJECT={project_name} \\
SORT_METRIC={sorting_metric} \\
jupyter nbconvert --to html --execute --no-input \\
  --HTMLExporter.mathjax_url="" \\
  --HTMLExporter.require_js_url="" \\
  --output $(cd {work_dir} && pwd)/{output_name} \\
  {notebook}""".format(
            work_dir=work_dir,
            project_name=project_name,
            sorting_metric=sorting_metric,
            output_name=os.path.basename(output_html),
            notebook=notebook
        )
    )
