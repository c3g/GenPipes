<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [BioDiscML Pipeline](#biodiscml-pipeline)
  - [Input files](#input-files)
  - [Configuration](#configuration)
  - [Output files](#output-files)
  - [Usage](#usage)
  - [biomarker_discovery](#biomarker_discovery)
  - [report](#report)
  - [References](#references)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


BioDiscML Pipeline
==================

Runs [BioDiscML](https://github.com/mickaelleclercq/BioDiscML) (Biomarker Discovery by Machine Learning)
on an arbitrary omics feature matrix (gene expression, proteomics, metabolomics, methylation,
clinico-pathological variables...) to discover a minimal biomarker signature and the best-performing
classifier for one or more contrasts.

BioDiscML automates the machine learning steps of biomarker discovery: data pre-processing, feature
ranking and selection, model selection and performance evaluation. It trains a large variety of
[Weka](https://www.cs.waikato.ac.nz/ml/weka/) classifiers across several feature search strategies
(forward/backward stepwise selection, "top k" features) and optimization criteria, evaluates each model
with several procedures (10-fold cross-validation, leave-one-out cross-validation, repeated holdout,
bootstrap, and an independent test set when sampling is enabled), and ranks them. It also retrieves
correlated features not included in the final model, to better understand the signature.

Unlike GenPipes' sequencing pipelines, this pipeline does not take FASTQ/BAM readsets as input: it takes a
features x samples matrix plus a GenPipes-style design file. It is deliberately standalone rather than a
step embedded in another pipeline, since BioDiscML itself is omics-agnostic: it can be pointed at the
output of rnaseq (`rawCountMatrix.tsv`, ideally after normalization and filtering), methylseq, or any
other feature table, GenPipes-produced or not.

Added as part of CFI Innovation Fund 2025 project #45950, Activity 1.4 ("New ML tools in GenPipes").

Input files
-----------

**Feature matrix** (`--feature-matrix`): tab-separated, features as rows and samples as columns. The first
column holds the feature IDs, the header row holds the sample names:

```text
gene	sample1	sample2	sample3	sample4
GENE_A	7.21	6.12	7.14	5.98
GENE_B	1.89	6.63	2.61	6.20
```

BioDiscML ranks features itself (information gain) and keeps at most the 1000 best, but its training time
grows quickly with the number of features and samples. Pre-filtering the matrix (e.g. removing lowly
expressed genes, or keeping the most variable/differential features) is strongly recommended for
large omics tables.

**Design file** (`--design`): the standard GenPipes
[design file](https://genpipes.readthedocs.io/en/latest/get-started/concepts/design_file.html): one row
per sample, one column per contrast, coded `1` (control), `2` (case) or `0` (excluded from that
contrast). Sample names must match the feature matrix header.

```text
sample	treated_vs_control
sample1	1
sample2	2
sample3	1
sample4	2
```

Each contrast is run independently, in its own `<contrast>/` output directory. Use `--contrasts` to run
only some of the design file's columns.

Configuration
-------------

BioDiscML is normally driven by a `config.conf` file. The pipeline generates it for each contrast from the
`[biodiscml]` section of the ini files:

| ini parameter | BioDiscML option | Description |
| --- | --- | --- |
| `class_column_name` | `classificationClassName` | Name of the class column built from the design file (`control`/`case`) |
| `sampling` | `sampling` | Split the data into a train set (2/3) and a test set (1/3) before training |
| `number_of_best_models` | `numberOfBestModels` | Number of models retained by `-bestmodel` |
| `number_of_best_models_sorting_metric` | `numberOfBestModelsSortingMetric` | Metric used to rank and retain the models (default `AVG_MCC`; any column of the results file, see [Output files](#output-files)) |
| `cpus` | `cpus` | Number of threads used by BioDiscML; keep equal to `[biomarker_discovery] cluster_cpu` |
| `other_options` | any | Extra `key=value` lines appended verbatim to `config.conf` (one per line) |

Any other BioDiscML option can be passed through `other_options`; see BioDiscML's
[config.conf template](https://github.com/mickaelleclercq/BioDiscML/blob/master/config.conf) for the
full list. For example, to restrict training to a few classifiers ("fast way") instead of BioDiscML's
full classifier collection:

```ini
[biodiscml]
other_options=
  classificationFastWay=true
  ccmd=bayes.NaiveBayes, AUC, top10
  ccmd=trees.J48, MCC, top10
```

`ram` sets the Java maximum heap (`-Xmx`); it must stay below `[biomarker_discovery] cluster_mem`, since
the Java 8 build used by BioDiscML does not detect the job's memory limit by itself.

BioDiscML's full training can take many hours on large matrices: adjust
`[biomarker_discovery] cluster_walltime`, `cluster_cpu` and `cluster_mem` to the size of the data.

Output files
------------

For each contrast, the pipeline writes into `<contrast>/`:

- `<contrast>.biodiscml_input.csv`: the training file given to BioDiscML (samples x features, plus the
  class column).
- `config.conf`: the generated BioDiscML configuration.
- `<contrast>_a.*`, `<contrast>_b.*`: BioDiscML's intermediate files (merged input data, then features
  ranked by information gain), in CSV and ARFF (Weka) formats.
- `<contrast>_c.classification.results.csv`: summary of all trained models with their evaluation metrics
  and selected features. Its main columns are:

  | Column | Description |
  | --- | --- |
  | `ID` | Model unique identifier |
  | `classifier`, `Options` | Weka classifier and its hyperparameters |
  | `OptimizedValue` | Criterion optimized during feature selection |
  | `SearchMode` | Feature search: forward (F), backward (B), forward-backward (FB), backward-forward (BF) stepwise selection, or "top k" features |
  | `nbrOfFeatures` | Number of features in the signature |
  | `TRAIN_10CV_*` | 10-fold cross-validation on the train set |
  | `TRAIN_LOOCV_*` | Leave-one-out cross-validation on the train set |
  | `TRAIN_RH_*`, `TRAIN_BS_*` | Repeated holdout and bootstrap on the train set |
  | `TEST_*` | Evaluation on the test set (empty when `sampling=false`) |
  | `TRAIN_TEST_RH_*`, `TRAIN_TEST_BS_*` | Repeated holdout and bootstrap on merged train and test sets |
  | `AVG_*`, `STD_*` | Average and standard deviation of the BER, MAE and MCC across all evaluations |
  | `AttributeList` | Selected features |

  Each evaluation reports accuracy (`ACC`), AUC, AUPRC, sensitivity (`SEN`), specificity (`SPE`),
  Matthews correlation coefficient (`MCC`), mean absolute error (`MAE`) and balanced error rate (`BER`),
  among others. See BioDiscML's
  [readme](https://github.com/mickaelleclercq/BioDiscML/blob/master/readme.md#output-files) for the
  complete list.
- For each retained model, `<contrast>_d.<model>.*`:
  - `*details.txt`: detailed information about the model, its performance, full signature and correlated
    features;
  - `*features.csv`: the values of the signature's features (`train_features.csv` for the train set,
    `all_features.csv` for train and test sets);
  - `*corrFeatures.csv`: the same, plus the features correlated with the signature;
  - `*.model`: the serialized Weka model, usable with BioDiscML's `-predict` mode on new data.
- `<contrast>.bestmodel.done`: marks the end of the `biomarker_discovery` step.
- `<contrast>.biodiscml_report.html`: the interactive HTML report (see [report](#report)).

Usage
-----

```text
usage: genpipes biodiscml [-h] [--clean] -c CONFIG [CONFIG ...]
                          [--container {wrapper, singularity} <IMAGE PATH>]
                          [-f] [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                          [--genpipes_file GENPIPES_FILE]
                          [-j {pbs,batch,daemon,slurm}] [--json-pt]
                          [-l {debug,info,warning,error,critical}]
                          [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                          [--wrap [WRAP]] --feature-matrix FEATURE_MATRIX
                          --design DESIGN_FILE [--contrasts CONTRASTS]

For more documentation, visit our website: https://genpipes.readthedocs.io

options:
  -h, --help            show this help message and exit
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -c, --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --genpipes_file, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -j, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  -l, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -o, --output-dir OUTPUT_DIR
                        output directory (default: current)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  -s, --steps STEPS     step range e.g. '1-5', '3,6,7', '2,4-8'
  --wrap [WRAP]         Path to the genpipes cvmfs wrapper script. Default is
                        genpipes/ressources/container/bin/container_wrapper.sh
                        . This is a convenience option for using genpipes in a
                        container
  --feature-matrix FEATURE_MATRIX
                        Feature matrix to mine for biomarkers: features x
                        samples, tab-separated, first column = feature ID,
                        header row = sample names (e.g. a GenPipes
                        rawCountMatrix.tsv, or any other omics feature table).
  --design DESIGN_FILE  GenPipes-style design file: one row per sample, one
                        column per contrast, coded 1 (control) / 2 (case) / 0
                        (excluded).
  --contrasts CONTRASTS
                        Comma-separated contrast column names (from the design
                        file) to run BioDiscML on. Default: all columns in the
                        design file.

Steps:

Protocol biodiscml
1 biomarker_discovery
2 report
```

Example:

```bash
genpipes biodiscml \
  --feature-matrix feature_matrix.tsv \
  --design design.tsv \
  -c $GENPIPES_INIS/biodiscml/biodiscml.base.ini $GENPIPES_INIS/common_ini/narval.ini \
  -o biodiscml_output \
  -g biodiscml_commands.sh
bash biodiscml_commands.sh
```

biomarker_discovery 
-------------------
 
For each requested contrast, reformats the feature matrix and design file into a BioDiscML
training file (samples x features, with a class column built from the design file's 1/2 coding),
generates a BioDiscML config.conf, trains models across BioDiscML's classifiers, then selects the
best model(s) based on the sorting metric configured in the [biodiscml] ini section.

report 
------
 
For each contrast, renders a self-contained interactive HTML report from the BioDiscML outputs
(`<contrast>.biodiscml_report.html`): ranked model leaderboard with the model(s) retained by
`-bestmodel` highlighted, and for each retained model a PCA with confidence ellipses, a
hierarchically clustered heatmap, the confusion matrix and per-feature boxplots. No external
dependency: Plotly is embedded, the file opens in any browser.

The report is rendered from a Jupyter notebook bundled with the pipeline
(`biodiscml_report.ipynb`); set `report_notebook` in the `[biodiscml]` ini section to use a custom one.

References
----------

Leclercq M, Vittrant B, Martin-Magniette ML, Scott Boyer MP, Perin O, Bergeron A, Fradet Y, Droit A.
Large-Scale Automatic Feature Selection for Biomarker Discovery in High-Dimensional OMICs Data.
*Front Genet.* 2019;10:452. [doi:10.3389/fgene.2019.00452](https://doi.org/10.3389/fgene.2019.00452)
([PMC6532608](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6532608/))
