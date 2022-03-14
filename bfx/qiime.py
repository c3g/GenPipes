################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def catenate(
    input_fastq,
    input_name,
    catenate_fasta
    ):

    inputs = input_fastq
    outputs = [catenate_fasta]
    
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_catenate', 'module_python']
        ],

        command="""\
$QIIME_HOME/split_libraries_fastq.py \\
  -i {input_files} \\
  --sample_id {sample_name} \\
  -o {dir_output} \\
  -r 30 \\
  -p 0.01 \\
  -n {sequence_max_n} \\
  --phred_offset {phred_offset} \\
  --barcode_type 'not-barcoded'""".format(
        input_files=','.join(input_fastq),
        sample_name=','.join(input_name),
        sequence_max_n=global_config_parser.param('qiime_catenate', 'sequence_max_n'),
        phred_offset=global_config_parser.param('qiime_catenate', 'phred_offset'),
        dir_output="catenate/",
        ),
        removable_files=[catenate_fasta]
    )

def otu_closed_ref_picking(
    input_without_chimer,
    output_directory,
    otu_file
    ):

    inputs = [input_without_chimer]
    outputs = [otu_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_picking', 'module_python'],
            ['qiime_otu_picking', 'module_vsearch']
        ],

        command="""\
$QIIME_HOME/pick_otus.py \\
  -i {input_without_chimer} \\
  -m {method} \\
  -r {reference_seqs_fp} \\
  -s {similarity_treshold} \\
  --suppress_new_clusters \\
  --threads {threads_number} \\
  -o {output_directory}""".format(
        input_without_chimer=input_without_chimer,
        method='usearch61_ref',
        reference_seqs_fp=global_config_parser.param('qiime_otu_picking', 'reference_seqs_fp'),
        similarity_treshold=global_config_parser.param('qiime_otu_picking', 'similarity'),
        threads_number=global_config_parser.param('qiime_otu_picking', 'threads'),
        output_directory=output_directory
        ),
        removable_files=[otu_file]
    )

def otu_open_ref_picking(
    input_without_chimer,
    output_directory,
    otu_file
    ):

    inputs = [input_without_chimer]
    outputs = [otu_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_picking', 'module_python'],
            ['qiime_otu_picking', 'module_vsearch']
        ],

        command="""\
$QIIME_HOME/pick_otus.py \\
  -i {input_without_chimer} \\
  -m {method} \\
  -r {reference_seqs_fp} \\
  -s {similarity_treshold} \\
  --threads {threads_number} \\
  -o {output_directory}""".format(
        input_without_chimer=input_without_chimer,
        method='usearch61_ref',
        reference_seqs_fp=global_config_parser.param('qiime_otu_picking', 'reference_seqs_fp'),
        similarity_treshold=global_config_parser.param('qiime_otu_picking', 'similarity'),
        threads_number=global_config_parser.param('qiime_otu_picking', 'threads'),
        output_directory=output_directory
        ),
        removable_files=[otu_file]
    )

def otu_denovo_picking(
    input_without_chimer,
    output_directory,
    output_otus
    ):

    inputs = [input_without_chimer]
    outputs = [output_otus]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_picking', 'module_python'],
            ['qiime_otu_picking', 'module_vsearch']
        ],

        command="""\
$QIIME_HOME/pick_otus.py \\
  -i {input_without_chimer} \\
  -m {method} \\
  -s {similarity_treshold} \\
  --threads {threads_number} \\
  -o {output_directory}""".format(
        input_without_chimer=input_without_chimer,
        method='usearch61',
        similarity_treshold=global_config_parser.param('qiime_otu_picking', 'similarity'),
        threads_number=global_config_parser.param('qiime_otu_picking', 'threads'),
        output_directory=output_directory
        ),
        removable_files=[output_otus]
    )

def otu_rep_picking(
    otu_file,
    filter_fasta,
    otu_rep_file
    ):

    inputs = [otu_file, filter_fasta]
    outputs = [otu_rep_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_rep_picking', 'module_python']
        ],

        command="""\
$QIIME_HOME/pick_rep_set.py \\
  -i {otu_file} \\
  -f {filter_fasta} \\
  -m {method} \\
  -o {otu_rep_file}""".format(
        otu_file=otu_file,
        filter_fasta=filter_fasta,
        method=global_config_parser.param('qiime_rep_picking', 'rep_set_picking_method'),
        otu_rep_file=otu_rep_file
        ),
        removable_files=[otu_rep_file]
    )

def otu_assigning(
    otu_rep_picking_fasta,
    output_directory,
    tax_assign_file
    ):

    inputs = [otu_rep_picking_fasta]
    outputs = [tax_assign_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_assigning', 'module_python']
        ],

        command="""\
$QIIME_HOME/parallel_assign_taxonomy_uclust.py \\
  -i {otu_rep_picking_fasta} \\
  -T \\
  -O {threads_number} \\
  -r {database_otus} \\
  -t {taxonomy_otus} \\
  -o {output_directory}""".format(
        otu_rep_picking_fasta=otu_rep_picking_fasta,
        threads_number=global_config_parser.param('qiime_otu_assigning', 'threads'),
        database_otus=global_config_parser.param('qiime_otu_assigning', 'reference_seqs_fp'),
        taxonomy_otus=global_config_parser.param('qiime_otu_assigning', 'id_to_taxonomy_fp'),
        output_directory=output_directory
        ),
        removable_files=[tax_assign_file]
    )

def otu_table(
    otu_file,
    tax_assign_file,
    otu_table_file,
    otu_table_summary
    ):

    inputs = [otu_file, tax_assign_file]
    outputs = [otu_table_file,otu_table_summary]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_table', 'module_python']
        ],

        command="""\
$QIIME_HOME/make_otu_table.py \\
  -i {otu_file} \\
  -t {tax_assign_file} \\
  -o {otu_table_file}""".format(
        otu_file=otu_file,
        tax_assign_file=tax_assign_file,
        otu_table_file=otu_table_file
        ),
        removable_files=[otu_table_file,otu_table_summary]
    )

def otu_alignment(
    otu_rep_picking_fasta,
    output_directory,
    align_seq_fasta
    ):

    inputs = [otu_rep_picking_fasta]
    outputs = [align_seq_fasta]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_otu_assigning', 'module_python']
        ],

        command="""\
$QIIME_HOME/parallel_align_seqs_pynast.py \\
  -i {otu_rep_picking_fasta} \\
  -T \\
  --jobs_to_start {threads_number} \\
  -o {output_directory}""".format(
        otu_rep_picking_fasta=otu_rep_picking_fasta,
        threads_number=global_config_parser.param('qiime_otu_alignment', 'threads'),
        output_directory=output_directory
        ),
        removable_files=[align_seq_fasta]
    )

def filter_alignment(
    align_seq_fasta,
    output_directory,
    filter_align_fasta
    ):

    inputs = [align_seq_fasta]
    outputs = [filter_align_fasta]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_filter_alignment', 'module_python']
        ],

        command="""\
$QIIME_HOME/filter_alignment.py \\
  -i {align_seq_fasta} \\
  -o {output_directory}""".format(
        align_seq_fasta=align_seq_fasta,
        output_directory=output_directory
        ),
        removable_files=[filter_align_fasta]
    )

def phylogeny(
    filter_align_fasta,
    phylo_file
    ):

    inputs = [filter_align_fasta]
    outputs = [phylo_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_phylogeny', 'module_python']
        ],

        command="""\
$QIIME_HOME/make_phylogeny.py \\
  -i {filter_align_fasta} \\
  -o {phylo_file}""".format(
        filter_align_fasta=filter_align_fasta,
        phylo_file=phylo_file
        ),
        removable_files=[phylo_file]
    )

def multiple_rarefaction(
    otus_input,
    rarefied_otu_directory
    ):

    inputs = otus_input
    outputs = [rarefied_otu_directory]
    
    if not isinstance(inputs, list):
        inputs=[inputs]
    
    return Job(
        inputs,
        outputs,
        [
            ['qiime_multiple_rarefaction', 'module_python']
        ],

        command="""\
$QIIME_HOME/multiple_rarefactions.py \\
  -i {otus_input} \\
  -m {multiple_rarefaction_min} \\
  -x {multiple_rarefaction_max} \\
  -s {multiple_rarefaction_step} \\
  -n 3 \\
  --output_path {rarefied_otu_directory}""".format(
        otus_input=otus_input[0],
        multiple_rarefaction_min=global_config_parser.param('qiime_multiple_rarefaction', 'multiple_rarefaction_min'),
        multiple_rarefaction_max=global_config_parser.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max'),
        multiple_rarefaction_step=global_config_parser.param('qiime_multiple_rarefaction', 'multiple_rarefaction_step'),
        rarefied_otu_directory=rarefied_otu_directory
        ),
        removable_files=[rarefied_otu_directory]
    )

def alpha_diversity(
    rarefied_otu_directory,
    alpha_diversity_directory
    ):

    inputs = [rarefied_otu_directory]
    outputs = [alpha_diversity_directory]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_alpha_diversity', 'module_python']
        ],

        command="""\
$QIIME_HOME/alpha_diversity.py \\
  -i {rarefied_otu_directory} \\
  -m {metrics} \\
  -o {alpha_diversity_directory}""".format(
        rarefied_otu_directory=rarefied_otu_directory,
        metrics="observed_species,chao1,shannon",
        alpha_diversity_directory=alpha_diversity_directory
        ),
        removable_files=[alpha_diversity_directory]
    )

def collate_alpha(
    alpha_diversity_directory,
    alpha_diversity_collated_merge_directory,
    chao1_stat,
    observed_species_stat,
    shannon_stat
    ):

    inputs = [alpha_diversity_directory]
    outputs = [chao1_stat,observed_species_stat,shannon_stat]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_collate_alpha', 'module_python']
        ],

        command="""\
$QIIME_HOME/collate_alpha.py \\
  -i {alpha_diversity_directory} \\
  -o {alpha_diversity_collated_merge_directory}""".format(
        alpha_diversity_directory=alpha_diversity_directory,
        alpha_diversity_collated_merge_directory=alpha_diversity_collated_merge_directory
        ),
        removable_files=[chao1_stat,observed_species_stat,shannon_stat]
    )

def sample_rarefaction_plot(
    chao1_stat,
    observed_species_stat,
    shannon_stat,
    sample_collated_directory,    
    sample_map,
    sample_rarefaction_directory,
    curve_sample
    ):

    inputs = [chao1_stat,observed_species_stat,shannon_stat]
    outputs = [sample_rarefaction_directory] + curve_sample

    return Job(
        inputs,
        outputs,
        [
            ['qiime_sample_rarefaction', 'module_python']
        ],

        command="""\
MPLBACKEND=Agg $QIIME_HOME/make_rarefaction_plots.py \\
  -i {sample_collated_directory} \\
  -m {sample_map} \\
  -o {sample_rarefaction_directory}""".format(
        sample_collated_directory=sample_collated_directory,
        sample_map=sample_map,
        sample_rarefaction_directory=sample_rarefaction_directory
        ),
        removable_files=[sample_rarefaction_directory]
    )

def single_rarefaction(
    table_otu,
    chao1_rarefied_stat,
    observed_species_rarefied_stat,
    shannon_rarefied_stat,
    otu_normalized_table,
    normalization_method
    ):

    inputs = [table_otu]
    outputs = [chao1_rarefied_stat,observed_species_rarefied_stat,shannon_rarefied_stat,otu_normalized_table,normalization_method]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_single_rarefaction', 'module_python']
        ],

        command="""\
$QIIME_HOME/single_rarefaction.py \\
  -i {table_otu} \\
  -o {otu_normalized_table} \\
  -d {depth}""".format(
        table_otu=table_otu,
        otu_normalized_table=otu_normalized_table,
        depth=global_config_parser.param('qiime_single_rarefaction', 'single_rarefaction_depth')
        ),
        removable_files=[otu_normalized_table]
    )

def css_normalization(
    table_otu,
    chao1_rarefied_stat,
    observed_species_rarefied_stat,
    shannon_rarefied_stat,
    otu_normalized_table,
    normalization_method
    ):

    inputs = [table_otu]
    outputs = [chao1_rarefied_stat,observed_species_rarefied_stat,shannon_rarefied_stat,otu_normalized_table,normalization_method]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_css_normalization', 'module_python'],
            ['qiime_css_normalization', 'module_R']
        ],

        command="""\
$QIIME_HOME/normalize_table.py \\
  -i {table_otu} \\
  -o {otu_normalized_table} \\
  -a {method}""".format(
        table_otu=table_otu,
        otu_normalized_table=otu_normalized_table,
        method="CSS"
        ),
        removable_files=[otu_normalized_table]
    )

def rarefaction_plot(
    alpha_diversity_collated_merge_rarefied_directory,
    chao1_stat,
    observed_species_stat,
    shannon_stat,
    map_file,
    alpha_diversity_rarefaction_file,
    alpha_diversity_rarefaction_rarefied_directory
    ):

    inputs = [chao1_stat,observed_species_stat,shannon_stat]
    outputs = [alpha_diversity_rarefaction_file]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_rarefaction_plot', 'module_python']
        ],

        command="""\
MPLBACKEND=Agg $QIIME_HOME/make_rarefaction_plots.py \\
  -i {alpha_diversity_collated_merge_rarefied_directory} \\
  -m {map_file} \\
  -o {alpha_diversity_rarefaction_rarefied_directory}""".format(
        alpha_diversity_collated_merge_rarefied_directory=alpha_diversity_collated_merge_rarefied_directory,
        map_file=map_file,
        alpha_diversity_rarefaction_rarefied_directory=alpha_diversity_rarefaction_rarefied_directory
        ),
        removable_files=[alpha_diversity_rarefaction_file]
    )

def summarize_taxa(
    otu_normalized_table,
    taxonomic_directory,
    taxonomic_phylum,
    taxonomic_class,
    taxonomic_order,
    taxonomic_family,
    taxonomic_genus
    ):

    inputs = [otu_normalized_table]
    outputs = [taxonomic_directory, taxonomic_phylum, taxonomic_class, taxonomic_order, taxonomic_family, taxonomic_genus]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_summarize_taxa', 'module_python']
        ],

        command="""\
$QIIME_HOME/summarize_taxa.py \\
  -i {otu_normalized_table} \\
  -a \\
  -o {taxonomic_directory}""".format(
        otu_normalized_table=otu_normalized_table,
        taxonomic_directory=taxonomic_directory
        ),
        removable_files=[taxonomic_directory, taxonomic_phylum, taxonomic_class, taxonomic_order, taxonomic_family, taxonomic_genus]
    )

def plot_taxa(
    taxonomic_input,
    alpha_diversity_taxonomy_bar_plot,
    taxonomic_directory
    ):

    inputs = taxonomic_input
    outputs = [alpha_diversity_taxonomy_bar_plot]

    if not isinstance(inputs, list):
        inputs=[inputs]
    
    return Job(
        inputs,
        outputs,
        [
            ['qiime_plot_taxa', 'module_python']
        ],

        command="""\
MPLBACKEND=Agg $QIIME_HOME/plot_taxa_summary.py \\
  -i {taxonomic_input} \\
  -l {label} \\
  -t png \\
  -x 20 \\
  -y 12 \\
  -c {chart_type} \\
  -o {taxonomic_directory}""".format(
        taxonomic_input=",".join(taxonomic_input),
        label="Phylum,Class,Order,Family,Genus",
        chart_type="bar",
        taxonomic_directory=taxonomic_directory
        ),
        removable_files=[alpha_diversity_taxonomy_bar_plot]
    )

def beta_diversity(
    metric,
    otu_normalized_table,
    phylogenetic_tree_file,
    dm_directory,
    dm_unweighted_file,
    dm_weighted_file,
    dm_euclidean_file
    ):

    if metric == 'unifrac':

        inputs = [otu_normalized_table, phylogenetic_tree_file]
        outputs = [dm_unweighted_file,dm_weighted_file]

        return Job(
            inputs,
            outputs,
            [
                ['qiime_beta_diversity', 'module_python']
            ],

            command="""\
$QIIME_HOME/beta_diversity.py \\
  -i {otu_normalized_table} \\
  -t {phylogenetic_tree_file} \\
  -o {dm_directory}""".format(
            otu_normalized_table=otu_normalized_table,
            phylogenetic_tree_file=phylogenetic_tree_file,
            dm_directory=dm_directory
            ),
            removable_files=[dm_unweighted_file,dm_weighted_file]
        )

    else:

        inputs = [otu_normalized_table]
        outputs = [dm_euclidean_file]

        return Job(
            inputs,
            outputs,
            [
                ['qiime_beta_diversity', 'module_python']
            ],

            command="""\
$QIIME_HOME/beta_diversity.py \\
  -i {otu_normalized_table} \\
  -m {metrics} \\
  -o {dm_directory}""".format(
            otu_normalized_table=otu_normalized_table,
            metrics='euclidean',
            dm_directory=dm_directory
            ),
            removable_files=[dm_euclidean_file]
        )

def pcoa(
    metric,
    dm_unweighted_file,
    dm_weighted_file,
    dm_euclidean_file,
    dm_directory,
    pcoa_directory,
    pcoa_unweighted_file,
    pcoa_weighted_file,
    pcoa_euclidean_file
    ):

    if metric == 'unifrac':

        inputs = [dm_unweighted_file, dm_weighted_file]
        outputs = [pcoa_unweighted_file,pcoa_weighted_file]

        return Job(
            inputs,
            outputs,
            [
                ['qiime_pcoa', 'module_python']
            ],

            command="""\
$QIIME_HOME/principal_coordinates.py \\
  -i {dm_directory} \\
  -o {pcoa_directory}""".format(
            dm_directory=dm_directory,
            pcoa_directory=pcoa_directory
            ),
            removable_files=[pcoa_directory]
        )

    else:

        inputs = [dm_euclidean_file]
        outputs = [pcoa_euclidean_file]

        return Job(
            inputs,
            outputs,
            [
                ['qiime_pcoa', 'module_python']
            ],

            command="""\
$QIIME_HOME/principal_coordinates.py \\
  -i {dm_directory} \\
  -o {pcoa_directory}""".format(
            dm_directory=dm_directory,
            pcoa_directory=pcoa_directory
            ),
            removable_files=[pcoa_directory]
        )

def pcoa_plot(
    pcoa_file,
    map_file,
    beta_diversity_pcoa,
    pcoa_plot_directory
    ):

    inputs = [pcoa_file]
    outputs = [beta_diversity_pcoa]

    return Job(
        inputs,
        outputs,
        [
            ['qiime_pcoa_plot', 'module_python']
        ],

        command="""\
$QIIME_HOME/make_2d_plots.py \\
  -i {pcoa_file} \\
  -m {map_file} \\
  -o {pcoa_plot_directory}""".format(
        pcoa_file=pcoa_file,
        map_file=map_file,
        pcoa_plot_directory=pcoa_plot_directory
        ),
        removable_files=[beta_diversity_pcoa]
    )

