# Copyright (C) 2026 Arnaud Belcour - Univ. Grenoble Alpes, Inria, Grenoble, France Microcosme
# Univ. Grenoble Alpes, Inria, Microcosme
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import os
import pandas as pd
import numpy as np

import json
import logging

import plotly.graph_objects as go

logger = logging.getLogger(__name__)

def diff_between_group(abundance_dataframe, positive_control_groups, negative_control_groups, group_abundance, threshold, output_file):
    """ Subselect organisms according to their presence/absence in subset of samples and then returns the sum of their abundance from a type of sample.

    Args:
        abundance_dataframe (pd.DataFrame): dataframe containing sample/ttime step abundances and taxon name
        positive_control_groups (list): list of groups that should contain organisms with positive abundance
        negative_control_groups (list):  list of groups that should contain organisms with zero abundance
        group_abundance (str): column to use to get the total abundace of organisms
        threshold (float): abundance threshold for positive control
        output_file (str): path to output file

    Returns:
        abundance_org (float): summed abundance of kept organisms
    """
    tmp_df = abundance_dataframe.copy()
    for group in positive_control_groups:
        tmp_df = tmp_df[tmp_df[group]>threshold]
    for group in negative_control_groups:
        tmp_df = tmp_df[tmp_df[group]==0]
    tmp_df.to_csv(output_file, sep='\t')
    abundance_org = tmp_df[group_abundance].sum()

    return abundance_org


def add_taxon_labels(selected_col, taxon_col, tmp_df, source_node, node_values, new_links, abundance_threshold=0):
    """ Add taxon name as nodes in Sankey diagram.

    Args:
        selected_col (str): column that will be used to compute the abundance of the organisms
        taxon_col (str): column containing taxon names
        tmp_df (pd.DataFrame): filtered dataframe containing abundance and taxon names
        source_node (str): source node number from which the edge will start
        node_values (dict): total value associated with a node
        new_links (dict): link between nodes with values
        abundance_threshold (int): abundance htreshold to show taxon name
    """
    for org in tmp_df[taxon_col].unique():
        org_df = tmp_df[tmp_df[taxon_col]==org]
        org_abund = org_df[selected_col].sum()
        if org_abund > abundance_threshold:
            if (source_node, org) not in new_links:
                new_links[(source_node, org)] = org_abund
            else:
                new_links[(source_node, org)] += org_abund
            if org not in node_values:
                node_values[org] = org_abund
            else:
                node_values[org] += org_abund


def generate_node_edges(all_edges, node_values):
    """ From dictionary of edges and nodes, filter edges and nodes to produce sankey inputs.
    Keep only nodes and edges with positive values.

    Args:
        all_edges (dict): dictionary containing edges as key (tuple with time step/organism name) and the associated abundance as value
        node_values (dict): dictionary containing node name as key and their abundance as values

    Returns:
        sources (list): list of node numbers that are source for the edge in Sankey diagram
        source_label_nodes (list): list of node labels matching sources list
        targets (list): list of node numbers that are target for the edge in Sankey diagram
        target_label_nodes (list): list of node labels matching targets list
        labels (list): labels of filtered node
        values (list): abundance values for the edges
    """
    sources = []
    source_label_nodes = []
    targets = []
    target_label_nodes = []
    labels = []
    values = []
    node_numbers = {}
    for link in all_edges:
        if all_edges[link] is not None and all_edges[link] > 0:
            source = link[0]
            target = link[1]
            if node_values[source] is not None and node_values[source] > 0 and node_values[target] is not None and node_values[target] > 0:
                if source not in node_numbers:
                    node_numbers[source] = len(node_numbers)
                    labels.append(source)
                if target not in node_numbers:
                    node_numbers[target] = len(node_numbers)
                    labels.append(target)
                sources.append(node_numbers[source])
                source_label_nodes.append(source)
                targets.append(node_numbers[target])
                target_label_nodes.append(target)
                values.append(all_edges[link])

    return sources, source_label_nodes, targets, target_label_nodes, labels, values


def generate_sankey_diagram(input_file, type_cols, taxon_col, output_folder, abundance_threshold=0, relative=None):
    """ Generate Sankey diagram from input file and input dictionary.

    Args:
        input_file (str): path to input abundance file
        type_cols (dict): dictionary associating samples and time steps
        taxon_col (str): column containing taxon names
        output_folder (str): path to output folder
        abundance_threshold (int): abundance htreshold to show taxon name
        relative (bool): if True, compute relative abundance instead of absolute
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    intermediary_folder = os.path.join(output_folder, 'intermediary_folder')
    if not os.path.exists(intermediary_folder):
        os.mkdir(intermediary_folder)

    df = pd.read_csv(input_file, sep='\t', index_col=0)
    df[taxon_col] = df[taxon_col].replace(np.nan, 'Unknown')
    samples_cols = [sample for res_type in type_cols for sample in type_cols[res_type]]
    cavern_name_to_ids = samples_cols + [taxon_col]
    df = df[cavern_name_to_ids]

    if relative is not None:
        for col in df.columns:
            if col in samples_cols:
                df[col] = df[col] / df[col].sum()

    threshold = 0
    for type_col in type_cols:
        df[type_col] = df[type_cols[type_col]].mean(axis=1)

    output_file = os.path.join(intermediary_folder, 'edge_kit_control_to_enrichment_T01.tsv')
    abund_org_unique_kc_T1 = diff_between_group(df, ['kit_control', 'enrichment_T01'], ['initial'], 'enrichment_T01', threshold, output_file)

    output_file = os.path.join(intermediary_folder, 'edge_initial_to_enrichment_T01.tsv')
    abund_org_unique_initial_T1 = diff_between_group(df, ['initial', 'enrichment_T01'], ['kit_control'], 'enrichment_T01', threshold, output_file)

    output_file = os.path.join(intermediary_folder, 'edge_kc_initial_to_enrichment_T01.tsv')
    abund_org_shared_kc_initial_T1 = diff_between_group(df, ['initial', 'kit_control', 'enrichment_T01'], [], 'enrichment_T01', threshold, output_file)

    negative_controls = ['kit_control', 'initial']
    if 'kit_control_T01' in type_cols:
        negative_controls.append('kit_control_T01')
    output_file = os.path.join(intermediary_folder, 'edge_appearing_to_enrichment_T01.tsv')
    abund_org_appearing_T1 = diff_between_group(df, ['enrichment_T01'], negative_controls, 'enrichment_T01', threshold, output_file)

    if 'kit_control_T01' in type_cols:
        output_file = os.path.join(intermediary_folder, 'edge_kit_control_T01_to_enrichment_T01.tsv')
        abund_org_kc_T1 = diff_between_group(df, ['kit_control_T01', 'enrichment_T01'], ['kit_control', 'initial'], 'enrichment_T01', threshold, output_file)
    else:
        abund_org_kc_T1 = None

    output_file = os.path.join(intermediary_folder, 'edge_kit_control_to_initial_dead.tsv')
    abund_org_dead_kc = diff_between_group(df, ['kit_control'], ['enrichment_T01'], 'kit_control', threshold, output_file)

    output_file = os.path.join(intermediary_folder, 'edge_initial_to_initial_dead.tsv')
    abund_org_dead_initial = diff_between_group(df, ['initial'], ['enrichment_T01'], 'initial', threshold, output_file)

    sankey_nodes = ['kit_control', 'initial', 'kc_initial', 'kit_control_T01', 'enrichment_T01', 'enrichment_T02']
    sankey_links = {('kit_control', 'enrichment_T01'): abund_org_unique_kc_T1, ('kc_initial', 'enrichment_T01'): abund_org_shared_kc_initial_T1,
                    ('initial', 'enrichment_T01'): abund_org_unique_initial_T1,
                    ('appearing', 'enrichment_T01'): abund_org_appearing_T1, ('kit_control_T01', 'enrichment_T01'): abund_org_kc_T1,
                   ('kit_control', 'initial_dead'): abund_org_dead_kc, ('initial', 'initial_dead'): abund_org_dead_initial}

    enrichment_times = sorted([type_col for type_col in type_cols if 'enrichment_' in type_col and type_col != 'enrichment_T01'])

    previous_time =  'enrichment_T01'
    for enrichment_time in enrichment_times:
        current_time = enrichment_time
        culture_dead = 'culture_dead_' + current_time
        kit_control = current_time.replace('enrichment_', 'kit_control_')
        appearing = 'appearing_' + current_time
        if enrichment_time in type_cols:
            if kit_control in type_cols:
                output_file = os.path.join(intermediary_folder, f'edge_{kit_control}_to_{current_time}.tsv')
                abund_org_kc_T1 = diff_between_group(df, [kit_control, current_time], ['kit_control', 'initial'], current_time, threshold, output_file)
            else:
                abund_org_kc_T1 = None
            output_file = os.path.join(intermediary_folder, f'edge_{previous_time}_to_{current_time}.tsv')
            abund_org_T1_T2 = diff_between_group(df, [previous_time, current_time], [], previous_time, threshold, output_file)
            output_file = os.path.join(intermediary_folder, 'abund_org_appearing_'+current_time+'.tsv')
            abund_org_appearing_T2 = diff_between_group(df, [current_time], [previous_time], current_time, threshold, output_file)
            output_file = os.path.join(intermediary_folder, f'edge_{previous_time}_to_{culture_dead}.tsv')
            abund_org_dead_T1 = diff_between_group(df, [previous_time], [current_time], previous_time, threshold, output_file)
            negative_controls = ['kit_control', 'initial', previous_time]
            output_file = os.path.join(intermediary_folder, f'edge_{appearing}_to_{current_time}.tsv')
            abund_org_appearing_T1 = diff_between_group(df, [current_time], negative_controls, current_time, threshold, output_file)
            sankey_links[(previous_time, current_time)] = abund_org_T1_T2
            sankey_links[(kit_control, current_time)] = abund_org_kc_T1
            sankey_links[(previous_time, culture_dead)] = abund_org_dead_T1
            sankey_links[(appearing, current_time)] = abund_org_appearing_T1
            previous_time = current_time

    node_values = {}
    already_processed_tuples = []
    all_nodes = [node for node_tuples in sankey_links for node in node_tuples]
    for sankey_node in all_nodes:
        if sankey_node in type_cols:
            node_links = [sankey_link for sankey_link in sankey_links if sankey_node in sankey_link]
            for node_link in node_links:
                if sankey_links[node_link] is not None and node_link not in already_processed_tuples:
                    already_processed_tuples.append(node_link)
                    if node_link[0] not in node_values: 
                        node_values[node_link[0]] = sankey_links[node_link]
                    else:
                        node_values[node_link[0]] += sankey_links[node_link]
                    if node_link[1] not in node_values: 
                        node_values[node_link[1]] = sankey_links[node_link]
                    else:
                        node_values[node_link[1]] += sankey_links[node_link]

    if len(enrichment_times) > 0:
        final_time = enrichment_times[-1]
    else:
        final_time = 'enrichment_T01'

    if final_time != 'enrichment_T01':
        tmp_df = df[df[final_time]>threshold]
        new_links= {}
        add_taxon_labels(final_time, taxon_col, tmp_df, final_time, node_values, new_links, abundance_threshold)
        all_culture_dead_tuples = [node_tuples for node_tuples in sankey_links for node in node_tuples if 'culture_dead_' in node]
        for all_culture_dead_tuple in all_culture_dead_tuples:
            previous_time = all_culture_dead_tuple[0]
            current_culture_dead_time = all_culture_dead_tuple[1]
            current_time = current_culture_dead_time.replace('culture_dead_', '')
            tmp_df = df[df[current_time]==0]
            tmp_df = tmp_df[tmp_df[previous_time]>threshold]
            add_taxon_labels(previous_time, taxon_col, tmp_df, current_culture_dead_time, node_values, new_links, abundance_threshold)

    else:
        tmp_df = df[df[final_time]>threshold]
        new_links = {}
        add_taxon_labels(previous_time, taxon_col, tmp_df, final_time, node_values, new_links, abundance_threshold)

    all_edges = sankey_links | new_links

    sources, source_label_nodes, targets, target_label_nodes, labels, values = generate_node_edges(all_edges, node_values)

    nodes = [node_nb for node_nb in range(len(labels))]
    node_values = [node_values[node_label] for node_label in labels]
    sankey_node_df = pd.DataFrame({'node': nodes, 'label': labels, 'value': node_values})
    node_output_file = os.path.join(output_folder, 'sankey_diagram_node.tsv')
    sankey_node_df.to_csv(node_output_file, sep='\t', index=False)

    sankey_df = pd.DataFrame({'source': sources, 'source_node': source_label_nodes, 'target': targets, 'target_node': target_label_nodes, 'value': values})
    edge_output_file = os.path.join(output_folder, 'sankey_diagram_edge.tsv')
    sankey_df.to_csv(edge_output_file, sep='\t', index=False)

    fig = go.Figure(data=[go.Sankey(
        node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label=labels
        ),
        link = dict(
        source = sources,
        target = targets,
        value = values
    ))])

    fig.update_layout(title_text="", font_size=16, width=1600, height=1000)
    diagram_output_file = os.path.join(output_folder, 'sankey_diagram.png')
    fig.write_image(diagram_output_file)


def handle_input(abundance_file, json_file, taxon_name, output_folder, abundance_threshold=0, relative=None):
    """ Parse json input file and gives inputs to generate_sankey_diagram.

    Args:
        abundance_file (str): path to input abundance file
        json_file (str): path to json file containing dictionary associating samples and time steps
        output_folder (str): path to output folder
        abundance_threshold (int): abundance htreshold to show taxon name
        relative (bool): if True, compute relative abundance instead of absolute
    """
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Read json file containing experiment of enrichment cultures.
    with open(json_file, 'r') as open_json_file:
        json_data = json.load(open_json_file)

    for experiment in json_data:
        logger.info(f'Generate Sankey diagram for experiments {experiment}')
        type_cols = json_data[experiment]
        experiment_output_folder = os.path.join(output_folder, experiment)
        generate_sankey_diagram(abundance_file, type_cols, taxon_name, experiment_output_folder, abundance_threshold, relative)
