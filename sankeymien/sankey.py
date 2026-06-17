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

def diff_between_group(abundance_dataframe, positive_control_groups, negative_control_groups, group_abundance, threshold):
    tmp_df = abundance_dataframe.copy()
    for group in positive_control_groups:
        tmp_df = tmp_df[tmp_df[group]>threshold]
    for group in negative_control_groups:
        tmp_df = tmp_df[tmp_df[group]==0]
    abundance_org = tmp_df[group_abundance].sum()
    return abundance_org


def add_taxon_labels(sources, targets, labels, values, new_nodes, len_intial_values, selected_col, taxon_col, tmp_df, source_node):
    """ Add taxon name as nodes in Sankey diagram.

    Args:
        sources (list): list of source node number in Sankey diagram
        targets (list): list of source target number in Snakey diagram
        values (list): values for the edge between source and target lists
        new_nodes (dict): dictionary linking step/organism names and node number
        len_intial_values (int): len of list of intial time steps/nodes
        selected_col (str): column that will be used to compute the abundance of the organisms
        taxon_col (str): column containing taxon names
        tmp_df (pd.DataFrame): filtered dataframe containing abundance and taxon names
        source_node (str): source node number from which the edge will start
    """
    for org in tmp_df[taxon_col].unique():
        org_df = tmp_df[tmp_df[taxon_col]==org]
        org_abund = org_df[selected_col].sum()
        if org_abund > 0:
            if org not in new_nodes:
                new_nodes[org] = len(new_nodes) + len_intial_values
            node_nb = new_nodes[org]
            sources.append(source_node)
            targets.append(node_nb)
            values.append(org_abund)
            labels.append(org)


def generate_sankey_diagram(input_file, type_cols, taxon_col, output_file):
    df = pd.read_csv(input_file, sep='\t', index_col=0)
    df[taxon_col] = df[taxon_col].replace(np.nan, 'Unknown')
    samples_cols = [sample for res_type in type_cols for sample in type_cols[res_type]]
    cavern_name_to_ids = samples_cols + [taxon_col]
    df = df[cavern_name_to_ids]

    for col in df.columns:
        if col in samples_cols:
            df[col] = df[col] / df[col].sum()

    threshold = 0
    for type_col in type_cols:
        df[type_col] = df[type_cols[type_col]].mean(axis=1)

    abund_org_unique_kc_T1 = diff_between_group(df, ['kit_control', 'enrichment_T01'], ['initial'], 'enrichment_T01', threshold)

    abund_org_unique_reservoir_T1 = diff_between_group(df, ['initial', 'enrichment_T01'], ['kit_control'], 'enrichment_T01', threshold)

    abund_org_shared_kc_reservoir_T1 = diff_between_group(df, ['initial', 'kit_control', 'enrichment_T01'], [], 'enrichment_T01', threshold)

    negative_controls = ['kit_control', 'initial']
    if 'kit_control_T01' in type_cols:
        negative_controls.append('kit_control_T01')
    abund_org_appearing_T1 = diff_between_group(df, ['enrichment_T01'], negative_controls, 'enrichment_T01', threshold)

    if 'kit_control_T01' in type_cols:
        abund_org_kc_T1 = diff_between_group(df, ['kit_control_T01', 'enrichment_T01'], ['kit_control', 'initial'], 'enrichment_T01', threshold)
    else:
        abund_org_kc_T1 = None

    abund_org_dead_kc = diff_between_group(df, ['kit_control'], ['enrichment_T01'], 'kit_control', threshold)

    abund_org_dead_reservoir = diff_between_group(df, ['initial'], ['enrichment_T01'], 'initial', threshold)

    sankey_nodes = ['kit_control', 'initial', 'kc_initial', 'kit_control_T01', 'enrichment_T01', 'enrichment_T02']
    sankey_links = {('kit_control', 'enrichment_T01'): abund_org_unique_kc_T1, ('kc_initial', 'enrichment_T01'): abund_org_shared_kc_reservoir_T1,
                    ('initial', 'enrichment_T01'): abund_org_unique_reservoir_T1,
                    ('appearing', 'enrichment_T01'): abund_org_appearing_T1, ('kit_control_T01', 'enrichment_T01'): abund_org_kc_T1,
                   ('kit_control', 'inoc_dead'): abund_org_dead_kc, ('initial', 'inoc_dead'): abund_org_dead_reservoir}

    enrichment_times = sorted([type_col for type_col in type_cols if 'enrichment_' in type_col and type_col != 'enrichment_T01'])

    previous_time =  'enrichment_T01'
    for enrichment_time in enrichment_times:
        current_time = enrichment_time
        culture_dead = 'culture_dead_' + current_time
        kit_control = current_time.replace('enrichment_', 'kit_control_')
        appearing = 'appearing_' + current_time
        if enrichment_time in type_cols:
            if kit_control in type_cols:
                abund_org_kc_T1 = diff_between_group(df, [kit_control, current_time], ['kit_control', 'initial'], current_time, threshold)
            else:
                abund_org_kc_T1 = None
            abund_org_T1_T2 = diff_between_group(df, [previous_time, current_time], [], previous_time, threshold)
            abund_org_appearing_T2 = diff_between_group(df, [current_time], [previous_time], current_time, threshold)
            abund_org_dead_T1 = diff_between_group(df, [previous_time], [current_time], previous_time, threshold)
            negative_controls = ['kit_control', 'initial', previous_time]
            abund_org_appearing_T1 = diff_between_group(df, [current_time], negative_controls, current_time, threshold)
            sankey_links[(previous_time, current_time)] = abund_org_T1_T2
            sankey_links[(kit_control, current_time)] = abund_org_kc_T1
            sankey_links[(previous_time, culture_dead)] = abund_org_dead_T1
            sankey_links[(appearing, current_time)] = abund_org_appearing_T1
            previous_time = current_time

    nb_nodes = {}
    labels = []
    sources = []
    targets = []
    values = []
    node_values = {}
    already_processed_tuples = []
    all_nodes = [node for node_tuples in sankey_links for node in node_tuples]
    for sankey_node in all_nodes:
        if sankey_node in type_cols:
            node_links = [sankey_link for sankey_link in sankey_links if sankey_node in sankey_link]
            for node_link in node_links:
                if sankey_links[node_link] is not None and node_link not in already_processed_tuples:
                    if node_link[0] not in nb_nodes:
                        nb_nodes[node_link[0]] = len(nb_nodes)
                        labels.append(node_link[0])
                    if node_link[1] not in nb_nodes:
                        nb_nodes[node_link[1]] = len(nb_nodes)
                        labels.append(node_link[1])
                    sources.append(nb_nodes[node_link[0]])
                    targets.append(nb_nodes[node_link[1]])
                    values.append(sankey_links[node_link])
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

    intial_values = values

    if final_time != 'enrichment_T01':
        tmp_df = df[df[final_time]>threshold]
        new_nodes = {}
        modified_values = intial_values.copy()
        len_intial_values = len(intial_values)+1
        add_taxon_labels(sources, targets, labels, modified_values, new_nodes, len_intial_values, final_time, taxon_col, tmp_df, nb_nodes[final_time])
        all_culture_dead_tuples = [node_tuples for node_tuples in sankey_links for node in node_tuples if 'culture_dead_' in node]
        for all_culture_dead_tuple in all_culture_dead_tuples:
            previous_time = all_culture_dead_tuple[0]
            current_culture_dead_time = all_culture_dead_tuple[1]
            current_time = current_culture_dead_time.replace('culture_dead_', '')
            tmp_df = df[df[current_time]==0]
            tmp_df = tmp_df[tmp_df[previous_time]>threshold]
            add_taxon_labels(sources, targets, labels, modified_values, new_nodes, len_intial_values, previous_time, taxon_col, tmp_df, nb_nodes[current_culture_dead_time])

    else:
        tmp_df = df[df[final_time]>threshold]
        new_nodes = {}
        modified_values = intial_values.copy()
        len_intial_values = len(intial_values)
        add_taxon_labels(sources, targets, labels, modified_values, new_nodes, len_intial_values, previous_time, taxon_col, tmp_df, nb_nodes[final_time])

    source_label_nodes = [labels[source] for source in sources]
    target_label_nodes = [labels[target] for target in targets]
    sankey_df = pd.DataFrame({'source': sources, 'source_node': source_label_nodes, 'target': targets, 'target_node': target_label_nodes, 'value': modified_values})
    sankey_df.to_csv(output_file.replace('.png', '.tsv'), sep='\t', index=False)

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
        value = modified_values
    ))])

    fig.update_layout(title_text="", font_size=16, width=1600, height=1000)
    fig.write_image(output_file)


def handle_input(abundance_file, json_file, taxon_name, output_folder):
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Read json file containing experiment of enrichment cultures.
    with open(json_file, 'r') as open_json_file:
        json_data = json.load(open_json_file)

    for experiment in json_data:
        logger.info(f'Generate Sankey diagram for experiments {experiment}')
        type_cols = json_data[experiment]
        output_file = os.path.join(output_folder, experiment+'.png')
        generate_sankey_diagram(abundance_file, type_cols, taxon_name, output_file)
