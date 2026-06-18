from sankeymien.sankey import handle_input, generate_sankey_diagram, generate_node_edges
import shutil
import os

def test_handle_input():
    abundance_file = 'abundance_table.tsv'
    json_file = 'experiments.json'
    taxon_col = 'taxon'
    output_folder = 'output_folder'
    handle_input(abundance_file, json_file, taxon_col, output_folder)

    sankey_diagram = os.path.join(output_folder, 'experiment_1', 'sankey_diagram.png')
    assert os.path.exists(sankey_diagram)
    sankey_diagram = os.path.join(output_folder, 'experiment_2', 'sankey_diagram.png')
    assert os.path.exists(sankey_diagram)

    shutil.rmtree(output_folder)


def test_generate_node_edges():
    all_edges = {('node_a', 'node_b'): 0.5, ('node_a', 'node_c'): 0.25, ('node_b', 'node_c'): 0}
    node_values = {'node_a': 0.75, 'node_b': 0.5, 'node_c': 0.25, 'node_d': 0}

    sources, source_label_nodes, targets, target_label_nodes, labels, values = generate_node_edges(all_edges, node_values)

    assert sources == [0, 0,]
    assert source_label_nodes == ['node_a', 'node_a']
    assert targets == [1, 2]
    assert target_label_nodes == ['node_b', 'node_c']
    assert labels == ['node_a', 'node_b', 'node_c']
    assert values == [0.5, 0.25]

