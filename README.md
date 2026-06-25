# Sankey diagram for microbial enrichment culture

Python package to make Sankey diagram for Microbial Enrichment cultures.

![](pictures/experiment_1.png)

## Installation

This package requires:
- [plotly](https://github.com/plotly/plotly.py) and [kaleido](https://github.com/plotly/Kaleido) (with google chrome installed) to generate Sankey diagrams.
- [pandas](https://github.com/pandas-dev/pandas) and [numpy](https://github.com/numpy/numpy) to handle dataframe and filter organisms according to abundance in samples.

Kaleido requires Google Chrome to be isntalled, which can be achieved with:

```
kaleido_get_chrome
```

Sankeymien and the python package dependencies can be installed with pip:
```
pip install sankeymien
```

It can be also installed by cloning this repository and with pip:
```
git clone https://github.com/ArnaudBelcour/sankeymien.git

cd sankeymien

pip install -e .
```

## Usage

The package can be launched using the command line:

```
sankeymien -a abundance_table.tsv -j experiments.json -o testoutput -t taxon
```

`-a abundance_table.tsv` is a tabulated file containing as a first column organism ID, several columns corresponding to samples anda column containing taxon name (and given with `-t taxon`).

`-j experiments.json` is a json file indicating the experiments and the association between time steps and samples.

`-o testoutput` is the output folder.

`--relative` to use relative abundance instead of the abundance from the input files. This takes the abundance of an organism in a sample and it divides it by the total abundance in the sample.

`--abundance-threshold` abundance threshold used to show the edges and the associated taxa in the enrichment cultures.

## Input

The abundance sample should look like this:

| organism | s1   | s2   | k1 | Ec1_1 | Ec1_2 | Ec1_3 | kec1 | taxon            |
|----------|------|------|----|-------|-------|-------|------|------------------|
| org_a    | 1    | 2    | 5  | 0     | 0     | 0     | 0    | Bos              |
| org_b    | 2000 | 3000 | 10 | 100   | 120   | 130   | 1    | Escherichia      |
| org_c    | 500  | 600  | 5  | 100   | 200   | 300   | 0    | Parcubacteria    |
| org_d    | 20   | 30   | 0  | 25    | 35    | 40    | 0    | Methanobacterium |

The json should correspond to this:

```json
{
    "experiment_1": {
        "initial": ["s1", "s2"],
        "kit_control": ["k1"],
        "enrichment_T01": ["Ec1_1", "Ec1_2", "Ec1_3"],
        "kit_control_T01": ["kec1"]
    }
}
```

## Output

It generates one output folder for each experiment in the json file. The output folder contains:

```
output_folder
├── experiment_1
│   ├── intermediary_folder
│   |   └── ...
│   |   └── edge_initial_to_enrichment_T01.tsv
│   |   └── sample_1_nitrogen_cycle.png
│   |   └── ...
│   ├── sankey_diagram_edge.tsv
│   ├── sankey_diagram_node.tsv
│   ├── sankey_diagram.png
├── sankeymien.log
├── sankeymien_metadata.json
```

`intermediary_folder` contains multiple tabulated file corresponding to sub part of the initial dataframe filtered to show only organisms ocurring in the associated edges of the Sankey diagram.

`sankey_diagram.png` a Sankey diagram showing the flow of relative abundance of organisms according to the different time steps of the enrichment cultures. The final enrichment is linked to taxon name to show the composition of the final communities.

`sankey_diagram_edge.tsv` indicates the different edges of the Sankey diagram.

`sankey_diagram_node.tsv` indicates the different nodes of the Sankey diagram.

`sankeymien.log` is a log file showing printed lines.

`sankeymien_metadata.json` contains metadata linked to the run: version of the different packages used.