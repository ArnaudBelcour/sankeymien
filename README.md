# Sankey diagram for microbial enrichment culture

Python package to make Sankey diagram for Microbial Enrichment cultures.

![](pictures/experiment_1.png)

## Installation

This package requires:
- plotly and kaleido (with google chrome installed).
- pandas and numpy.

It can be installed by cloning this repository and with pip:
```
pip isntall -e .
```

## Usage

It can be used as a command line:

```
sankeymien -a abundance_table.tsv -j experiments.json -o testoutput -t taxon
```

`-a abundance_table.tsv` is a tabulated file containing as a first column organism ID, several columns corresponding to samples anda column containing taxon name (and given with `-t taxon`).

`-j experiments.json` is a json file indicating the experiments and the association between time steps and samples.

`-o testoutput` is the output folder.

## Usage

It generates a Sankey diagram showing the flow of relative abundance of organisms according to the different time steps of the enrichment cultures. The final enrichment is linked to taxon name to show the composition of the final communities.
