#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convert FamPlex equivalences to BEL Statements."""

import sys

import click
import pandas as pd

import pybel
from bio2bel_famplex.constants import EQUIVALENCES_URL
from pybel import BELGraph
from pybel.dsl import named_complex_abundance, protein


def get_df() -> pd.DataFrame:
    """Get FamPlex relations as a Pandas dataframe."""
    return pd.read_csv(EQUIVALENCES_URL, header=None)


famplex_to_identifiers = {
    "ECCODE": ("eccode", r"^\d+\.-\.-\.-|\d+\.\d+\.-\.-|\d+\.\d+\.\d+\.-|\d+\.\d+\.\d+\.(n)?\d+$"),
    "GO": ("go", "^GO:\d{7}$"),
    "IP": ("interpro", "^IPR\d{6}$"),
    "MESH": ("mesh", "^(C|D)\d{6}$"),
    "NCIT": ("ncit", "^C\d+$"),
    "PF": ("pfam", "^PF\d{5}$"),
    "RE": ("reactome", "(^R-[A-Z]{3}-\d+(-\d+)?(\.\d+)?$)|(^REACT_\d+(\.\d+)?$)")
}

famplex_to_belns = {
    "SFAM": "https://arty.scai.fraunhofer.de/artifactory/bel/namespace/selventa-protein-families/selventa-protein-families-20170725.belns",
    "SCOMP": "https://arty.scai.fraunhofer.de/artifactory/bel/namespace/selventa-named-complexes/selventa-named-complexes-20170725.belns",
    "FPLX": "https://raw.githubusercontent.com/sorgerlab/famplex/1b7e14ec0fd02ee7ed71514c6e267f57d5641a4b/export/famplex.belns"
}


def build_graph(df: pd.DataFrame) -> BELGraph:
    """Build a BEL Graph from a famplex dataframe.

    :param df: A Pandas dataframe representing famplex equivalences
    """
    graph = BELGraph(name="Famplex_Equivalences", version="0.0.1", authors="Kristian Kolpeja")
    graph.namespace_url.update(famplex_to_belns)

    graph.namespace_pattern.update({
        identifiers_key: pattern
        for _, (identifiers_key, pattern) in famplex_to_identifiers.items()
    })

    for index, (namespace, name, famplex_name) in df.iterrows():
        is_complex = "complex" in name.lower()

        if is_complex:
            func = named_complex_abundance
        else:
            func = protein

        if namespace in famplex_to_belns:
            if is_complex:
                namespace = "SCOMP"
            else:
                namespace = "SFAM"

        elif namespace in famplex_to_identifiers:
            namespace = famplex_to_identifiers[namespace][0]
        else:
            continue

        external = func(namespace, name)
        famplex = func('FPLX', famplex_name)

        graph.add_equivalence(external, famplex)

    return graph


@click.command()
@click.option("-f", "--file", type=click.File("w"), help="A file to write famplex equivalences to BEL",
              default=sys.stdout)
def main(file):
    """Provide BEL Statements from BEL graph."""
    df = get_df()
    graph = build_graph(df)
    pybel.to_bel(graph, file)


if __name__ == '__main__':
    main()
