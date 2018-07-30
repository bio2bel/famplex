#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convert famplex relations to BEL Statements."""

import pandas as pd
import pybel
from pybel.constants import HAS_MEMBER
from pybel.dsl import protein, named_complex_abundance

FPLX = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'


def get_df():
    """Get famplex relations as a Pandas dataframe.

    :rtype: pandas.DataFrame
    """
    return pd.read_csv(FPLX, header=None)


def build_graph(df):
    """Build a BEL Graph from a famplex dataframe.

    :param df: A Pandas dataframe representing famplex relations
    :type df: pandas.DataFrame
    :rtype: pybel.BELGraph
    """
    graph = pybel.BELGraph(name="Famplex", version="0.01", authors="Kristian Kolpeja")
    graph.namespace_url.update({
        "HGNC": "https://arty.scai.fraunhofer.de/artifactory/bel/namespace/hgnc/hgnc-20180215.belns",
        "FPLX": "https://raw.githubusercontent.com/sorgerlab/famplex/1b7e14ec0fd02ee7ed71514c6e267f57d5641a4b/export/famplex.belns"
    })

    for index, (ns1, n1, rel, ns2, n2) in df.iterrows():
        if ns1 == "UP":
            continue

        sub = protein(ns1, n1)

        if rel == 'isa':
            obj = protein(ns2, n2)
            graph.add_is_a(sub, obj)

        else:  # partof
            obj = named_complex_abundance(ns2, n2)
            graph.add_unqualified_edge(obj, sub, HAS_MEMBER)
    return graph


def main():
    """Provide BEL Statements from BEL graph."""
    df = get_df()
    graph = build_graph(df)
    with open("famplex.bel", "w") as file:
        pybel.to_bel(graph, file)


if __name__ == '__main__':
    main()
