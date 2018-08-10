#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convert famplex relations to BEL Statements."""

import pandas as pd
import pybel
from pybel.constants import HAS_MEMBER, NAMESPACE, NAME
from pybel.dsl import named_complex_abundance, protein

from bio2bel_famplex.constants import RELATIONS_URL

__all__ = ["enrich_graph", "build_graph"]

def get_df():
    """Get famplex relations as a Pandas dataframe.

    :rtype: pandas.DataFrame
    """
    return pd.read_csv(RELATIONS_URL, header=None)


def enrich_graph(graph):
    """ Find famplexes in the graph and enrich them"""
    df = get_df()
    for node in list(graph):
        if is_famplex_node(graph, node, df):
            subdf = look_up(graph, node, df)
            append_graph(subdf, graph)


def is_famplex_node(graph, node, df):
    data = graph.node[node]
    namespace = data.get(NAMESPACE)
    if namespace is None:
        return False

    return namespace == "FPLX"


def look_up(graph, node, df):
    """ Get a subset of a dataframe relevant to this node"""
    data = graph.node[node]
    name = data.get(NAME)
    subdf = df[(df[1] == name) | (df[4] == name)]
    return subdf


def build_graph(df):
    """Build a BEL Graph from a famplex dataframe.

    :param df: A Pandas dataframe representing famplex relations
    :type df: pandas.DataFrame
    :rtype: pybel.BELGraph
    """
    graph = pybel.BELGraph(name="Famplex", version="0.01", authors="Kristian Kolpeja")
    append_graph(df, graph)
    return graph


def append_graph(df, graph):
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


def main():
    """Provide BEL Statements from BEL graph."""
    df = get_df()
    graph = build_graph(df)
    with open("famplex.bel", "w") as file:
        pybel.to_bel(graph, file)


if __name__ == '__main__':
    main()
