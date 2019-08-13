# -*- coding: utf-8 -*-

"""Manager for Bio2BEL FamPlex."""

from typing import Mapping

import networkx as nx
from tqdm import tqdm

from bio2bel.manager.bel_manager import BELManagerMixin
from pybel import BELGraph
from pybel.dsl import Protein
from .equivalences import append_equivalences_graph, get_equivalences_df
from .relations import build_relations_graph, get_relations_df

__all__ = [
    'Manager',
]


class Manager(BELManagerMixin):
    """Protein family and complex hierarchy."""

    def __init__(self, *args, **kwargs):  # noqa:D107
        relations_df = get_relations_df()
        self.graph = build_relations_graph(relations_df)
        equivalences_df = get_equivalences_df()
        append_equivalences_graph(equivalences_df, self.graph)

    @classmethod
    def _get_connection(cls):
        pass

    @staticmethod
    def is_populated() -> bool:
        """Return if the database is populated."""
        return True

    def summarize(self) -> Mapping[str, int]:
        """Summarize the database."""
        return dict(
            relations=self.count_relations(),
            entries=self.count_entries(),
        )

    def count_entries(self) -> int:
        """Count the number of entries in the database."""
        return self.graph.number_of_nodes()

    def count_relations(self) -> int:
        """Count the number of relationships in the database."""
        return self.graph.number_of_edges()

    def to_bel(self) -> BELGraph:
        """Generate a BEL graph."""
        return self.graph

    def normalize_terms(self, graph: BELGraph, use_tqdm: bool = False) -> None:
        """Normalize FamPlex nodes in the graph."""
        it = graph.nodes()
        if use_tqdm:
            it = tqdm(it, total=graph.number_of_nodes())
        m = {
            node: Protein(namespace=node.namespace, name=node.name, identifier=node.name)
            for node in it
            if isinstance(node, Protein) and node.namespace.upper() == 'FPLX'
        }
        nx.relabel_nodes(graph, m, copy=False)
