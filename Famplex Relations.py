import pybel
import urllib.request
import pandas as pd


######Importing proteins and protein complexes######

from pybel.dsl import protein, named_complex_abundance
from pybel.constants import HAS_MEMBER


######Accesing famplex relations file######

fplx = 'https://raw.githubusercontent.com/sorgerlab/famplex/master/relations.csv'
df = pd.read_csv(fplx, header=None)


######Representing relations as BEL statements######

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


with open("famplex.bel", "w") as file:
    pybel.to_bel(graph, file)
