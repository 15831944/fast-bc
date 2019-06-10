import geog
import networkx as nx
import osmgraph

# By default any way with a highway tag will be loaded
gmap = osmgraph.parse_file('/home/davide/Downloads/map(2).osm') # or .osm or .pbf
g = nx.convert_node_labels_to_integers(gmap)
for n1, n2 in g.edges():
    c1, c2 = osmgraph.tools.coordinates(g, (n1, n2))   
    print n1, n2, (int(round(geog.distance(c1, c2))) + 1)

