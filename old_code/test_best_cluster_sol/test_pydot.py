import pydot

graph = pydot.Dot("my_graph", graph_type="graph")

# Add nodes
my_node = pydot.Node("a", label="Root", xlabel="t=1") # Instead of foo use the node numbers in the graph
graph.add_node(my_node)
# Or, without using an intermediate variable:
graph.add_node(pydot.Node("b", shape="circle", xlabel="t=2", fixedsize="true", width=0.3)) # use the xlabels to set the timeline. Change node sizes according to # of cells.

# Add edges
my_edge = pydot.Edge("a", "b", color="blue", label="1-12") # Here label can be used as mutations
graph.add_edge(my_edge)

graph.add_node(pydot.Node("c", shape="circle", xlabel="t=3"))
graph.add_node(pydot.Node("d", shape="circle"))

# Or, without using an intermediate variable:
graph.add_edge(pydot.Edge("b", "c", color="blue"))

graph.add_edge(pydot.Edge("b", "d", style="dotted"))

graph.write_png("output.png")
