#!/usr/bin/env python3
import sys
from os.path import dirname, abspath

sys.path.append(abspath(dirname(__file__) + "/.."))

if __name__ == "__main__":
    from graph.graph_structs import MetaGraph, EdgeTag

    # We will first create a graph with structure but without periodicity information
    graph = MetaGraph()

    # Edge tag for no period crossing for simple checks
    bare_tag = EdgeTag(0, 0, 0)

    # Add some edges to create loops and leaves/twigs
    graph.add_edge(0, 2, bare_tag)
    graph.add_edge(0, 1, bare_tag)
    graph.add_edge(0, 3, bare_tag)
    graph.add_edge(2, 3, bare_tag)
    graph.add_edge(1, 4, bare_tag)
    graph.add_edge(1, 5, bare_tag)
    graph.add_edge(5, 6, bare_tag)
    graph.add_edge(7, 8, bare_tag)
    graph.add_edge(8, 9, bare_tag)
    graph.add_edge(7, 9, bare_tag)
    graph.add_edge(0, 10, bare_tag)
    graph.add_edge(11, 10, bare_tag)
    graph.add_edge(12, 11, bare_tag)
    graph.add_edge(10, 9, bare_tag)

    # Test node-type detection
    graph.mark_states()
    graph.dump("test.dot")

    # Test copying of graph
    other_graph = graph.copy()
    other_graph.dump("cloned.dot")

    # Test simplification of graph
    other_graph.reduce()
    other_graph.dump("reduced.dot")

    # Test component detection (detects whole molecules in your scenario :) )
    other_graph.find_components()
    other_graph.dump("components.dot")

    """# Detect, how many components this graph has and return an array, which of this percolate in in x-dimension (parameter = 0)
    num_components, has_loop = other_graph.find_loops(0)
    print("No loop expected:", num_components, has_loop)"""

    # # Create a graph with actual directed edges incorporating dimension crossings
    print("Construct a realistic graph for loop detection")
    graph = MetaGraph()

    # Edge tag for an edge crossing the x-dimension period upwards
    ux_tag = EdgeTag(1, 0, 0)

    # Edge tag for an edge not crossing x-dimension period at all
    nx_tag = EdgeTag(0, 0, 0)

    # Edge tag for an edge crossing the x-dimension period downwards
    dx_tag = EdgeTag(-1, 0, 0)

    # Create more elaborate period-crossing patterns
    graph.add_edge(0, 2, ux_tag)
    graph.add_edge(0, 1, ux_tag)
    graph.add_edge(0, 3, nx_tag)
    graph.add_edge(2, 3, dx_tag)
    graph.add_edge(1, 4, nx_tag)
    graph.add_edge(1, 5, nx_tag)
    graph.add_edge(5, 6, dx_tag)
    graph.add_edge(7, 8, dx_tag)
    graph.add_edge(8, 9, ux_tag)
    graph.add_edge(7, 9, dx_tag)
    graph.add_edge(0, 10, dx_tag)
    graph.add_edge(11, 10, dx_tag)
    graph.add_edge(12, 11, ux_tag)
    graph.add_edge(10, 9, ux_tag)

    graph.add_edge(13, 14, dx_tag)
    graph.add_edge(14, 15, nx_tag)
    graph.add_edge(15, 13, ux_tag)

    graph.add_metadata(0, "hello :)")

    graph.dump("raw_tagged.dot", True)

    # Now simplify the more elaborate graph
    graph.reduce()
    graph.dump("tagged.dot", True)

    """# Detect y loops. There shouldn't be any
    num_components, has_loop = graph.find_loops(1)
    print("No loop expected:", num_components, has_loop)

    # Detect x loops. There should be one
    num_components, has_loop = graph.find_loops(0)
    print("Loop expected (7,8,9):", num_components, has_loop)

    # Below this is the sample output for loop detection in all 3 dimensions.
    print(
        "To detect whether a component percolates in all three dimensions, run all three find-loop operations and combine the results for individual components"
    )
    print("Find loop in x: ", graph.find_loops(0))
    print("Find loop in y: ", graph.find_loops(1))
    print("Find loop in z: ", graph.find_loops(2))"""

    # New rotationally-invariant loop check:

    print(
        "To detect whether a component percolates in all three dimensions in a rotationally invariant/stable way, check it using find_stable_loops."
    )
    print("Find stable loop dimension: ", graph.find_stable_loops())
    print("Get metadata: ", graph.get_components())

    print("Check duplicate elimination...")

    graph = MetaGraph()

    graph.add_edge(0, 1, ux_tag)
    graph.add_edge(0, 1, ux_tag)
    graph.add_edge(1, 2, nx_tag)
    graph.add_edge(0, 1, nx_tag)
    graph.add_edge(0, 2, nx_tag)
    graph.add_edge(0, 2, nx_tag)

    graph.dump("duplicates_removed.dot", True)

    graph_comp = graph.get_component_graph()

    graph_comp.dump("duplicates_removed_comp_graph.dot", True)