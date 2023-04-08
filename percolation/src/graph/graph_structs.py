import typing
from copy import deepcopy

from enum import Enum
import numpy as np
from fractions import gcd


class EdgeTag:
    """Class representing periodicity information for an associated edge.
    Provides the edge with information, in which dimension it crosses the periodicity in which direction.

    Attributes:
        dx (int): Direction and number of crossings of the x-periodicity
        dy (int): Direction and number of crossings of the y-periodicity
        dz (int): Direction and number of crossings of the z-periodicity

    """

    def __init__(self, dx, dy, dz):
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def get_vec(self):
        """
        Obtain a representation of the Edge tag as a numpy array
        """
        return np.array([self.dx, self.dy, self.dz])

    def append(self, other):
        """
        Create a new EdgeTag from the concatenation of this and the other edge

        Args:
            self (EdgeTag): This Edge Tag
            other (EdgeTag): The other Edge Tag to append to this

        """
        return EdgeTag(other.dx + self.dx, other.dy + self.dy, other.dz + self.dz)

    def calculate_connecting(self, target):
        """
        Create a new EdgeTag that could connect the two endpoints of this and the other edge tag (subtraction)

        Args:
            self (EdgeTag): This Edge Tag
            target (EdgeTag): The other Edge Tag that the connecting tag should be calculated to

        """
        return EdgeTag(target.dx - self.dx, target.dy - self.dy, target.dz - self.dz)

    def invert(self):
        """
        Create an EdgeTag with the inverse direction of the current Tag

        Args:
            self (EdgeTag): This Edge Tag

        """
        return EdgeTag(-self.dx, -self.dy, -self.dz)

    def normalized(self):
        """
        Create an EdgeTag with the dimensions' respective values reduced by possible common factors

        Args:
            self (EdgeTag): This Edge Tag

        """
        normalizer = self.dx
        normalizer = gcd(normalizer, self.dy)
        normalizer = gcd(normalizer, self.dz)
        return EdgeTag(
            int(self.dx / normalizer),
            int(self.dy / normalizer),
            int(self.dz / normalizer),
        )

    def __getitem__(self, key):
        """
        Allow for array-like access to this tag's entries.

        Args:
            self (EdgeTag): This Edge Tag
            key (int): The key of the dimension to be accessed

        """
        if key < 0 or key > 2:
            return None

        data = {0: self.dx, 1: self.dy, 2: self.dz}

        return data[key]

    def __setitem__(self, key, item):
        """
        Allow for array-like setting of this tag's entries.

        Args:
            self (EdgeTag): This Edge Tag
            key (int): The key of the dimension to be accessed
            item (int): value to assign to that dimension

        """
        if key < 0 or key > 2:
            raise Exception("invalid key {0}".format(key))

        if key == 0:
            self.dx = item
        elif key == 1:
            self.dy = item
        else:
            self.dz = item

    def __eq__(self, obj):
        """
        Custom equal comparison routine

        Args:
            self (EdgeTag): This Edge Tag
            obj (undefined): The other object to compare this to

        """
        return (
            isinstance(obj, EdgeTag)
            and self.dx == obj.dx
            and self.dy == obj.dy
            and self.dz == obj.dz
        )

    def __ne__(self, obj):
        """
        Custom not-equal comparison routine

        Args:
            self (EdgeTag): This Edge Tag
            obj (undefined): The other object to compare this to

        """
        return not self.__eq__(obj)

    def __str__(self):
        """
        Obtain string representation of tag

        Args:
            self (EdgeTag): This Edge Tag

        """
        return "[[{0}|{1}|{2}]]".format(self.dx, self.dy, self.dz)

    def __repr__(self):
        """
        Obtain string representation of tag

        Args:
            self (EdgeTag): This Edge Tag

        """
        return "[[{0}|{1}|{2}]]".format(self.dx, self.dy, self.dz)


class EdgeEntry:
    """
    Structure to hold edge index, neighbor index and edge tag together. The neighbor is the other node this edge points to.
    Used as an entry in an adjacency list representing one edge

    Attributes:
        index (int): The index of this edge in the graph
        neighbor (int): The index of the neighbor in the graph that is accessible via this edge
        tag (EdgeTag): The tag associated with this connecting edge

    Args:
        index (int): The index of this edge in the graph
        neighbor (int): The index of the neighbor in the graph that is accessible via this edge
        tag (EdgeTag): The tag associated with this connecting edge

    """

    def __init__(self, index, neighbor, tag):
        self.index = index
        self.neighbor = neighbor
        self.tag = tag

    def __eq__(self, value):
        return (
            isinstance(value, EdgeEntry)
            and self.neighbor == value.neighbor
            and self.tag == value.tag
        )

    def __hash__(self):
        return hash(self.neighbor)


class NodeState(Enum):
    """
    An enumerable type that denotes a type or sate of node

    The possible values are:
    DEFAULT: All nodes start in this state
    DISABLED: This state may be obtained after reduce() is performed and a node has been marked for removal on copy
    LEAF: The actual leaf nodes get this label
    TWIG: Linear sequences of nodes connecting a leaf to another component (a strongly connected component or another leaf)
    ISOLATED: Node not connected to any other node in the graph

    Inheritance:
        Enum:

    """

    DEFAULT = 0
    DISABLED = 1
    LEAF = 2
    TWIG = 3
    ISOLATED = 4


# Function to check linear independence and create a normalized basis set
def make_basis(basis, new_entry):
    """
    Check if a new_entry is linearly independent of previous entries to a basis and extend the basis set if necessary/possible.
    Returns the resulting basis set with a possible normalized version of the new entry candidate appended

    Args:
        basis (set): Previous/current basis set (may be empty)
        new_entry (EdgeTag): The EdgeTag whose independence is to be checked

    """
    zero_tag = EdgeTag(0, 0, 0)
    # Zero vector is not independent:
    if new_entry == zero_tag:
        return basis

    # normalize vectors. May be bad?
    # TODO: FIXME if required
    new_entry = new_entry.normalized()

    is_independent = len(basis) < 3

    if is_independent:
        # Build matrix for rank check:

        # Obtain numpy representation
        data = []
        for old in basis:
            data.append(old.get_vec())
        data.append(new_entry.get_vec())

        # Transform to matrix
        A = np.row_stack(data)
        # print(A)

        # get SVD
        U, s, V = np.linalg.svd(A)

        TOLERANCE = 1e-14

        # get number of non-zero eigenvalues
        rank = int(np.sum(s > TOLERANCE))
        # print(rank, len(basis) + 1)

        is_independent = rank == len(basis) + 1

    if is_independent:
        basis.append(new_entry)

    return basis


class MetaGraph:
    """Class encapsulating the meta-graph of a system
    As far as graph theory is concerned, it is an undirected graph denoted by a directed representation using adjacency lists for each node

    Attributes:
        _adjacency (int): Adjacency information for each node
        _vertex_states (list): Array to keep track of node-specific metadata
        _meta_data (list): Array to keep track of node-specific metadata
        _vertex_component (list): Array to keep track of graph components a vertex is in
        _next_vertex_index (int): Number of currently managed vertices for resizing
        _next_edge_index (int): Number of currently managed edges for loop detection
    """

    # This does some basic setup without initialization

    def __init__(self):
        # Adjacency information
        self._adjacency = []

        # Array to keep track of specific node states
        self._vertex_states = []

        # Array to keep track of node-specific metadata
        self._meta_data = []

        # Array to keep track of graph components a vertex is in
        self._vertex_component = []

        # Index of maximum managed vertex for resizing
        self._next_vertex_index = 0

        # Index of maximum managed edge for loop detection
        self._next_edge_index = 0

    def reserve(self, new_size):
        """
        Function to reserve memory for nodes and edges

        Args:
            self (MetaGraph): The current graph
            new_size (int): The desired minimum number of managed nodes (indexing starts at 0)

        """
        # Pad the adjacency list to account for possibly new vertices
        if self._next_vertex_index < new_size:
            required = new_size - self._next_vertex_index + 1

            for i in range(required):
                self._adjacency.append(set([]))
                self._vertex_states.append(NodeState.DEFAULT)
                self._vertex_component.append(-1)
                self._meta_data.append(None)

            self._next_vertex_index = new_size

    #
    def add_edge(self, from_index, to_index, tag):
        """
        This function creates a new edge from the first (from_index) to the second argument (to_index).
        The tag specifies, which dimensions are crossed.
        The method automatically generates new nodes to accommodate for a possible new highest node index and will
        introduce both orientations of the edge in order to make an effectively undirected graph despite the graph
        internally using a directed representation

        Args:
            self (MetaGraph): The current graph
            from_index (int): Index of the base node of the edge
            to_index (int): Index of the target node of the edge
            tag (EdgeTag): The tag associated with the base to target orientation of the edge

        """
        max_added = max(from_index, to_index)

        min_size = max_added + 1
        self.reserve(min_size)

        # Get new edge index
        new_edge_index = self._next_edge_index

        self._next_edge_index += 1

        # Insert edge in both directions
        self._adjacency[from_index].add(
            EdgeEntry(new_edge_index, to_index, tag))
        self._adjacency[to_index].add(
            EdgeEntry(new_edge_index, from_index, tag.invert())
        )

    def add_metadata(self, node_index, data):
        """
        Method to add metadata to a node of the graph
        Can be retrieved by the get_components() method

        Args:
            self (MetaGraph): The current graph
            node_index (int): Index of the node whose metadata is supposed to change
            data (undefined): The data to assign to the node

        """
        min_size = node_index + 1
        self.reserve(min_size)

        self._meta_data[node_index] = data

    def dump(self, path, show_edge_labels=False):
        """
        Outputs the graph to a file in the DOT format (https://en.wikipedia.org/wiki/DOT_(graph_description_language))

        Args:
            self (MetaGraph): The current graph
            path (string): The path to the desired output file
            show_edge_labels=False (boolean): A boolean flag which will enable the outputting of edge tags as labels in the graph

        """
        with open(path, "w") as out:
            if show_edge_labels:
                out.write("digraph MetaGraph {\n")
            else:
                out.write("graph MetaGraph {\n")
            # Print all vertices
            for i in range(self._next_vertex_index):

                comp_postfix = ""
                if self._vertex_component[i] != -1:
                    comp_postfix = " - C#{0}".format(self._vertex_component[i])

                # Add vertex to graph
                if self._vertex_states[i] == NodeState.DEFAULT:
                    out.write(
                        'P{0} [label="P{0}{1}" color=black];\n'.format(
                            i, comp_postfix)
                    )
                elif self._vertex_states[i] == NodeState.LEAF:
                    out.write(
                        'P{0} [label="P{0}{1}" color=green];\n'.format(
                            i, comp_postfix)
                    )
                elif self._vertex_states[i] == NodeState.TWIG:
                    out.write(
                        'P{0} [label="P{0}{1}" color=brown];\n'.format(
                            i, comp_postfix)
                    )
                elif self._vertex_states[i] == NodeState.ISOLATED:
                    out.write(
                        'P{0} [label="P{0}{1}" color=blue];\n'.format(
                            i, comp_postfix)
                    )
                elif self._vertex_states[i] == NodeState.DISABLED:
                    out.write(
                        'P{0} [label="P{0}{1}" color=red];\n'.format(
                            i, comp_postfix)
                    )

            # Print all edges
            for i in range(self._next_vertex_index):
                if self._vertex_states[i] != NodeState.DISABLED:
                    # Get all connecting edges
                    for entry in self._adjacency[i]:
                        # Filter due to edge being present in both directions
                        if (
                            entry.neighbor
                            >= i
                            # and self._vertex_states[entry.neighbor]
                            #!= NodeState.DISABLED
                        ):
                            # Write edge to output
                            if show_edge_labels:
                                out.write(
                                    'P{0} -> P{1} [label="{2}; {3}; {4}"];\n'.format(
                                        i,
                                        entry.neighbor,
                                        entry.tag[0],
                                        entry.tag[1],
                                        entry.tag[2],
                                    )
                                )
                            else:
                                out.write(
                                    "P{0} -- P{1};\n".format(i, entry.neighbor))

            out.write("}\n")

    def mark_states(self):
        """
        Detect leaf, twig and isolated nodes
        Only marks the nodes as those types. Does not restructure the graph

        Args:
            self (MetaGraph): The current graph

        """
        if self._next_vertex_index < 0:
            return

        neighbor_count = [0] * self._next_vertex_index
        queue = []

        # Detect leaves or singled-out nodes
        for i in range(self._next_vertex_index):
            neighbor_count[i] = len(self._adjacency[i])
            if self._vertex_states[i] != NodeState.DISABLED:
                self._vertex_states[i] = NodeState.DEFAULT

        for i in range(self._next_vertex_index):
            if self._vertex_states[i] != NodeState.DISABLED:
                if neighbor_count[i] == 0:
                    self._vertex_states[i] = NodeState.ISOLATED
                elif neighbor_count[i] == 1:
                    self._vertex_states[i] = NodeState.LEAF
                    # Get neighbor of leaf, reduce neighbor count and add to queue
                    neighbor = None
                    for entry in self._adjacency[i]:
                        neighbor = entry.neighbor
                    queue.append(neighbor)
                    neighbor_count[neighbor] -= 1

        # Iterate over new potential twigs
        while queue:
            node = queue.pop(0)
            # Check if node not yet ignored or
            if self._vertex_states[node] == NodeState.DEFAULT:
                if neighbor_count[node] < 2:
                    self._vertex_states[node] = NodeState.TWIG
                    for edge in self._adjacency[node]:
                        if self._vertex_states[edge.neighbor] == NodeState.DEFAULT:
                            neighbor_count[edge.neighbor] -= 1
                            queue.append(edge.neighbor)

    def reduce(self):
        """
        Mark all twig, leaf and isolated nodes as disabled and remove connecting edges
        This simplifies the graph structure but also changes the graph overall

        Args:
            self (MetaGraph): The current graph

        """
        self.mark_states()
        new_adj = []

        # Detect leaves or singled-out nodes
        for i in range(self._next_vertex_index):
            new_adj.append(set([]))
            if self._vertex_states[i] != NodeState.DEFAULT:
                self._vertex_states[i] = NodeState.DISABLED
                continue

        for i in range(self._next_vertex_index):
            for entry in self._adjacency[i]:
                if self._vertex_states[entry.neighbor] != NodeState.DISABLED:
                    new_adj[i].add(entry)

        self._adjacency = new_adj

    def get_purged_graph(self):
        """
        Remove all disabled nodes to simplify the overall graph

        Args:
            self (MetaGraph): The current graph

        """
        new_index = []

        curr_index = 0

        # Detect leaves or singled-out nodes
        for i in range(self._next_vertex_index):
            if self._vertex_states[i] == NodeState.DISABLED:
                new_index.append(-1)
            else:
                new_index.append(curr_index)
                curr_index += 1

        new_graph = MetaGraph()
        new_graph.reserve(curr_index)

        for curr in range(self._next_vertex_index):
            if self._vertex_states[curr] != NodeState.DISABLED:
                base_index = new_index[curr]
                for edge in self._adjacency[curr]:
                    neighbor = edge.neighbor
                    tag = edge.tag

                    if self._vertex_states[neighbor] != NodeState.DISABLED:
                        neighbor_index = new_index[neighbor]
                        new_graph.add_edge(base_index, neighbor_index, tag)

        return new_graph

    def get_simplified_graph(self):

        is_branching = []
        new_index = []

        curr_index = 0

        neutral_tag = EdgeTag(0, 0, 0)

        encountered_components = set([])

        num_components, members = self.get_components()

        queue = []
        # Detect leaves or singled-out nodes
        for i in range(self._next_vertex_index):
            if self._vertex_states[i] == NodeState.DISABLED:
                new_index.append(-1)
                is_branching.append(False)
            else:
                has_loop = False
                for edge in self._adjacency[i]:
                    if edge.neighbor == i:
                        has_loop = True
                        break

                if has_loop or len(self._adjacency[i]) != 2:
                    new_index.append(curr_index)
                    curr_index += 1
                    is_branching.append(True)
                    encountered_components.add(self._vertex_component[i])
                    queue.append((i, i, -1, neutral_tag))
                else:
                    new_index.append(-1)
                    is_branching.append(False)

        for i in range(num_components):
            if i in encountered_components:
                continue
            else:
                member_index, data = members[i][0]

                # Do not consider disabled components
                if self._vertex_states[i] == NodeState.DISABLED:
                    continue

                new_index[member_index] = curr_index
                curr_index += 1
                is_branching[member_index] = True
                queue.append((member_index, member_index, -1, neutral_tag))

        new_graph = MetaGraph()
        new_graph.reserve(curr_index)

        # Iterate over graph
        while queue:
            current, origin, previous, curr_distance = queue.pop(0)
            for edge in self._adjacency[current]:
                neighbor = edge.neighbor
                neighbor_tag = edge.tag

                # do not allow for loopbacks
                if neighbor == previous:
                    continue

                neighbor_dist = curr_distance.append(neighbor_tag)
                if is_branching[neighbor]:
                    new_graph.add_edge(
                        new_index[origin], new_index[neighbor], neighbor_dist
                    )
                else:
                    queue.append((neighbor, origin, current, neighbor_dist))
        return new_graph

    def get_component_graph(self):
        """
        Unify all nodes within the same copy of the periodicity cell

        Args:
            self (MetaGraph): The current graph

        """

        periodic_comp = []

        for i in range(self._next_vertex_index):
            periodic_comp.append(-1)

        curr_comp = 0
        no_skip_tag = EdgeTag(0, 0, 0)
        for i in range(self._next_vertex_index):
            if periodic_comp[i] == -1 and self._vertex_states[i] != NodeState.DISABLED:
                # Explore the entire component
                queue = [i]
                while queue:
                    curr = queue.pop(0)
                    # Only process nodes not yet visited
                    if periodic_comp[curr] == -1:
                        periodic_comp[curr] = curr_comp
                        # Look into neighboring vertices
                        for edge in self._adjacency[curr]:
                            if (
                                edge.tag == no_skip_tag
                                and periodic_comp[edge.neighbor] == -1
                                and self._vertex_states[edge.neighbor] != NodeState.DISABLED
                            ):
                                queue.append(edge.neighbor)

                # Iterate to next component
                curr_comp += 1

        new_graph = MetaGraph()
        new_graph.reserve(curr_comp)

        for curr in range(self._next_vertex_index):
            if self._vertex_states[curr] != NodeState.DISABLED:
                base_comp = periodic_comp[curr]
                for edge in self._adjacency[curr]:
                    neighbor = edge.neighbor
                    tag = edge.tag

                    if (
                        tag != no_skip_tag
                        and self._vertex_states[neighbor] != NodeState.DISABLED
                    ):
                        neighbor_comp = periodic_comp[neighbor]
                        new_graph.add_edge(base_comp, neighbor_comp, tag)

        return new_graph

    def find_components(self):
        """
        This simply detects connected components within the graph by doing a flood-fill algorithm.

        Args:
            self (MetaGraph): The current graph

        """
        for i in range(self._next_vertex_index):
            self._vertex_component[i] = -1

        curr_comp = 0
        for i in range(self._next_vertex_index):
            if (
                self._vertex_component[i] == -1
                and self._vertex_states[i] != NodeState.DISABLED
            ):
                # Explore the entire component
                queue = [i]
                while queue:
                    curr = queue.pop(0)
                    # Only process nodes not yet visited
                    if self._vertex_component[curr] == -1:
                        self._vertex_component[curr] = curr_comp
                        # Look into neighboring vertices
                        for edge in self._adjacency[curr]:
                            if self._vertex_component[edge.neighbor] == -1:
                                queue.append(edge.neighbor)

                # Iterate to next component
                curr_comp += 1

        return curr_comp

    def get_components(self):
        """
        Function to retrieve the number of components and the nodes constituting them; together with their metadata

        Args:
            self (MetaGraph): The current graph

        """
        num_components = self.find_components()
        component_data = []

        for i in range(num_components):
            component_data.append([])

        for curr in range(self._next_vertex_index):
            # Check for presence in component
            if self._vertex_component[curr] != -1:
                component_data[self._vertex_component[curr]].append(
                    (curr, self._meta_data[curr])
                )

        return num_components, component_data

    def find_stable_loops(self):
        """
        Function to detect rotationally invariant loops for each individual component.
        It uses a BFS to explore the components for periodicity percolation.
        For the periodicity vectors, it then checks the dimensionality of the basis.
        The return value is the number of components and an array with the dimension of the percolation pattern.
        (0: no percolation, 1: line, 2: sheet, 3: rigid grid)

        Args:
            self (MetaGraph): The current graph

        """
        num_components = self.find_components()

        # Array to keep track, which component contains a loop in this dimension
        component_percolation_dimension = []
        # Bitlist, which component has been explored already
        component_explored = []
        for i in range(num_components):
            component_percolation_dimension.append(0)
            component_explored.append(False)

        # Array to keep track, which node has already been visited
        visited = []
        # Array holding the required dimension-crossing distance for a node
        distance = []

        # Filter out disabled nodes from visitation list
        for i in range(self._next_vertex_index):
            distance.append(-1)
            if self._vertex_states[i] == NodeState.DISABLED:
                visited.append(True)
            else:
                visited.append(False)

        # Explore graph
        for i in range(self._next_vertex_index):
            if (
                visited[i]
                or self._vertex_component[i] == -1
                or component_explored[self._vertex_component[i]]
            ):
                continue

            # Mark component as being explored
            curr_component = self._vertex_component[i]
            component_explored[curr_component] = True

            # Basis of percolation vectors
            basis = []

            # Do a BFS on the component. If we encounter two different distances for a node, we have a loop
            start_tag = EdgeTag(0, 0, 0)
            queue = [(i, start_tag)]
            distance[i] = start_tag

            while queue:
                curr, curr_dist = queue.pop(0)

                if visited[curr]:
                    # Found two different crossing numbers to this node, Loop detected
                    if distance[curr] != curr_dist:
                        basis = make_basis(
                            basis, curr_dist.calculate_connecting(
                                distance[curr])
                        )
                        if len(basis) >= 3:
                            queue.clear()
                            break
                    # Otherwise just skip this one
                else:
                    visited[curr] = True
                    distance[curr] = curr_dist
                    # print(curr_dist)

                    for edge in self._adjacency[curr]:
                        next_node = edge.neighbor
                        # Calculate next distance
                        next_dist = curr_dist.append(edge.tag)

                        if visited[next_node]:
                            # Found a different length path => loop
                            if distance[next_node] != next_dist:
                                dist = next_dist.calculate_connecting(
                                    distance[next_node])
                                basis = make_basis(
                                    basis,
                                    dist
                                )
                                if len(basis) >= 3:
                                    queue.clear()
                                    break
                            else:
                                # Don't bother finding another path of the same length
                                continue
                        else:
                            queue.append((next_node, next_dist))
            component_percolation_dimension[curr_component] = len(basis)

        return num_components, component_percolation_dimension

    def copy(self):
        """
        Create a copy of this graph with no shared data

        Args:
            self (MetaGraph): The current graph

        """
        new_graph = MetaGraph()
        new_graph._adjacency = deepcopy(self._adjacency)
        new_graph._vertex_states = deepcopy(self._vertex_states)
        new_graph._vertex_component = deepcopy(self._vertex_component)
        new_graph._meta_data = deepcopy(self._meta_data)
        new_graph._next_vertex_index = self._next_vertex_index
        new_graph._next_edge_index = self._next_edge_index
        return new_graph
