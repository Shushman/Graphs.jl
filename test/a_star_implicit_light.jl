# Similar to test for standard A*
using Graphs
using Test

# The "touring Romania" example from Russell and Norvig
g1_wedges = [
    (1, 20, 75.),
    (1, 16, 140.), #2, Arad -- Sibiu
    (1, 17, 118.),
    (20, 13, 71.),
    (13, 16, 151.),
    (16, 6, 99.),
    (16, 15, 80.), #7, Sibiu -- Rimnicu Vilcea
    (6, 2, 211.),
    (17, 10, 111.),
    (10, 11, 70.),
    (11, 4, 75.),
    (4, 3, 120.),
    (3, 15, 146.),
    (3, 14, 138.),
    (15, 14, 97.), #15, Rimnicu Vilcea -- Pitesti
    (14, 2, 101.), #16, Pitesti -- Bucharest
    (2, 7, 90.),
    (2, 18, 85.),
    (18, 8, 98.),
    (8, 5, 86.),
    (18, 19, 142.),
    (19, 9, 92.),
    (9, 12, 87.) ]

g1_heuristics = [
    366, 0, 160, 242, 161, 176, 77, 151, 226, 244, 241, 234, 380, 100, 193,
    253, 329, 80, 199, 374 ]

# Create a dictionary to map indices to neighbours and edge weights
nbrs_wt_dict = Dict{Int64, Dict{Int64, Float64}}()

for (u, v, wt) in g1_wedges
    if haskey(nbrs_wt_dict, u)
        nbrs_wt_dict[u][v] = wt
    else
        nbrs_wt_dict[u] = Dict(v => wt)
    end
end

# Setup edge weight function based on dictionary
edge_wt_fn(u::Int64, v::Int64) = nbrs_wt_dict[u][v]

struct TestVisitorImplicit <: AbstractDijkstraVisitor
    nbrs_wt_dict::Dict{Int64, Dict{Int64, Float64}}
    test_graph::SimpleVListGraph
    goal_vtx::Int64
end

function Graphs.include_vertex!(vis::TestVisitorImplicit, u::Int64, v::Int64,
                                d::Float64, nbrs::Vector{Int64})

    # Terminate if goal vertex reached
    if v == vis.goal_vtx
        return false
    end

    # Otherwise, insert all neighbours from dict
    # AND add to graph if not already there
    for n in collect(keys(vis.nbrs_wt_dict[v]))
        if ~(n in vis.test_graph.vertices)
            add_vertex!(test_graph, n)
        end

        # Neighbors will have INDICES of vertices, not vertices themselves
        push!(nbrs, vertex_index(vis.test_graph, n))
    end

    return nbrs
end

# Create initial graph with only start
start = 1
goal = 2
test_graph = SimpleVListGraph([start])

# Create visitor with goal vertex 2 and heuristic
vis = TestVisitorImplicit(nbrs_wt_dict, test_graph, goal)
heur(n) = g1_heuristics[n]

a_star_sp_states = a_star_light_shortest_path_implicit!(test_graph, edge_wt_fn, 1, vis, heur)

# Now actually extract the shortest path indices
sp_idxs = shortest_path_indices(a_star_sp_states, test_graph, start, goal)
sp = [test_graph.vertices[v] for v in sp_idxs]

@test sp == [1, 16, 15, 14, 2]
