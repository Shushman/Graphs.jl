# A* shortest-path search with two generalizations over a_star_spath:
# 1. The graph is implicit; edges are not explicitly stored but instead generated
# during vertex expansion by a user-defined visitor function
# 2. The include_vertex! visitor function can also ADD new vertices, i.e. the size
# of the graph can increase during the search
#
# We do not bother with deleting vertices as it is typically not worth the effort
# for implicit search

"""
Data structure for heap entry associated with each expanded vertex.

Attributes:
    - `v_idx::Int64` The integer index of the vertex being expanded
    - `gvalue::D` The cost-to-come of the vertex from the source
    - `fvalue::D` The cost-to-come + heuristic cost-to-go to the goal
"""
struct AStarHEntry{D <: Number}
    v_idx::Int64
    gvalue::D
    fvalue::D
end

Base.isless(e1::AStarHEntry, e2::AStarHEntry) = (e1.fvalue, e1.gvalue) < (e2.fvalue, e2.gvalue) # Lexicographic ordering


"""
Compared with a_star_spath, the attributes of AStarStates are all Dicts instead of
fixed-length Vectors, as arbitrarily many vertices can be added by the visitors.
"""
@with_kw mutable struct AStarStates{D<:Number}
    parent_indices::Dict{Int64,Int64} = Dict{Int64,Int64}()
    dists::Dict{Int64,D} = Dict{Int64,D}()
    colormap::Dict{Int64,Int64} = Dict{Int64,Int64}()
    heap::MutableBinaryMinHeap{AStarHEntry{D}} = MutableBinaryMinHeap{AStarHEntry{D}}()
    hmap::Dict{Int64,Int64} = Dict{Int64,Int64}()
end


function set_source!(state::AStarStates{D}, s::Int64) where {D <: Number, V}
    state.parent_indices[s] = s
    state.dists[s] = zero(D)
end

"""
Execute expand operation on the chosen node, while implicitly generating its neighbors based on a
visitor method, and populating `neighbors` with the neighboring vertices.
"""
function process_neighbors_implicit!(
    state::AStarStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    neighbors::Vector{Int64},
    u::Int64, du::D, visitor::AbstractDijkstraVisitor,
    heuristic::Function) where {V, D <: Number}

    dv = zero(D)

    for iv in neighbors

        # Default color 0
        v_color = get(state.colormap, iv, 0)

        if v_color == 0

            # Inserting for the first time
            state.dists[iv] = dv = du + edge_wt_fn(graph.vertices[u], graph.vertices[iv])
            state.parent_indices[iv] = u
            state.colormap[iv] = 1
            Graphs.discover_vertex!(visitor, graph.vertices[u], graph.vertices[iv], dv)
            state.hmap[iv] = push!(state.heap, AStarHEntry(iv, dv, dv + heuristic(graph.vertices[iv])))

        elseif v_color == 1

            dv = du + edge_wt_fn(graph.vertices[u], graph.vertices[iv])

            # Update cost-to-come if cheaper from current parent
            if dv < state.dists[iv]

                state.dists[iv] = dv
                state.parent_indices[iv] = u

                Graphs.update_vertex!(visitor, graph.vertices[u], graph.vertices[iv], dv)
                update!(state.heap, state.hmap[iv], AStarHEntry(iv, dv, dv + heuristic(graph.vertices[iv])))
            end
        end
    end
end


function a_star_light_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::Function, # distances associated with edges
    source::Int64,             # the source index
    visitor::AbstractDijkstraVisitor,# visitor object
    heuristic::Function,      # Heuristic function for vertices
    ::Type{D} = Float64) where {V, D <: Number}

    state = AStarStates{D}()
    set_source!(state, source)
    source_entry = AStarHEntry(source, zero(D), heuristic(graph.vertices[source]))
    state.hmap[source] = push!(state.heap, source_entry)

    while ~(isempty(state.heap))

        # pick next vertex to include
        entry = pop!(state.heap)
        ui = entry.v_idx
        du = entry.gvalue

        state.colormap[ui] = 2

        # Will be populated by include_vertex!
        nbrs = Vector{Int64}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[state.parent_indices[ui]], graph.vertices[ui], du, nbrs) == false
            return state
        end

        # process u's neighbors
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs, ui, du, visitor, heuristic)
        Graphs.close_vertex!(visitor, graph.vertices[ui])
    end

    state
end

"""
Given the AStarStates result, extract the shortest path
"""
function shortest_path_indices(parent_indices::Dict{Int64, Int64}, graph::AbstractGraph{V},
                               source_idx::Int64, target_idx::Int64) where {V}

    @assert haskey(parent_indices, target_idx) == true "Target has no parent!"

    sp = [target_idx]

    # Walk back from target to source using parent indices
    curr_idx = target_idx
    while curr_idx != source_idx
        curr_idx = parent_indices[curr_idx]
        pushfirst!(sp, curr_idx)
    end

    return sp
end
