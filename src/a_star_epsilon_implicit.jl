# TODO: Can just reuse AStarHEntry??
# struct AStarHEntry{D <: Number}
#     vIdx::Int
#     gvalue::D
#     fvalue::D
# end
#
# Base.isless(e1::AStarHEntry, e2::AStarHEntry) = (e1.fvalue < e2.fvalue) || (e1.fvalue == e2.fvalue && e1.gvalue < e2.gvalue)

import DataStructure: compare

struct AStarEpsilonHEntry{D <: Number}
    vIdx::Int
    gvalue::D
    fvalue::D
    focal_heuristic::D
end

struct CompareHeap
end
compare(comp::CompareHeap, e1::AStarEpsilonHEntry, e2::AStarEpsilonHEntry) = (e1.fvalue, e1.gvalue) < (e2.fvalue, e2.gvalue)

struct CompareFocalHeap
end
compare(comp::CompareFocalHeap, e1::AStarEpsilonHEntry, e2::AStarEpsilonHEntry) = (e1.focal_heuristic, e1.fvalue, -e1.gvalue) < (e2.focal_heuristic, e2.fvalue, -e2.gvalue)



"""
Compared with a_star_spath, the attributes of AStarStates are all Dicts instead of
fixed-length Vectors, as arbitrarily many vertices can be added by the visitors.
"""
@with_kw mutable struct AStarEpsilonStates{D<:Number}
    parent_indices::Dict{Int,Int} = Dict{Int,Int}()
    dists::Dict{Int,D} = Dict{Int,D}()
    colormap::Dict{Int,Int} = Dict{Int,Int}()
    heap::MutableBinaryHeap{AStarEpsilonHEntry{D},CompareHeap} = MutableBinaryHeap{AStarEpsilonHEntry{D},CompareHeap}()
    hmap::Dict{Int,Int} = Dict{Int,Int}()
    # Subset of heap
    focal_heap::MutableBinaryHeap{AStarEpsilonHEntry{D},CompareFocalHeap} = MutableBinaryHeap{AStarEpsilonHEntry{D},CompareFocalHeap}()
    focal_hmap::Dict{Int,Int} = Dict{Int,Int}()
    best_fvalue::Float64 = 0.0
end


function set_source!(state::AStarStates{D}, s::Int) where {D <: Number, V}
    state.parent_indices[s] = s
    state.dists[s] = zero(D)
    state.colormap[s] = 2
end

"""
Execute expand operation on the chosen node, while implicitly generating its neighbors based on a
visitor method, and populating `neighbors` with the neighboring vertices.
"""
function process_neighbors_implicit!(
    state::AStarStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    neighbors::Vector{Int},
    parent_entry::AStarEpsilonHEntry{D},
    visitor::AbstractDijkstraVisitor,
    weight::Float64,
    admissible_heuristic::Function,
    focal_state_heuristic::Function,
    focal_transition_heuristic::Function) where {V, D <: Number}

    dv = zero(D)
    u = parent_entry.vIdx
    du = parent_entry.gvalue

    for iv in neighbors

        # Default color 0
        v_color = get(state.colormap, iv, 0)

        if v_color == 0

            # Inserting for the first time
            state.dists[iv] = dv = du + edge_wt_fn(graph.vertices[u], graph.vertices[iv])
            state.parent_indices[iv] = u
            state.colormap[iv] = 1
            Graphs.discover_vertex!(visitor, graph.vertices[u], graph.vertices[iv], dv)

            new_fvalue = dv + admissible_heuristic(graph.vertices[iv])
            new_focal_heur = parent_entry.focal_heuristic + focal_state_heuristic(graph.vertices[iv]) +
                                                            focal_transition_heuristic(graph.vertices[u], graph.vertices[iv])

            new_entry = AStarEpsilonHEntry(iv, dv, new_fvalue, new_focal_heur))

            state.hmap[iv] = push!(state.heap, new_entry)

            # Only insert in focal list if condition satisfied
            if new_fvalue <= weight * state.best_fvalue
                state.focal_hmap[iv] = push!(state.focal_heap, new_entry)
            end

        elseif v_color == 1
            # Must also

            dv = du + edge_wt_fn(graph.vertices[u], graph.vertices[iv])

            # Update cost-to-come if cheaper from current parent
            if dv < state.dists[iv]

                state.dists[iv] = dv
                state.parent_indices[iv] = u

                Graphs.update_vertex!(visitor, graph.vertices[u], graph.vertices[iv], dv)

                updated_fvalue = dv + admissible_heuristic(graph.vertices[iv])
                updated_entry = AStarEpsilonHEntry(iv, dv, updated_fvalue, parent_entry.focal_heuristic)

                update!(state.heap, state.hmap[iv], updated_entry)

                # Enter into focal list if new fvalue is good enough
                # but it was not before
                if updated_fvalue <= weight * state.best_fvalue &&
                    parent_entry.fvalue > weight * state.best_fvalue
                    state.focal_hmap[iv] = push!(state.focal_heap, updated_entry)
                end

            end
        end
    end
end


function a_star_light_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::Function, # distances associated with edges
    source::Int,             # the source index
    visitor::AbstractDijkstraVisitor,# visitor object\
    weight::Float64,
    admissible_heuristic::Function,      # Heuristic function for vertices
    focal_state_heuristic::Function,
    focal_transition_heuristic::Function,
    state::AStarStates{D}) where {V, D <: Number}

    d0 = zero(D)
    set_source!(state, source)

    # Will be populated by include_vertex!
    source_nbrs = Vector{Int}(undef,0)

    # Call the user-defined visitor function for expanding the source

    if Graphs.include_vertex!(visitor, graph.vertices[source], graph.vertices[source], d0, source_nbrs) == false
        return state
    end

    # Create dummy source entry
    source_entry = AStarEpsilonHEntry(source, d0, d0 + admissible_heuristic(graph.vertices[source]), d0)
    process_neighbors_implicit!(state, graph, edge_wt_fn, source_nbrs, source_entry, visitor,
                                weight, admissible_heuristic, focal_state_heuristic, focal_transition_heuristic)
    Graphs.close_vertex!(visitor, graph.vertices[source])

    while ~(isempty(state.heap))

        # Enter open-set items that are now valid into focal list
        old_best_fvalue = state.best_fvalue
        best_fvalue = top(state.heap).fvalue

        if best_fvalue > old_best_fvalue

            # Iterate over open set  in increasing order of fvalue and insert in focal list if valid
            for node in sort(state.heap.nodes, by = x->x.value.fvalue)
                fvalue = node.value.fvalue
                if fvalue > old_best_fvalue && fvalue <= weight * best_fvalue
                    state.focal_hmap[node.value.vIdx] = push!(state.focal_heap, node.value)
                end
                if fvalue > weight * old_best_fvalue
                    break
                end
            end
        end

        # pick next vertex to include
        # TODO: Pop focal heap
        focal_entry, focal_handle = top_with_handle(state.focal_heap)
        heap_handle = state.hmap[focal_entry.vIdx]

        ui = focalentry.vIdx
        du = focalentry.gvalue
        state.colormap[ui] = 2

        # Will be populated by include_vertex!
        nbrs = Vector{Int}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[state.parent_indices[ui]], graph.vertices[ui], du, nbrs) == false
            return state
        end

        # Delete from open list before considering neighbors
        # TODO: Check this!!!
        pop!(state.focal_heap)
        delete!(state.heap, heap_handle)

        # process u's neighbors
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs, focal_entry, visitor,
                                    weight, admissible_heuristic, focal_state_heuristic, focal_transition_heuristic)
        Graphs.close_vertex!(visitor, graph.vertices[ui])
    end

    state
end

function a_star_light_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::Function, # distances associated with edges
    source::Int,
    visitor::AbstractDijkstraVisitor,
    weight::Float64,
    admissible_heuristic::Function = n -> 0,
    focal_state_heuristic::Function = n -> 0,
    focal_transition_heuristic::Function = (u, v) -> 0,
    ::Type{D} = Float64) where {V, D <: Number}
    state = AStarEpsilonStates{D}()
    a_star_light_shortest_path_implicit!(graph, edge_wt_fn, source, visitor, heuristic, state)
end

"""
Given the AStarStates result, extract the shortest path
"""
function shortest_path_indices(state::AStarStates{D}, graph::AbstractGraph{V},
                       source::V, target::V) where {D, V}

    source_idx = vertex_index(graph, source)
    target_idx = vertex_index(graph, target)

    @assert haskey(state.parent_indices, target_idx) == true "Target has no parent!"

    sp = [target_idx]

    # Walk back from target to source using parent indices
    curr_idx = target_idx
    while curr_idx != source_idx
        curr_idx = state.parent_indices[curr_idx]
        pushfirst!(sp, curr_idx)
    end

    return sp
end
