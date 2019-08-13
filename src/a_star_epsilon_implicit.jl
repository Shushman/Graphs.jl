struct AStarEpsilonHEntry{D <: Number}
    v_idx::Int64
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
The hmaps map the vertex index (which is unique) to the heap handle; and can be used to cross-reference entries
across heaps. This is our state-to-heap
"""
@with_kw mutable struct AStarEpsilonStates{D<:Number}
    parent_indices::Dict{Int64,Int64} = Dict{Int64,Int64}()
    dists::Dict{Int64,D} = Dict{Int64,D}()
    colormap::Dict{Int64,Int64} = Dict{Int64,Int64}()
    heap::MutableBinaryHeap{AStarEpsilonHEntry{D},CompareHeap} = MutableBinaryHeap{AStarEpsilonHEntry{D},CompareHeap}()
    hmap::Dict{Int64,Int64} = Dict{Int64,Int64}()
    # Subset of heap
    focal_heap::MutableBinaryHeap{AStarEpsilonHEntry{D},CompareFocalHeap} = MutableBinaryHeap{AStarEpsilonHEntry{D},CompareFocalHeap}()
    focal_hmap::Dict{Int64,Int64} = Dict{Int64,Int64}()
    best_fvalue::D = zero(D)
end


function set_source!(state::AStarEpsilonStates{D}, s::Int64) where {D <: Number, V}
    state.parent_indices[s] = s
    state.dists[s] = zero(D)
    state.colormap[s] = 2
end

"""
Execute expand operation on the chosen node, while implicitly generating its neighbors based on a
visitor method, and populating `neighbors` with the neighboring vertices.
"""
function process_neighbors_implicit!(
    state::AStarEpsilonStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    neighbors::Vector{Int64},
    parent_entry::AStarEpsilonHEntry{D},
    visitor::AbstractDijkstraVisitor,
    eps_weight::Float64,
    admissible_heuristic::Function,
    focal_state_heuristic::Function,
    focal_transition_heuristic::Function) where {V, D <: Number}

    dv = zero(D)
    u = parent_entry.v_idx
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

            new_entry = AStarEpsilonHEntry(iv, dv, new_fvalue, new_focal_heur)

            # @show new_entry

            state.hmap[iv] = push!(state.heap, new_entry)

            # Only insert in focal list if condition satisfied
            if new_fvalue <= eps_weight * state.best_fvalue
                state.focal_hmap[iv] = push!(state.focal_heap, new_entry)
            end

        elseif v_color == 1

            dv = du + edge_wt_fn(graph.vertices[u], graph.vertices[iv])

            # Update cost-to-come if cheaper from current parent
            if dv < state.dists[iv]

                state.dists[iv] = dv
                state.parent_indices[iv] = u

                Graphs.update_vertex!(visitor, graph.vertices[u], graph.vertices[iv], dv)

                updated_fvalue = dv + admissible_heuristic(graph.vertices[iv])
                updated_entry = AStarEpsilonHEntry(iv, dv, updated_fvalue, parent_entry.focal_heuristic)

                old_fvalue = state.heap.nodes[state.heap.node_map[state.hmap[iv]]].value.fvalue
                # @show old_fvalue

                # @show updated_entry
                update!(state.heap, state.hmap[iv], updated_entry)

                # Enter into focal list if new fvalue is good enough
                # but it was not before
                if updated_fvalue <= eps_weight * state.best_fvalue &&
                    old_fvalue > eps_weight * state.best_fvalue
                    state.focal_hmap[iv] = push!(state.focal_heap, updated_entry)
                end

            end
        end
    end
end


function a_star_light_epsilon_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::Function, # distances associated with edges
    source::Int64,             # the source index
    visitor::AbstractDijkstraVisitor,# visitor object\
    eps_weight::Float64,
    admissible_heuristic::Function,      # Heuristic function for vertices
    focal_state_heuristic::Function,
    focal_transition_heuristic::Function,
    state::AStarEpsilonStates{D}) where {V, D <: Number}

    d0 = zero(D)
    set_source!(state, source)
    source_heur = admissible_heuristic(graph.vertices[source])
    state.best_fvalue = source_heur

    # Will be populated by include_vertex!
    source_nbrs = Vector{Int64}(undef, 0)

    # Call the user-defined visitor function for expanding the source
    if Graphs.include_vertex!(visitor, graph.vertices[source], graph.vertices[source], d0, source_nbrs) == false
        return state
    end

    # Create dummy source entry
    source_entry = AStarEpsilonHEntry(source, d0, source_heur, d0)
    process_neighbors_implicit!(state, graph, edge_wt_fn, source_nbrs, source_entry, visitor,
                                eps_weight, admissible_heuristic, focal_state_heuristic, focal_transition_heuristic)
    Graphs.close_vertex!(visitor, graph.vertices[source])

    while ~(isempty(state.heap))

        # Enter open-set items that are now valid into focal list
        old_best_fvalue = state.best_fvalue
        state.best_fvalue = top(state.heap).fvalue

        # @show state.heap
        # @show state.hmap
        # @show state.focal_heap
        # @show state.focal_hmap

        if state.best_fvalue > old_best_fvalue

            # Iterate over open set  in increasing order of fvalue and insert in focal list if valid
            for node in sort(state.heap.nodes, by = x->x.value.fvalue)
                fvalue = node.value.fvalue

                if fvalue > eps_weight * old_best_fvalue && fvalue <= eps_weight * state.best_fvalue
                    state.focal_hmap[node.value.v_idx] = push!(state.focal_heap, node.value)
                end

                if fvalue > eps_weight * state.best_fvalue
                    break
                end
            end
        end

        # Temporary - check consistency
        # TODO: Remove later

        mismatch = false
        best_fvalue = top(state.heap).fvalue
        for node in sort(state.heap.nodes, by = x->x.value.fvalue)
            fvalue = node.value.fvalue

            if fvalue <= eps_weight * best_fvalue
                # Check for entry in focal heap
                if ~(haskey(state.focal_hmap, node.value.v_idx))
                    @info "Focal set missing $(node.value) for f-val $(best_fvalue)"
                    mismatch = true
                end
            else
                if haskey(state.focal_hmap, node.value.v_idx)
                    @info "Focal set should not have ",node.value.v_idx
                end
            end
        end
        @assert mismatch == false


        # pick next vertex to include
        focal_entry, focal_handle = top_with_handle(state.focal_heap)
        # @show graph.vertices[focal_entry.v_idx]
        # @show focal_entry
        heap_handle = state.hmap[focal_entry.v_idx]

        ui = focal_entry.v_idx
        du = focal_entry.gvalue
        state.colormap[ui] = 2

        # Will be populated by include_vertex!
        nbrs = Vector{Int64}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[state.parent_indices[ui]], graph.vertices[ui], du, nbrs) == false
            return state
        end
        # @show [graph.vertices[n] for n in nbrs]

        # Delete from open list before considering neighbors
        # TODO: Check this!!!
        pop!(state.focal_heap)
        # delete!(state.focal_heap, focal_handle)
        delete!(state.focal_hmap, focal_entry.v_idx)
        delete!(state.heap, heap_handle)
        delete!(state.hmap, focal_entry.v_idx)

        # process u's neighbors
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs, focal_entry, visitor,
                                    eps_weight, admissible_heuristic, focal_state_heuristic, focal_transition_heuristic)
        Graphs.close_vertex!(visitor, graph.vertices[ui])
    end

    state
end

function a_star_light_epsilon_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::Function, # distances associated with edges
    source::Int64,
    visitor::AbstractDijkstraVisitor,
    eps_weight::Float64,
    admissible_heuristic::Function = s -> 0,
    focal_state_heuristic::Function = (s, gs) -> 0,
    focal_transition_heuristic::Function = (s1, s2, gs1, gs2) -> 0,
    ::Type{D} = Float64) where {V, D <: Number}
    state = AStarEpsilonStates{D}()
    a_star_light_epsilon_shortest_path_implicit!(graph, edge_wt_fn, source, visitor, eps_weight,
                                         admissible_heuristic, focal_state_heuristic, focal_transition_heuristic, state)
end
