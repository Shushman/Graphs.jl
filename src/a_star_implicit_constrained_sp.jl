## A-star implicit with one or more additive weight constraints on feasible paths
## Basically, find the shortest path (w.r.t cost) that satisfies a bunch of
## constraints that depend on other weight functions
struct AStarMCSPHEntry{D <: Number}
    v_idx::Int64
    parent_idx::Int64
    gvalue::D
    fvalue::D
    weights::Vector{D}
    parent_handle::Int64    # refers to the heap handle in CLOSED LIST that tracks this path
end

Base.isless(e1::AStarMCSPHEntry, e2::AStarMCSPHEntry) = (e1.fvalue, e1.gvalue) < (e2.fvalue, e2.gvalue) # Lexicographic

@with_kw mutable struct AStarMCSPStates{D <: Number}
    open_list::MutableBinaryMinHeap{AStarMCSPHEntry{D}}    = MutableBinaryMinHeap{AStarMCSPHEntry{D}}()
    closed_list::MutableBinaryMinHeap{AStarMCSPHEntry{D}}  = MutableBinaryMinHeap{AStarMCSPHEntry{D}}()
    open_list_hmap::Dict{Int64,Set{Int64}}                 = Dict{Int64,Set{Int64}}()
    closed_list_hmap::Dict{Int64,Set{Int64}}               = Dict{Int64,Set{Int64}}()
end

function is_dominated(state::AStarMCSPStates{D}, s::AStarMCSPHEntry{D}) where {D <: Number}

    for ol_idx in state.open_list_hmap[s.v_idx]
        ols_entry = state.open_list.nodes[state.open_list.node_map[ol_idx]]
        if isless(ols_entry, s)
            return true
        end
    end

    for cl_idx in state.closed_list_hmap[s.v_idx]
        cls_entry = state.closed_list.nodes[state.closed_list.node_map[ol_idx]]
        if isless(cls_entry, s)
            return true
        end
    end

    return false
end


function is_eligible(entry::AStarMCSPHEntry{D},
                     nbr::Int64,
                     weight_functions::Vector{Function},
                     weight_heuristics::Vector{Function},
                     weight_constraints::Vector{D}) where {D <: Number}

    n_weights = length(weight_functions)
    new_weight_vector = zeros(D, n_weights)

    for i = 1:n_weights
        new_constr_wt = entry.weights[i] +  weight_functions[i](entry.v_idx, nbr) +
                        weight_heuristics[i](nbr)
        if new_constr_wt  > weight_constraints[i]
            return false, new_weight_vector
        else
            new_weight_vector[i] = new_constr_wt
        end
    end

    return true, new_weight_vector
end

function process_neighbors_implicit!(
    state::AStarMCSPStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    neighbors::Vector{Int64},
    entry::AStarMCSPHEntry{D},
    entry_cl_handle::Int64,    # For closed list
    visitor::AbstractDijkstraVisitor,
    heuristic::Function,
    weight_functions::Vector{Function},
    weight_heuristics::Vector{Function},
    weight_constraints::Vector{D}
    )

    for nbr_idx in neighbors

        is_elig, new_wt_vector =  is_eligible(entry, nbr_idx, weight_functions, weight_heuristics, weight_constraints)

        if is_elig

            # Create new entry
            new_gvalue = entry.gvalue + edge_wt_fn(entry.v_idx, nbr_idx)
            new_fvalue = new_gvalue + heuristic(nbr_idx)
            new_entry = AStarMCSPHEntry(nbr_idx, entry.v_idx, new_gvalue, new_fvalue,
                                        new_weight_vector, entry_cl_handle)

            if ~(is_dominated(state, new_entry))
                new_entry_handle = push!(state.open_list, new_entry)
                push!(state.open_list_hmap, new_entry_handle)
            end

        end
    end
end



function a_star_constrained_shortest_path_implicit!(
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    source::Int64,
    visitor::AbstractDijkstraVisitor,
    heuristic::Function,
    weight_heuristics::Vector{Function},
    weight_constraints::Vector{D},
    ::Type{D} = Float64) where {V, D <: Number}

    state = AStarMCSPStates{D}()

    n_weights = length(weight_constraints)

    source_entry = AStarMCSPHEntry(source, zero(D), heuristic(source), zeros(D, n_weights), 0, 0)
    push!(state.open_list, source_entry)

    while ~(isempty(state.open_list))

        # Remove from open_list and open_list_hmap
        entry, entry_handle = top_with_handle(state.open_list)
        pop!(state.open_list)
        delete!(state.open_list_hmap[entry.v_idx], entry_handle)

        nbrs = Vector{Int64}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[entry.parent_idx],
                                  graph.vertices[entry.v_idx], entry.gvalue, nbrs) == false
            return state, entry
        end

        # Insert into closed list to get new handle
        cl_entry_handle = push!(state.closed_list, entry)
        push!(state.closed_list_hmap[entry.v_idx], cl_entry_handle)

        # Process nbrs
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs,
                                    entry, cl_entry_handle, visitor,
                                    heuristic, weight_functions, weight_heuristics, weight_constraints)
        Graphs.close_vertex!(visitor, graph.vertices[entry.v_idx])
    end

    return state, source_entry
end

function shortest_path_cost_weights(state::AStarMCSPStates{D}, graph::AbstractGraph{V},
                                    source::Int64, target_entry::AStarMCSPHEntry{D}) where {D <: Number}

    target_idx = target_entry.v_idx

    sp = [target_idx]

    # Walk back USING closed list heap
    curr_idx = target_idx
    curr_par_handle = target_entry.parent_handle

    while curr_idx != source

        cl_list_entry = state.closed_list.nodes[state.closed_list.node_map[curr_par_handle]]

        curr_idx = cl_list_entry.v_idx
        curr_par_handle = cl_list_entry.parent_handle
        pushfirst!(sp, curr_idx)
    end

    return sp, target_entry.gvalue, target_entry.weights
end
