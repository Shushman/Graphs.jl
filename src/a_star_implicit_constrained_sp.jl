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

# A path is dominated if cost and weight are both worse than an existing entry
function is_dominated(state::AStarMCSPStates{D}, s::AStarMCSPHEntry{D}) where {D <: Number}

    if haskey(state.open_list_hmap, s.v_idx)
        for ol_idx in state.open_list_hmap[s.v_idx]
            ols_entry = state.open_list.nodes[state.open_list.node_map[ol_idx]].value
            if ols_entry.gvalue <= s.gvalue && old_entry.weights <= s.weights
                return true
            end
        end
    end

    if haskey(state.closed_list_hmap, s.v_idx)
        for cl_idx in state.closed_list_hmap[s.v_idx]
            cls_entry = state.closed_list.nodes[state.closed_list.node_map[cl_idx]].value
            if cls_entry.gvalue <= s.gvalue && cls_entry.weights <= s.weights
                return true
            end
        end
    end

    return false
end

# A new entry is eligble if the weight-to-go based on this transition
# plus the heuristic weight-to-go is within the weight constraint
function is_eligible(graph::AbstractGraph{V},
                     entry::AStarMCSPHEntry{D},
                     nbr::Int64,
                     weight_functions::Vector{F1},
                     weight_heuristics::Vector{F2},
                     weight_constraints::Vector{D}) where {V, D <: Number, F1 <: Function, F2 <: Function}

    n_weights = length(weight_functions)
    new_weight_vector = zeros(D, n_weights)

    for i = 1:n_weights
        new_constr_wt = entry.weights[i] +  weight_functions[i](graph.vertices[entry.v_idx], graph.vertices[nbr])

        if new_constr_wt + weight_heuristics[i](graph.vertices[nbr])  > weight_constraints[i]
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
    edge_wt_fn::F1,
    neighbors::Vector{Int64},
    parent_entry::AStarMCSPHEntry{D},
    entry_cl_handle::Int64,    # For closed list
    visitor::AbstractDijkstraVisitor,
    heuristic::F2,
    weight_functions::Vector{F3},
    weight_heuristics::Vector{F4},
    weight_constraints::Vector{D}
    ) where {V, D <: Number, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    for nbr_idx in neighbors


        is_elig, new_weight_vector =  is_eligible(graph, parent_entry, nbr_idx, weight_functions, weight_heuristics, weight_constraints)

        if is_elig
            # Create new entry
            new_gvalue = parent_entry.gvalue + edge_wt_fn(graph.vertices[parent_entry.v_idx], graph.vertices[nbr_idx])
            new_fvalue = new_gvalue + heuristic(graph.vertices[nbr_idx])
            new_entry = AStarMCSPHEntry(nbr_idx, parent_entry.v_idx, new_gvalue, new_fvalue,
                                        new_weight_vector, entry_cl_handle)

            if ~(is_dominated(state, new_entry))
                new_entry_handle = push!(state.open_list, new_entry)
                if ~(haskey(state.open_list_hmap, nbr_idx))
                    Graphs.discover_vertex!(visitor, graph.vertices[parent_entry.v_idx], graph.vertices[nbr_idx], new_gvalue)
                    state.open_list_hmap[nbr_idx] = Set{Int64}(new_entry_handle)
                else
                    push!(state.open_list_hmap[nbr_idx], new_entry_handle)
                end
            end

        end
    end
end



function a_star_constrained_shortest_path_implicit!(
    graph::AbstractGraph{V},
    edge_wt_fn::F1,
    source::Int64,
    visitor::AbstractDijkstraVisitor,
    heuristic::F2,
    weight_functions::Vector{F3},
    weight_heuristics::Vector{F4},
    weight_constraints::Vector{D},
    ::Type{D} = Float64) where {V, D <: Number, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function}

    state = AStarMCSPStates{D}()

    n_weights = length(weight_constraints)

    source_entry = AStarMCSPHEntry(source, source, zero(D), heuristic(source), zeros(D, n_weights), 0)
    source_handle = push!(state.open_list, source_entry)
    state.open_list_hmap[source] = Set{Int64}(source_handle)

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
        if ~(haskey(state.closed_list_hmap, entry.v_idx))
            state.closed_list_hmap[entry.v_idx] = Set{Int64}(cl_entry_handle)
        else
            push!(state.closed_list_hmap[entry.v_idx], cl_entry_handle)
        end

        # Process nbrs
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs,
                                    entry, cl_entry_handle, visitor,
                                    heuristic, weight_functions, weight_heuristics, weight_constraints)
        Graphs.close_vertex!(visitor, graph.vertices[entry.v_idx])
    end

    return state, source_entry
end

# Given the solution state for the A*MCSP Heap, walk back using the closed list
# to actually get the shortest feasible path
function shortest_path_cost_weights(state::AStarMCSPStates{D}, graph::AbstractGraph{V},
                                    source::Int64, target_entry::AStarMCSPHEntry{D}) where {V, D <: Number}

    target_idx = target_entry.v_idx

    sp = [target_idx]

    # Walk back USING closed list heap
    curr_idx = target_idx
    curr_par_handle = target_entry.parent_handle

    while curr_idx != source

        cl_list_entry = state.closed_list.nodes[state.closed_list.node_map[curr_par_handle]].value

        curr_idx = cl_list_entry.v_idx
        curr_par_handle = cl_list_entry.parent_handle
        pushfirst!(sp, curr_idx)
    end

    return sp, target_entry.gvalue, target_entry.weights
end
