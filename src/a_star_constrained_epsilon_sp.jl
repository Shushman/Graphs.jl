struct AStarEpsilonMCSPHEntry{D <: Number}
    v_idx::Int64
    parent_idx::Int64
    gvalue::D
    fvalue::D
    focal_heuristic::D
    weights::Vector{D}
    parent_handle::Int64
end

Base.isless(e1::AStarEpsilonMCSPHEntry, e2::AStarEpsilonMCSPHEntry) = (e1.fvalue, e1.gvalue) < (e2.fvalue, e2.gvalue)

compare(comp::CompareFocalHeap, e1::AStarEpsilonMCSPHEntry, e2::AStarEpsilonMCSPHEntry) = (e1.focal_heuristic, e1.fvalue, -e1.gvalue) < (e2.focal_heuristic, e2.fvalue, -e2.gvalue)


@with_kw mutable struct AStarEpsilonMCSPStates{D <: Number}
    open_list::MutableBinaryMinHeap{AStarEpsilonMCSPHEntry{D}}              = MutableBinaryMinHeap{AStarEpsilonMCSPHEntry{D}}()
    closed_list::MutableBinaryMinHeap{AStarEpsilonMCSPHEntry{D}}            = MutableBinaryMinHeap{AStarEpsilonMCSPHEntry{D}}()
    open_list_hmap::Dict{Int64,Set{Int64}}                                  = Dict{Int64,Set{Int64}}()
    closed_list_hmap::Dict{Int64,Set{Int64}}                                = Dict{Int64,Set{Int64}}()
    # Subset of heap
    focal_heap::MutableBinaryHeap{AStarEpsilonMCSPHEntry{D},CompareFocalHeap}   = MutableBinaryHeap{AStarEpsilonMCSPHEntry{D},CompareFocalHeap}()
    ol_handles_in_focal::Set{Int64}                                         = Set{Int64}()
    focal_to_ol_handle::Dict{Int64,Int64}                                   = Dict{Int64,Int64}()
    best_fvalue::D = zero(D)
end

# A path is dominated if cost and weight are both worse than an existing entry
function is_dominated(state::AStarEpsilonMCSPStates{D}, s::AStarEpsilonMCSPHEntry{D}) where {D <: Number}

    if haskey(state.open_list_hmap, s.v_idx)
        for ol_idx in state.open_list_hmap[s.v_idx]
            ols_entry = state.open_list.nodes[state.open_list.node_map[ol_idx]].value
            if ols_entry.gvalue <= s.gvalue && ols_entry.weights <= s.weights
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
                     entry::AStarEpsilonMCSPHEntry{D},
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
    state::AStarEpsilonMCSPStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::F1,
    neighbors::Vector{Int64},
    parent_entry::AStarEpsilonMCSPHEntry{D},
    entry_cl_handle::Int64,    # For closed list
    visitor::AbstractDijkstraVisitor,
    eps_weight::Float64,
    admissible_heuristic::F2,
    focal_state_heuristic::F3,
    focal_transition_heuristic::F4,
    weight_functions::Vector{F5},
    weight_heuristics::Vector{F6},
    weight_constraints::Vector{D}
    ) where {V, D <: Number, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function, F5 <: Function, F6 <: Function}


    for nbr_idx in neighbors

        is_elig, new_weight_vector =  is_eligible(graph, parent_entry, nbr_idx, weight_functions, weight_heuristics, weight_constraints)

        if is_elig
            # Create new entry for open list
            # Updated gvalues, fvalues, weights AND focal stuff
            new_gvalue = parent_entry.gvalue + edge_wt_fn(graph.vertices[parent_entry.v_idx], graph.vertices[nbr_idx])
            new_fvalue = new_gvalue + admissible_heuristic(graph.vertices[nbr_idx])
            new_focal_heur = parent_entry.focal_heuristic + focal_state_heuristic(graph.vertices[nbr_idx]) +
                             focal_transition_heuristic(graph.vertices[parent_entry.v_idx], graph.vertices[nbr_idx])

            new_entry = AStarEpsilonMCSPHEntry(nbr_idx, parent_entry.v_idx, new_gvalue, new_fvalue,
                                                new_focal_heur, new_weight_vector, entry_cl_handle)

            # Insert into open list if not dominated
            # And accordingly into focal list if valid
            if ~(is_dominated(state, new_entry))
                new_entry_handle = push!(state.open_list, new_entry)
                if ~(haskey(state.open_list_hmap, nbr_idx))
                    Graphs.discover_vertex!(visitor, graph.vertices[parent_entry.v_idx], graph.vertices[nbr_idx], new_gvalue)
                    state.open_list_hmap[nbr_idx] = Set{Int64}(new_entry_handle)
                else
                    push!(state.open_list_hmap[nbr_idx], new_entry_handle)
                end

                # Now enter to focal list if valid according to weight
                if new_fvalue <= eps_weight * state.best_fvalue
                    new_focal_handle = push!(state.focal_heap, new_entry)
                    push!(state.ol_handles_in_focal, new_entry_handle)
                    state.focal_to_ol_handle[new_focal_handle] = new_entry_handle
                end # end if
            end # end if
        end
    end
end


function a_star_epsilon_constrained_shortest_path_implicit!(
    graph::AbstractGraph{V},                # the graph
    edge_wt_fn::F1, # distances associated with edges
    source::Int64,             # the source index
    visitor::AbstractDijkstraVisitor,# visitor object\
    eps_weight::Float64,
    admissible_heuristic::F2,      # Heuristic function for vertices
    focal_state_heuristic::F3,
    focal_transition_heuristic::F4,
    weight_functions::Vector{F5},
    weight_heuristics::Vector{F6},
    weight_constraints::Vector{D},
    ::Type{D} = Float64
    ) where {V, D <: Number, F1 <: Function, F2 <: Function, F3 <: Function, F4 <: Function, F5 <: Function, F6 <: Function}


    state = AStarEpsilonMCSPStates{D}()

    n_weights = length(weight_constraints)

    source_heur = admissible_heuristic(graph.vertices[source])
    state.best_fvalue = source_heur

    source_entry = AStarEpsilonMCSPHEntry(source, source, zero(D), source_heur,
                                         focal_state_heuristic(source), zeros(D, n_weights), 0)
    source_handle = push!(state.open_list, source_entry)
    state.open_list_hmap[source] = Set{Int64}(source_handle)
    focal_handle = push!(state.focal_heap, source_entry)
    state.focal_to_ol_handle[focal_handle] = source_handle

    while ~(isempty(state.open_list))

        old_best_fvalue = state.best_fvalue
        state.best_fvalue = top(state.open_list).fvalue

        # Update focal list
        if state.best_fvalue > old_best_fvalue
            # Iterate over open set  in increasing order of fvalue and insert in focal list if valid
            for node in sort(state.open_list.nodes, by = x->x.value.fvalue)
                fvalue = node.value.fvalue

                if fvalue > eps_weight * old_best_fvalue && fvalue <= eps_weight * state.best_fvalue &&
                    ~(node.handle in state.ol_handles_in_focal)
                    # Need to insert open list entry to focal list
                    # And assign node handle to focal handle
                    focal_handle = push!(state.focal_heap, node.value)
                    push!(state.ol_handles_in_focal, node.handle)
                    state.focal_to_ol_handle[focal_handle] = node.handle
                end

                if fvalue > eps_weight * state.best_fvalue
                    break
                end
            end
        end

        focal_entry, focal_handle = top_with_handle(state.focal_heap)

        # Get corresponding open list handle and entry
        ol_handle = state.focal_to_ol_handle[focal_handle]

        # Remove from open list and focal list
        pop!(state.focal_heap)
        delete!(state.focal_to_ol_handle, focal_handle) # Don't really need this?
        # @show focal_handle
        # @show ol_handle
        delete!(state.open_list, ol_handle)
        delete!(state.open_list_hmap[focal_entry.v_idx], ol_handle)

        nbrs = Vector{Int64}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[focal_entry.parent_idx],
                                  graph.vertices[focal_entry.v_idx], focal_entry.gvalue, nbrs) == false
            return state, focal_entry
        end


        # Insert into closed list to get new handle
        cl_entry_handle = push!(state.closed_list, focal_entry)
        if ~(haskey(state.closed_list_hmap, focal_entry.v_idx))
            state.closed_list_hmap[focal_entry.v_idx] = Set{Int64}(cl_entry_handle)
        else
            push!(state.closed_list_hmap[focal_entry.v_idx], cl_entry_handle)
        end

        # Process nbrs
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs,
                                    focal_entry, cl_entry_handle, visitor, eps_weight,
                                    admissible_heuristic, focal_state_heuristic, focal_transition_heuristic,
                                    weight_functions, weight_heuristics, weight_constraints)
        Graphs.close_vertex!(visitor, graph.vertices[focal_entry.v_idx])
    end

    return state, source_entry
end


# Given the solution state for the A*EpsilonMCSP Heap, walk back using the closed list
# to actually get the shortest feasible path
function shortest_path_cost_weights(state::AStarEpsilonMCSPStates{D}, graph::AbstractGraph{V},
                                    source::Int64, target_entry::AStarEpsilonMCSPHEntry{D}) where {V, D <: Number}

    target_idx = target_entry.v_idx

    sp = [target_idx]
    weights = [target_entry.weights]
    costs = [target_entry.gvalue]

    # Walk back USING closed list heap
    curr_idx = target_idx
    curr_par_handle = target_entry.parent_handle

    while curr_idx != source

        cl_list_entry = state.closed_list.nodes[state.closed_list.node_map[curr_par_handle]].value

        curr_idx = cl_list_entry.v_idx
        curr_par_handle = cl_list_entry.parent_handle
        pushfirst!(sp, curr_idx)
        pushfirst!(weights, cl_list_entry.weights)
        pushfirst!(costs, cl_list_entry.gvalue)
    end

    return sp, costs, weights
end
