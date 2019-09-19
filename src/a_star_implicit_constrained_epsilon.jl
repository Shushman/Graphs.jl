## NOTE - Uses code from src/a_star_implicit_constrained.jl

struct AStarEpsilonMCSPHEntry{D <: Number}
    a_star_entry::AStarMCSPHEntry{D}
    focal_heuristic::D
end

# Base.isless(e1::AStarEpsilonMCSPHEntry, e2::AStarEpsilonMCSPHEntry) = (e1.a_star_entry.fvalue, e1.a_star_entry.gvalue) < (e2.a_star_entry.fvalue, e2.a_star_entry.gvalue)

compare(comp::CompareFocalHeap, e1::AStarEpsilonMCSPHEntry, e2::AStarEpsilonMCSPHEntry) = (e1.focal_heuristic, e1.a_star_entry.fvalue, -e1.a_star_entry.gvalue) < (e2.focal_heuristic, e2.a_star_entry.fvalue, -e2.a_star_entry.gvalue)


@with_kw mutable struct AStarEpsilonMCSPStates{D <: Number}
    a_star_states::AStarMCSPStates{D}                                           = AStarMCSPStates{D}()
    # Subset of heap
    focal_heap::MutableBinaryHeap{AStarEpsilonMCSPHEntry{D},CompareFocalHeap}   = MutableBinaryHeap{AStarEpsilonMCSPHEntry{D},CompareFocalHeap}()
    ol_handles_in_focal::Set{Int64}                                             = Set{Int64}()
    focal_to_ol_handle::Dict{Int64,Int64}                                       = Dict{Int64,Int64}()
    best_fvalue::D = zero(D)
    focal_heur_value::Dict{Int64,D}                                         = Dict{Int64,D}()
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

        is_elig, new_weight_vector =  is_eligible(graph, parent_entry.a_star_entry, nbr_idx, weight_functions, weight_heuristics, weight_constraints)

        if is_elig
            # Create new entry for open list
            # Updated gvalues, fvalues, weights AND focal stuff
            new_gvalue = parent_entry.a_star_entry.gvalue + edge_wt_fn(graph.vertices[parent_entry.a_star_entry.v_idx], graph.vertices[nbr_idx])
            new_fvalue = new_gvalue + admissible_heuristic(graph.vertices[nbr_idx])
            new_focal_heur = parent_entry.focal_heuristic + focal_state_heuristic(graph.vertices[nbr_idx]) +
                             focal_transition_heuristic(graph.vertices[parent_entry.a_star_entry.v_idx], graph.vertices[nbr_idx])


            new_entry = AStarEpsilonMCSPHEntry(AStarMCSPHEntry(nbr_idx, parent_entry.a_star_entry.v_idx, new_gvalue, new_fvalue, new_weight_vector, entry_cl_handle),
                                                new_focal_heur)

            # Insert into open list if not dominated
            # And accordingly into focal list if valid
            if ~(is_dominated(state.a_star_states, new_entry.a_star_entry))
                new_entry_handle = push!(state.a_star_states.open_list, new_entry.a_star_entry)
                state.focal_heur_value[new_entry_handle] = new_focal_heur

                if ~(haskey(state.a_star_states.open_list_hmap, nbr_idx))
                    Graphs.discover_vertex!(visitor, graph.vertices[parent_entry.a_star_entry.v_idx], graph.vertices[nbr_idx], new_gvalue)
                    state.a_star_states.open_list_hmap[nbr_idx] = Set{Int64}(new_entry_handle)
                else
                    push!(state.a_star_states.open_list_hmap[nbr_idx], new_entry_handle)
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

"""
Runs the Focal-MCSP (Focal search with Multiple Constraints bounded-suboptimal Shortest Path) algorithm on the given list-of-indices graph.
As with all A* versions (except a_star_spath), the graph is implicit and the graph visitor must generate the neighbors of the expanded
vertex just-in-time.

Arguments:
    - `graph::AbstractGraph{V} where {V}`
    - `edge_wt_fn::F1 where {F1 <: Function}` Maps (u,v) to the edge weight of the u -> v edge
    - `source::Int64`
    - `visitor::AbstractDijkstraVisitor` The graph visitor that implements include_vertex!, which is called when
                                         a vertex is expanded
    - `eps_weight::Float64` w sub-optimality factor for focal search
    - `admissible_heuristic::F2 where {F2 <: Function}` Maps u to h(u); an admissible heuristic for the cost-to-go
    - `focal_state_heuristic::F3 where {F3 <: Function}`
    - `focal_transition_heuristic::F4 where {F4 <: Function}`
    - `weight_functions::Vector{F5} where {F5 <: Function}` A list of constraint weight functions for edges
    - `weight_heuristics::Vector{F6} where {F6 <: Function}` A list of admissible heuristics for the weight constraint for edges
    - `weight_constraints::Vector{D} where {D <: Number}` A list of correspoding weight constraint values
    - `::Type{D} = Float64` The type of the edge weight value

Returns:
    - `::AStarMCSPEpsilonStates` The result of the A* search
    - `::AStarMCSP` The MCSP heap entry for the target node (or the source if no path found)
"""
function a_star_implicit_constrained_epsilon_path!(
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

    source_entry = AStarEpsilonMCSPHEntry(AStarMCSPHEntry(source, source, zero(D), source_heur, zeros(D, n_weights), 0),
                                         focal_state_heuristic(source))

    source_handle = push!(state.a_star_states.open_list, source_entry.a_star_entry)
    state.a_star_states.open_list_hmap[source] = Set{Int64}(source_handle)
    state.focal_heur_value[source_handle] = focal_state_heuristic(source)

    focal_handle = push!(state.focal_heap, source_entry)
    push!(state.ol_handles_in_focal, source_handle)
    state.focal_to_ol_handle[focal_handle] = source_handle

    while ~(isempty(state.a_star_states.open_list))

        old_best_fvalue = state.best_fvalue
        state.best_fvalue = top(state.a_star_states.open_list).fvalue

        # Update focal list
        if state.best_fvalue > old_best_fvalue
            # Iterate over open set  in increasing order of fvalue and insert in focal list if valid
            for node in sort(state.a_star_states.open_list.nodes, by = x->x.value.fvalue)
                fvalue = node.value.fvalue

                if fvalue > eps_weight * old_best_fvalue && fvalue <= eps_weight * state.best_fvalue &&
                    ~(node.handle in state.ol_handles_in_focal)
                    # Need to insert open list entry to focal list
                    # And assign node handle to focal handle
                    new_focal_entry = AStarEpsilonMCSPHEntry(node.value, state.focal_heur_value[node.handle])
                    focal_handle = push!(state.focal_heap, new_focal_entry)
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

        # @show state.focal_heap

        # Remove from open list and focal list
        pop!(state.focal_heap)
        delete!(state.focal_to_ol_handle, focal_handle) # Don't really need this?
        delete!(state.a_star_states.open_list, ol_handle)
        delete!(state.a_star_states.open_list_hmap[focal_entry.a_star_entry.v_idx], ol_handle)

        nbrs = Vector{Int64}(undef, 0)

        if Graphs.include_vertex!(visitor, graph.vertices[focal_entry.a_star_entry.parent_idx],
                                  graph.vertices[focal_entry.a_star_entry.v_idx], focal_entry.a_star_entry.gvalue, nbrs) == false
            return state, focal_entry.a_star_entry
        end


        # Insert into closed list to get new handle
        cl_entry_handle = push!(state.a_star_states.closed_list, focal_entry.a_star_entry)
        if ~(haskey(state.a_star_states.closed_list_hmap, focal_entry.a_star_entry.v_idx))
            state.a_star_states.closed_list_hmap[focal_entry.a_star_entry.v_idx] = Set{Int64}(cl_entry_handle)
        else
            push!(state.a_star_states.closed_list_hmap[focal_entry.a_star_entry.v_idx], cl_entry_handle)
        end

        # Process nbrs
        process_neighbors_implicit!(state, graph, edge_wt_fn, nbrs,
                                    focal_entry, cl_entry_handle, visitor, eps_weight,
                                    admissible_heuristic, focal_state_heuristic, focal_transition_heuristic,
                                    weight_functions, weight_heuristics, weight_constraints)
        Graphs.close_vertex!(visitor, graph.vertices[focal_entry.a_star_entry.v_idx])
    end

    return state, source_entry.a_star_entry
end
