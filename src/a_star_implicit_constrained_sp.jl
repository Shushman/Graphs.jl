## A-star implicit with one or more additive weight constraints on feasible paths
## Basically, find the shortest path (w.r.t cost) that satisfies a bunch of
## constraints that depend on other weight functions
struct AStarMCSPHEntry{D <: Number}
    v_idx::Int64
    gvalue::D
    fvalue::D
    weights::Vector{D}
    self_handle::Int # refers to
    parent_handle::Int64    # refers to the heap handle in CLOSED LIST that tracks this path
end

Base.isless(e1::AStarMCSPHEntry, e2::AStarMCSPHEntry) = (e1.fvalue, e1.gvalue) < (e2.fvalue, e2.gvalue) # Lexicographic

@with_kw mutable struct AStarMCSPStates{D <: Number}
    open_list::BinaryMinHeap{AStarMCSPHEntry{D}}    = BinaryMinHeap{AStarMCSPHEntry{D}}()
    closed_list::BinaryMinHeap{AStarMCSPHEntry{D}}  = BinaryMinHeap{AStarMCSPHEntry{D}}()
    open_list_hmap::Dict{Int64,Set{Int64}}          = Dict{Int64,Set{Int64}}()
    closed_list_hmap::Dict{Int64,Set{Int64}}        = Dict{Int64,Set{Int64}}()
end

function is_dominated(state::AStarMCSPStates{D}, s::AStarMCSPHEntry{D}) where {D <: Number}

    for ol_idx in state.open_list_hmap[s.v_idx]
        ols_entry = state.open_list[ol_idx]
        if isless(ols_entry, s)
            return true
        end
    end

    for cl_idx in state.closed_list_hmap[s.v_idx]
        cls_entry = state.closed_list[cl_idx]
        if isless(cls_entry, s)
            return true
        end
    end

    return false
end

function process_neighbors_implicit!(
    state::AStarMCSPStates{D},
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    neighbors::Vector{Int64},
    entry::AStarMCSPHEntry{D},
    entry_handle::Int64,    # For closed list
    visitor::AbstractDijkstraVisitor,
    heuristic::Function,
    weight_heuristics::Vector{Function},
    weight_constraints::Vector{D}
    )

    for iv in neighbors

        # Iterate over neighbors and add if eligible and not dominated
        # Pay attention to parent handle etc.




function a_star_constrained_shortest_path_implicit!(
    graph::AbstractGraph{V},
    edge_wt_fn::Function,
    source::Int64,
    visitor::AbstractDijkstraVisitor,
    heuristic::Function,
    weight_heuristics::Vector{Function},
    weight_constraints::Vector{D},
    state::AStarMCSPStates{D}) where {V, D <: Number}

    n_weights = length(weight_constraints)

    source_entry = AStarMCSPHEntry(source, zero(D), heuristic(source), zeros(D, n_weights), 0, 0)
    push!(state.open_list, source_entry)

    while ~(isempty(state.open_list))
