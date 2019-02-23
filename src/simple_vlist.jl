"""
Custom graph structure which only has a list of vertices; edges generated implicitly.
"""
mutable struct SimpleVListGraph{V} <: AbstractGraph{V,Edge{V}}
    vertices::Vector{V}
end

function Graphs.add_vertex!(g::SimpleVListGraph{V}, v::V) where V
    push!(g.vertices, v)
end

function Graphs.num_vertices(g::SimpleVListGraph{V}) where V
    return length(g.vertices)
end

"""
Pop the last vertex of the graph. Not a Graphs function.
"""
function remove_last_vertex!(g::SimpleVListGraph{V}) where V
    pop!(g.vertices)
end