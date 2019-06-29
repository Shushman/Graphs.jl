"""
Custom graph structure which only has a list of vertices; edges generated implicitly.
Particularly useful for a_star_implicit_light.jl
"""
mutable struct SimpleVListGraph{V} <: AbstractGraph{V,Edge{V}}
    vertices::Vector{V}
end

function SimpleVListGraph{V}() where V
    return SimpleVListGraph(Vector{V}(undef,0))
end

function Graphs.is_directed(g::SimpleVListGraph{V}) where V
    return false
end

function Graphs.num_edges(g::SimpleVListGraph{V}) where V
    return 0
end

function Graphs.add_vertex!(g::SimpleVListGraph{V}, v::V) where V
    push!(g.vertices, v)
end

function Graphs.num_vertices(g::SimpleVListGraph{V}) where V
    return length(g.vertices)
end

function Graphs.vertex_index(g::SimpleVListGraph{V}, v::V) where V
    return findfirst(isequal(v), g.vertices)
end

"""
    remove_last_vertex!(g::SimpleVListGraph{V}) where V

Pop the last vertex of the graph. Not a Graphs function.
"""
function remove_last_vertex!(g::SimpleVListGraph{V}) where V
    pop!(g.vertices)
end
