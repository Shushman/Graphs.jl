using Test

tests = [
    # "edgelist",
    # "adjlist",
    # "bellman_test",
    # "inclist",
    # "graph",
    # "gmatrix",
    # "bfs",
    # "dfs",
    # "conn_comp",
    # "dijkstra",
    # "a_star_spath",
    # "a_star_implicit_test",
    # "a_star_implicit_constrained_test",
    "a_star_implicit_constrained_epsilon_test",
    # "mst",
    # "floyd",
    # "dot",
    # "dot2",
    # "cliques",
    # "random",
    # "generators",
    # "maximum_adjacency_visit",
    # "issue_related_tests",
    # "inclist_dict_delete"
    ]

testdir = joinpath(dirname(@__DIR__), "test")

for t in tests
    @testset "$(t)" begin
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
