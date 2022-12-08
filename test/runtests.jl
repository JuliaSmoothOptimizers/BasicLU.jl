using LinearAlgebra, SparseArrays, MAT
using Printf, Test
using BasicLU

@info("BASICLU_INSTALLATION: $(BasicLU.BASICLU_INSTALLATION)")

include("utils.jl")

test_small()
test_factorize("./data/", false)
test_factorize("./data/", true)
test_update("./data/")
test_maxvolume("./data/")

@testset "Basic test" begin
    m = 1000
    obj = LUFactor(m)
    B = sprand(m, m, 5e-3) + I
    err = factorize(obj, B)
    rhs = randn(m)
    lhs = solve(obj, rhs, 'N')
    res = norm(B * lhs - rhs, Inf)
    @test res ≤ √eps(Float64)
    col = sparsevec([1], [1.0], m)
    lhs = solve_for_update(obj, col, getsol=true)
    vmax, j = findmax(abs.(lhs))
    piv = lhs[j]
    solve_for_update(obj, j)
    piverr = update(obj, piv)
    lhs = solve(obj, rhs, 'N')
    B[:,j] = col
    res = norm(B * lhs - rhs, Inf)
    @test res ≤ √eps(Float64)
end

@testset "Integrated rank-one update" begin
    m = 500
    for _ in 1:100
        F = LUFactor(m)
        B = sprand(m, m, 5e-3) + 3I
        err = factorize(F, B)
        rhs = randn(m)
        lhs = solve(F, rhs, 'N')
        res = norm(B * lhs - rhs, Inf)
        @test res ≤ sqrt(eps())
        ridx1, ridx2 = rand(1:m÷2), rand(m÷2+1:m)
        v1, v2 = rand(), rand()
        col = sparsevec([ridx1], [v1], m)
        _, piverr = BasicLU.update(F, ridx1, col)
        lhs = solve(F, rhs, 'N')
        B[:,ridx1] .= col
        @test norm(B * lhs - rhs, Inf) <= sqrt(eps())
    end
end
