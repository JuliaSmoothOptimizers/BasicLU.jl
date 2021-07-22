using LinearAlgebra, SparseArrays, MAT
using Printf, Test
using BasicLU

include("utils.jl")

test_factorize("./data/", false)
test_factorize("./data/", true)
test_update("./data/")

m = 1000
this = BasicLU.initialize(m)
B = sprand(m, m, 5e-3) + I
err = BasicLU.factorize(this, B)
if err ≠ BasicLU.BASICLU_OK
    error("factorization failed or singular")
end
rhs = randn(m)
lhs = solve(this, rhs, 'N')
res = norm(B * lhs - rhs, Inf)
@test res ≤ √eps(cdbl)
col = sparsevec([1], [1.0], m)
lhs = solve4update(this, col, true)
(vmax,j) = findmax(abs.(lhs))
xtbl = lhs[j]
solve4update(this, j)
piverr = update(this, xtbl)
lhs = solve(this, rhs, 'N')
B[:,j] = col
res = norm(B * lhs - rhs, Inf)
@test res ≤ √eps(cdbl)
