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

m = 1000
obj = basiclu_object(m)
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
