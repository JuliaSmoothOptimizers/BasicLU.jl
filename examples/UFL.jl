using LinearAlgebra, SparseArrays
using LinearOperators, Printf, Krylov
using DelimitedFiles, MatrixMarket

dict = Dict("ge" => true)

atol = 0.0
rtol = 1e-12
itmax = 0
verbose = 0

using Plots
pgfplotsx()
using LaTeXStrings

for pb in keys(dict)
  if dict[pb]
    println("Problème $pb")
  
    println("Lecture des données: ")
    A = MatrixMarket.mmread("./data/" * pb * ".mtx")
    m, n = size(A)
    b = randn(n)
    println("Dimension de A: ($m, $n)")
    println("Nombre de nnz: ", length(A.nzval))
    basis = readdlm("./data/" * pb * "-basis.txt", Int64)[:]
    notbasis = setdiff(1:n, basis)
    B = sparse(A[:, basis]')
    bB = b[basis]
    N = sparse(A[:, notbasis]')
    bN = b[notbasis]
    println("✔")

    println("Résolution du système linéaire: ")
    _, stats  = lsqr(A', b, atol=atol, rtol=rtol, history=true, itmax=itmax, verbose=verbose)
    _, stats2 = lsmr(A', b, atol=atol, rtol=rtol, history=true, itmax=itmax, verbose=verbose)

    F = lu(B)
    M = LinearOperator(Float64, m, m, false, false, (y, v, α, β) -> (y .= v; ldiv!(F, y) ; ldiv!(F', y)))
    _, _, stats3 = tricg(N, bN, -B' * bB, N=M, atol=atol, rtol=rtol, history=true, itmax=itmax, verbose=verbose)
    _, _, stats4 = trimr(N, bN, -B' * bB, N=M, atol=atol, rtol=rtol, history=true, itmax=itmax, verbose=verbose)
    println("✔")

    print("Génération des graphiques: ")
#=    plot(0:length(stats.residuals)-1, stats.residuals, label="LSQR", yaxis=:log10, lw=0.5, color=2, xlabel=L"k", ylabel=L"\|r_k\|", title=pb, legend=:best)
    plot!(0:length(stats2.residuals)-1, stats2.residuals, label="LSMR", yaxis=:log10, lw=0.5, color=3)
    savefig("graphics/residuals_$(pb)" * ".svg")=#

    plot(0:length(stats.Aresiduals)-1, stats.Aresiduals, label="LSQR", yaxis=:log10, lw=0.5, color=2, xlabel=L"k", ylabel=L"\|A r_k\|", title=pb, legend=:best)
    plot!(0:length(stats2.Aresiduals)-1, stats2.Aresiduals, label="LSMR", yaxis=:log10, lw=0.5, color=3)
    savefig("graphics/Aresiduals_$(pb)" * ".svg")

    plot(0:length(stats3.residuals)-1, stats3.residuals, label="TriCG", yaxis=:log10, lw=0.5, color=2, xlabel=L"k", ylabel=L"\|r_k\|", title=pb, legend=:best)
    plot!(0:length(stats4.residuals)-1, stats4.residuals, label="TriMR", yaxis=:log10, lw=0.5, color=3)
    savefig("graphics/residuals_$(pb)" * ".svg")
    println("✔")
  end
end
