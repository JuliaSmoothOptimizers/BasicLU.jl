var documenterSearchIndex = {"docs":
[{"location":"#BasicLU.jl","page":"Documentation","title":"BasicLU.jl","text":"","category":"section"},{"location":"","page":"Documentation","title":"Documentation","text":"Modules = [BasicLU]","category":"page"},{"location":"#BasicLU.LUFactor","page":"Documentation","title":"BasicLU.LUFactor","text":"F = LUFactor(dim::Int64)\n\n\n\n\n\n","category":"type"},{"location":"#BasicLU.getfactors-Tuple{LUFactor}","page":"Documentation","title":"BasicLU.getfactors","text":"L, U, p, q = getfactors(F::LUFactor) -> L::SparseMatrixCSC{Float64, Int64},\n                                              U::SparseMatrixCSC{Float64, Int64},\n                                              p::Vector{Int64},\n                                              q::Vector{Int64}\n\nExtract LU factors after fresh factorization. L is unit lower triangular, U is upper triangular and p and q are permutation vectors such that (ignoring round-off errors) B[p,q] = L * U when matrix B was factorizied.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.maxvolbasis-Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}}","page":"Documentation","title":"BasicLU.maxvolbasis","text":"basis, F = maxvolbasis(A::SparseMatrixCSC{Float64, Int64}; lindeptol::Float64=1e-8,\n                       volumetol::Float64=2.0, maxpass::Int64=2, verbose::Bool=true) -> Vector{Int64}, LUFactor\n\nFind a set of column indices for the matrix AI = [A I] such that AI[:,basis] is square and nonsingular and the number of slack columns in the basis is minimum (this is the row rank deficiency of A). Return the vector of column indices of AI which form the basis matrix and a LUFactor which holds a factorization of the basis matrix.\n\nMethod: Scale the slack columns of AI by lindeptol and try to find a maximum volume basis for this matrix by making at most maxpass calls to maxvolume. If verbose is true, then print the number of basis updates after each call.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.maxvolume","page":"Documentation","title":"BasicLU.maxvolume","text":"nupdate = maxvolume(F::LUFactor, A::SparseMatrixCSC{Float64, Int64}, basis::Vector{Int64}, volumetol::Float64=2.0) -> Int64\n\nGiven an initial basis such that A[:,basis] is square and nonsingular, make one pass over the nonbasic columns of A and pivot each column into the basis when it increases the absolute value of the determinant of the basis matrix by more than a factor volumetol. On return basis has been updated. Return the number of basis updates performed.\n\n\n\n\n\n","category":"function"},{"location":"#BasicLU.solve!-Tuple{LUFactor, SparseArrays.SparseVector{Float64, Int64}, Char}","page":"Documentation","title":"BasicLU.solve!","text":"solve!(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char) -> SparseVector{Float64, Int64}\n\nSolve linear system with sparse right-hand side. Solution overwrites rhs. trans must be 'T' for transposed solve or 'N' for forward solve.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.solve!-Tuple{LUFactor, Vector{Float64}, Char}","page":"Documentation","title":"BasicLU.solve!","text":"solve!(F::LUFactor, rhs::Vector{Float64}, trans::Char) -> Vector{Float64}\n\nSolve linear system with dense right-hand side. Solution overwrites rhs. trans must be 'T' for transposed solve or 'N' for forward solve.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.solve-Tuple{LUFactor, SparseArrays.SparseVector{Float64, Int64}, Char}","page":"Documentation","title":"BasicLU.solve","text":"x = solve(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char) -> SparseVector{Float64, Int64}\n\nSolve linear system with sparse right-hand side. rhs is not modified. trans must be 'T' for transposed solve or 'N' for forward solve.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.solve-Tuple{LUFactor, Vector{Float64}, Char}","page":"Documentation","title":"BasicLU.solve","text":"x = solve(F::LUFactor, rhs::Vector{Float64}, trans::Char) -> Vector{Float64}\n\nSolve linear system with dense right-hand side. rhs is not modified. trans must be 'T' for transposed solve or 'N' for forward solve.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.solve_for_update-Tuple{LUFactor, Int64}","page":"Documentation","title":"BasicLU.solve_for_update","text":"solve_for_update(F::LUFactor, pos::Int64; getsol::Bool=false) -> SparseVector{Float64, Int64}\n\nSolve transposed system in preparation to update the factorization. pos holds the column index of the factorized matrix to be replaced in the next call to update. When getsol = true, then the solution from the transposed solve with a unit vector as right-hand side is returned. Otherwise only the update is prepared.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.solve_for_update-Tuple{LUFactor, SparseArrays.SparseVector{Float64, Int64}}","page":"Documentation","title":"BasicLU.solve_for_update","text":"solve_for_update(F::LUFactor, newcol::SparseVector{Float64, Int64}; getsol::Bool=false) -> SparseVector{Float64, Int64}\n\nSolve forward system in preparation to update the factorization. newcol holds the column to be inserted into the factorized matrix in the next call to update. When getsol = true, then the solution from the forward solve with right-hand side newcol is returned. Otherwise only the update is prepared.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.update-Tuple{LUFactor, Float64}","page":"Documentation","title":"BasicLU.update","text":"update(F::LUFactor, pivot::Float64) -> Float64\n\nUpdate the factorization after a column modification. The column position and the new column must have been set in previous calls to solve_for_update.\n\npivot is the pivot element corresponding to the update operation; i.e. when column j of B is to be replaced by vector v, then pivot = (B\\v)[j]. The absolute difference between pivot and a recomputed version can be obtained with getinfo(F, :pivotError); this is also the return value. A pivot error larger than 1e-8, say, indicates numerical instability and suggests refactorization.\n\nAn error is thrown when the recomputed pivot element is below the absolute pivot tolerance. In this case no update is performed and the old factorization remains valid.\n\n\n\n\n\n","category":"method"},{"location":"#BasicLU.update-Tuple{LUFactor, Int64, SparseArrays.SparseVector{Float64, Int64}}","page":"Documentation","title":"BasicLU.update","text":"update(F::LUFactor, pos::Int, newcol::SparseVector) -> (lhs, piverr)\n\nUpdate the factorization by inserting column newcol at index pos.\n\n\n\n\n\n","category":"method"},{"location":"#LinearAlgebra.factorize-Tuple{LUFactor, SparseArrays.SparseMatrixCSC{Float64, Int64}}","page":"Documentation","title":"LinearAlgebra.factorize","text":"factorize(F::LUFactor, B::SparseMatrixCSC{Float64, Int64}; check::Bool=true)\n\nFactorize sparse matrix B, which must be square and have the dimension for which F was created.\n\nThe factorization method stops when all elements of the active submatrix are below the absolute pivot tolerance. In this case the L and U matrix are padded with columns of the identity matrix, yielding the factorization of a matrix B* which equals B except that the dependent columns are replaced by columns of the identity matrix. The number of actual pivot steps performed can be obtained with getinfo(F, :nPivot).\n\nWhen check = true, an error is thrown if the number of pivot steps is less than the dimension of B.\n\n\n\n\n\n","category":"method"}]
}
