#
# Julia interface to BASICLU
#

module BasicLU

using LinearAlgebra
using LinearAlgebra: checksquare
using Printf
using SparseArrays
using SparseArrays: getcolptr
using Libdl

import Base.@kwdef
import LinearAlgebra.factorize

if haskey(ENV, "JULIA_BASICLU_LIBRARY_PATH")
  const libbasiclu = joinpath(ENV["JULIA_BASICLU_LIBRARY_PATH"], "libbasiclu.$dlext")
  const BASICLU_INSTALLATION = "CUSTOM"
else
  using basiclu_jll
  const BASICLU_INSTALLATION = "YGGDRASIL"
end

export Param, Info, LUFactor
export paramnumber, infonumber, getdim, getparam, getinfo, setparam!
export factorize, getfactors, solve, solve!, solve_for_update, update, maxvolume, maxvolbasis

include("deprecated.jl")

# status codes
const BASICLU_OK = 0
const BASICLU_WARNING_singular_matrix = 2
const BASICLU_ERROR_invalid_store = -1
const BASICLU_ERROR_invalid_call = -2
const BASICLU_ERROR_argument_missing = -3
const BASICLU_ERROR_invalid_argument = -4
const BASICLU_ERROR_maximum_updates = -5
const BASICLU_ERROR_singular_update = -6
const BASICLU_ERROR_invalid_object = -8
const BASICLU_ERROR_out_of_memory = -9

# user parameters in xstore (zero based)
const BASICLU_DROP_TOLERANCE = 4
const BASICLU_ABS_PIVOT_TOLERANCE = 5
const BASICLU_REL_PIVOT_TOLERANCE = 6
const BASICLU_BIAS_NONZEROS = 7
const BASICLU_MAXN_SEARCH_PIVOT = 8
const BASICLU_PAD = 9
const BASICLU_STRETCH = 10
const BASICLU_COMPRESSION_THRESHOLD = 11
const BASICLU_SPARSE_THRESHOLD = 12
const BASICLU_REMOVE_COLUMNS = 13
const BASICLU_SEARCH_ROWS = 14

# user readable from xstore (zero based)
const BASICLU_DIM = 64
const BASICLU_STATUS = 65
const BASICLU_NUPDATE = 70
const BASICLU_NFORREST = 71
const BASICLU_NFACTORIZE = 72
const BASICLU_NUPDATE_TOTAL = 73
const BASICLU_NFORREST_TOTAL = 74
const BASICLU_NSYMPERM_TOTAL = 75
const BASICLU_LNZ = 76
const BASICLU_UNZ = 77
const BASICLU_RNZ = 78
const BASICLU_MIN_PIVOT = 79
const BASICLU_MAX_PIVOT = 80
const BASICLU_UPDATE_COST = 81
const BASICLU_TIME_FACTORIZE = 82
const BASICLU_TIME_SOLVE = 83
const BASICLU_TIME_UPDATE = 84
const BASICLU_TIME_FACTORIZE_TOTAL = 85
const BASICLU_TIME_SOLVE_TOTAL = 86
const BASICLU_TIME_UPDATE_TOTAL = 87
const BASICLU_LFLOPS = 88
const BASICLU_UFLOPS = 89
const BASICLU_RFLOPS = 90
const BASICLU_CONDEST_L = 91
const BASICLU_CONDEST_U = 92
const BASICLU_MAX_ETA = 93
const BASICLU_NORM_L = 94
const BASICLU_NORM_U = 95
const BASICLU_NORMEST_LINV = 96
const BASICLU_NORMEST_UINV = 97
const BASICLU_MATRIX_ONENORM = 98
const BASICLU_MATRIX_INFNORM = 99
const BASICLU_RESIDUAL_TEST = 111
const BASICLU_MATRIX_NZ = 100
const BASICLU_RANK = 101
const BASICLU_BUMP_SIZE = 102
const BASICLU_BUMP_NZ = 103
const BASICLU_NSEARCH_PIVOT = 104
const BASICLU_NEXPAND = 105
const BASICLU_NGARBAGE = 106
const BASICLU_FACTOR_FLOPS = 107
const BASICLU_TIME_SINGLETONS = 108
const BASICLU_TIME_SEARCH_PIVOT = 109
const BASICLU_TIME_ELIM_PIVOT = 110
const BASICLU_PIVOT_ERROR = 120

const basiclu_object_type = Cvoid

@kwdef mutable struct Param
    # Factorization parameters
    absPivotTol::Float64 = 0.0
    relPivotTol::Float64 = 0.0
    dropTol::Float64 = 0.0
    biasNz::Int64 = 0
    maxSearch::Int64 = 0
    searchRows::Int64 = 0
    removeCols::Int64 = 0
    memPad::Int64 = 0
    memStretch::Float64 = 0.0
    memCompress::Float64 = 0.0
    # Solve parameters
    sparseThres::Float64 = 0.0
end

@kwdef mutable struct Info
    dim::Int64 = 0
    status::Int64 = 0
    # Info from last factorization
    nPivot::Int64 = 0
    nMatElem::Int64 = 0
    nBumpCol::Int64 = 0
    nBumpElem::Int64 = 0
    nSearchPivot::Int64 = 0
    nExpand::Int64 = 0
    nGarbage::Int64 = 0
    nFactorFlop::Int64 = 0
    timeFactorize::Float64 = 0.0
    timeSingletons::Float64 = 0.0
    timeSearchPivot::Float64 = 0.0
    timeElimPivot::Float64 = 0.0
    residualTest::Float64 = 0.0
    onenormMatrix::Float64 = 0.0
    infnormMatrix::Float64 = 0.0
    normL::Float64 = 0.0
    normU::Float64 = 0.0
    normestLinv::Float64 = 0.0
    normestUinv::Float64 = 0.0
    condestL::Float64 = 0.0
    condestU::Float64 = 0.0
    # Info from solves and updates since last factorization
    nUpdate::Int64 = 0
    nForrest::Int64 = 0
    nElemL::Int64 = 0
    nElemU::Int64 = 0
    nElemR::Int64 = 0
    nFlopL::Int64 = 0
    nFlopU::Int64 = 0
    nFlopR::Int64 = 0
    updateCost::Float64 = 0.0
    minPivot::Float64 = 0.0
    maxPivot::Float64 = 0.0
    maxEta::Float64 = 0.0
    pivotError::Float64 = 0.0
    timeSolve::Float64 = 0.0
    timeUpdate::Float64 = 0.0
    # Accumulated info since object was created
    nFactorize::Int64 = 0
    nUpdateTotal::Int64 = 0
    nForrestTotal::Int64 = 0
    nSympermTotal::Int64 = 0
    timeFactorizeTotal::Float64 = 0.0
    timeSolveTotal::Float64 = 0.0
    timeUpdateTotal::Float64 = 0.0
end

function paramnumber(name::Symbol)
    if name == :absPivotTol
        return BASICLU_ABS_PIVOT_TOLERANCE
    elseif name == :relPivotTol
        return BASICLU_REL_PIVOT_TOLERANCE
    elseif name == :dropTol
        return BASICLU_DROP_TOLERANCE
    elseif name == :biasNz
        return BASICLU_BIAS_NONZEROS
    elseif name == :maxSearch
        return BASICLU_MAXN_SEARCH_PIVOT
    elseif name == :searchRows
        return BASICLU_SEARCH_ROWS
    elseif name == :removeCols
        return BASICLU_REMOVE_COLUMNS
    elseif name == :memPad
        return BASICLU_PAD
    elseif name == :memStretch
        return BASICLU_STRETCH
    elseif name == :memCompress
        return BASICLU_COMPRESSION_THRESHOLD
    elseif name == :sparseThres
        return BASICLU_SPARSE_THRESHOLD
    end
    throw(ErrorException("unknown parameter name"))
end

function infonumber(name::Symbol)
    if name == :dim
        return BASICLU_DIM
    elseif name == :status
        return BASICLU_STATUS
    elseif name == :nPivot
        return BASICLU_RANK
    elseif name == :nMatElem
        return BASICLU_MATRIX_NZ
    elseif name == :nBumpCol
        return BASICLU_BUMP_SIZE
    elseif name == :nBumpElem
        return BASICLU_BUMP_NZ
    elseif name == :nSearchPivot
        return BASICLU_NSEARCH_PIVOT
    elseif name == :nExpand
        return BASICLU_NEXPAND
    elseif name == :nGarbage
        return BASICLU_NGARBAGE
    elseif name == :nFactorFlop
        return BASICLU_FACTOR_FLOPS
    elseif name == :timeFactorize
        return BASICLU_TIME_FACTORIZE
    elseif name == :timeSingletons
        return BASICLU_TIME_SINGLETONS
    elseif name == :timeSearchPivot
        return BASICLU_TIME_SEARCH_PIVOT
    elseif name == :timeElimPivot
        return BASICLU_TIME_ELIM_PIVOT
    elseif name == :residualTest
        return BASICLU_RESIDUAL_TEST
    elseif name == :onenormMatrix
        return BASICLU_MATRIX_ONENORM
    elseif name == :infnormMatrix
        return BASICLU_MATRIX_INFNORM
    elseif name == :normL
        return BASICLU_NORM_L
    elseif name == :normU
        return BASICLU_NORM_U
    elseif name == :normestLinv
        return BASICLU_NORMEST_LINV
    elseif name == :normestUinv
        return BASICLU_NORMEST_UINV
    elseif name == :condestL
        return BASICLU_CONDEST_L
    elseif name == :condestU
        return BASICLU_CONDEST_U
    elseif name == :nUpdate
        return BASICLU_NUPDATE
    elseif name == :nForrest
        return BASICLU_NFORREST
    elseif name == :nElemL
        return BASICLU_LNZ
    elseif name == :nElemU
        return BASICLU_UNZ
    elseif name == :nElemR
        return BASICLU_RNZ
    elseif name == :nFlopL
        return BASICLU_LFLOPS
    elseif name == :nFlopU
        return BASICLU_UFLOPS
    elseif name == :nFlopR
        return BASICLU_RFLOPS
    elseif name == :updateCost
        return BASICLU_UPDATE_COST
    elseif name == :minPivot
        return BASICLU_MIN_PIVOT
    elseif name == :maxPivot
        return BASICLU_MAX_PIVOT
    elseif name == :maxEta
        return BASICLU_MAX_ETA
    elseif name == :pivotError
        return BASICLU_PIVOT_ERROR
    elseif name == :timeSolve
        return BASICLU_TIME_SOLVE
    elseif name == :timeUpdate
        return BASICLU_TIME_UPDATE
    elseif name == :nFactorize
        return BASICLU_NFACTORIZE
    elseif name == :nUpdateTotal
        return BASICLU_NUPDATE_TOTAL
    elseif name == :nForrestTotal
        return BASICLU_NFORREST_TOTAL
    elseif name == :nSympermTotal
        return BASICLU_NSYMPERM_TOTAL
    elseif name == :timeFactorizeTotal
        return BASICLU_TIME_FACTORIZE_TOTAL
    elseif name == :timeSolveTotal
        return BASICLU_TIME_SOLVE_TOTAL
    elseif name == :timeUpdateTotal
        return BASICLU_TIME_UPDATE_TOTAL
    end
    throw(ErrorException("unknown info name"))
end

"""
    F = LUFactor(dim::Int64)
"""
mutable struct LUFactor
    istore::Ptr{Int64}
    xstore::Ptr{Float64}
    Li::Ptr{Int64}
    Ui::Ptr{Int64}
    Wi::Ptr{Int64}
    Lx::Ptr{Float64}
    Ux::Ptr{Float64}
    Wx::Ptr{Float64}
    lhs::Ptr{Float64}
    ilhs::Ptr{Int64}
    nzlhs::Int64
    realloc_factor::Float64
    function LUFactor(dim::Int64)
        F = new()
        retcode = ccall((:basiclu_obj_initialize, libbasiclu), Int64,
                        (Ptr{basiclu_object_type}, Int64),
                        pointer_from_objref(F), dim)
        checkretcode("basiclu_obj_initialize", retcode)
        finalizer(F) do x
            ccall((:basiclu_obj_free, libbasiclu), Cvoid,
                  (Ptr{basiclu_object_type},), pointer_from_objref(x))
        end
    end
end

function getdim(F::LUFactor)
    ccall((:basiclu_obj_get_dim, libbasiclu), Int64,
          (Ptr{basiclu_object_type},), pointer_from_objref(F))
end

function getparam(F::LUFactor, name::Symbol)
    val = 0.0
    if F.xstore != C_NULL
        val = unsafe_load(F.xstore, paramnumber(name) + 1)
    end
    convert(fieldtype(Param, name), val)
end

function getparam(F::LUFactor)
    param = Param()
    for f in fieldnames(Param)
        setfield!(param, f, getparam(F, f))
    end
    param
end

function setparam!(F::LUFactor, name::Symbol, val)
    if F.xstore != C_NULL
        unsafe_store!(F.xstore, val, paramnumber(name) + 1)
    end
end

function setparam!(F::LUFactor, param::Param)
    for f in fieldnames(Param)
        setparam!(F, f, getfield(param, f))
    end
end

function getinfo(F::LUFactor, name::Symbol)
    val = 0.0
    if F.xstore != C_NULL
        val = unsafe_load(F.xstore, infonumber(name) + 1)
    end
    convert(fieldtype(Info, name), val)
end

function getinfo(F::LUFactor)
    info = Info()
    for f in fieldnames(Info)
        setfield!(info, f, getinfo(F, f))
    end
    info
end

"""
    factorize(F::LUFactor, B::SparseMatrixCSC{Float64, Int64}; check::Bool=true)

Factorize sparse matrix `B`, which must be square and have the dimension for
which `F` was created.

The factorization method stops when all elements of the active submatrix are
below the absolute pivot tolerance. In this case the `L` and `U` matrix are
padded with columns of the identity matrix, yielding the factorization of a
matrix `B*` which equals `B` except that the dependent columns are replaced by
columns of the identity matrix. The number of actual pivot steps performed can
be obtained with `getinfo(F, :nPivot)`.

When `check = true`, an error is thrown if the number of pivot steps is less
than the dimension of `B`.
"""
function factorize(F::LUFactor, B::SparseMatrixCSC{Float64, Int64}; check::Bool=true)
    dim = getdim(F)
    m = checksquare(B)
    if dim != m
        throw(DimensionMismatch("matrix dimension does not match basiclu object"))
    end
    if dim != 0
        Bp = getcolptr(B) .- 1
        Bi = rowvals(B) .- 1
        Bx = nonzeros(B)        # don't need a copy
        retcode = ccall((:basiclu_obj_factorize, libbasiclu), Int64,
                        (Ptr{basiclu_object_type}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
                        pointer_from_objref(F), Bp, pointer(Bp, 2), Bi, Bx)
        if retcode == BASICLU_WARNING_singular_matrix && !check
            retcode = BASICLU_OK
        end
        checkretcode("basiclu_obj_factorize", retcode)
    end
end

"""
    L, U, p, q = getfactors(F::LUFactor) -> L::SparseMatrixCSC{Float64, Int64},
                                                  U::SparseMatrixCSC{Float64, Int64},
                                                  p::Vector{Int64},
                                                  q::Vector{Int64}

Extract LU factors after fresh factorization. `L` is unit lower triangular, `U`
is upper triangular and `p` and `q` are permutation vectors such that (ignoring
round-off errors) `B[p,q] = L * U` when matrix `B` was factorizied.
"""
function getfactors(F::LUFactor)
    if getinfo(F, :nUpdate) != 0
        throw(ErrorException("cannot extract LU factors after factorization has been updated"))
    end
    dim = getdim(F)
    if dim == 0
        L = spzeros(0, 0)
        U = spzeros(0, 0)
        rowperm = zeros(Int64, 0)
        colperm = zeros(Int64, 0)
        return L, U, rowperm, colperm
    end
    info = getinfo(F)
    rowperm = Vector{Int64}(undef, dim)
    colperm = Vector{Int64}(undef, dim)
    Lcolptr = Vector{Int64}(undef, dim + 1)
    Lrowidx = Vector{Int64}(undef, info.nElemL + dim)
    Lvalue = Vector{Float64}(undef, info.nElemL + dim)
    Ucolptr = Vector{Int64}(undef, dim + 1)
    Urowidx = Vector{Int64}(undef, info.nElemU + dim)
    Uvalue = Vector{Float64}(undef, info.nElemU + dim)
    retcode = ccall((:basiclu_obj_get_factors, libbasiclu), Int64,
                    (Ptr{basiclu_object_type}, Ptr{Int64}, Ptr{Int64},
                     Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
                     Ptr{Int64}, Ptr{Int64}, Ptr{Float64}),
                    pointer_from_objref(F), rowperm, colperm,
                    Lcolptr, Lrowidx, Lvalue,
                    Ucolptr, Urowidx, Uvalue)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no fresh factorization available"))
    end
    checkretcode("basiclu_get_factors", retcode)
    L = SparseMatrixCSC{Float64, Int64}(dim, dim, Lcolptr .+ 1, Lrowidx .+ 1, Lvalue)
    U = SparseMatrixCSC{Float64, Int64}(dim, dim, Ucolptr .+ 1, Urowidx .+ 1, Uvalue)
    return L, U, rowperm .+ 1, colperm .+ 1
end

"""
    x = solve(F::LUFactor, rhs::Vector{Float64}, trans::Char) -> Vector{Float64}

Solve linear system with dense right-hand side. `rhs` is not modified.
`trans` must be `'T'` for transposed solve or `'N'` for forward solve.
"""
function solve(F::LUFactor, rhs::Vector{Float64}, trans::Char)
    lhs = copy(rhs)
    solve!(F, lhs, trans)
end

"""
    solve!(F::LUFactor, rhs::Vector{Float64}, trans::Char) -> Vector{Float64}

Solve linear system with dense right-hand side. Solution overwrites `rhs`.
`trans` must be `'T'` for transposed solve or `'N'` for forward solve.
"""
function solve!(F::LUFactor, rhs::Vector{Float64}, trans::Char)
    checktrans(trans)
    dim = getdim(F)
    if length(rhs) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim != 0
        retcode = ccall((:basiclu_obj_solve_dense, libbasiclu), Int64,
                        (Ptr{basiclu_object_type},  Ptr{Float64}, Ptr{Float64}, Cchar),
                        pointer_from_objref(F), rhs, rhs, trans)
        if retcode == BASICLU_ERROR_invalid_call
            throw(ErrorException("no factorization available"))
        end
        checkretcode("basiclu_obj_solve_dense", retcode)
    end
    rhs
end

"""
    x = solve(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char) -> SparseVector{Float64, Int64}

Solve linear system with sparse right-hand side. `rhs` is not modified.
`trans` must be `'T'` for transposed solve or `'N'` for forward solve.
"""
function solve(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char)
    lhs = copy(rhs)
    solve!(F, lhs, trans)
end

"""
    solve!(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char) -> SparseVector{Float64, Int64}

Solve linear system with sparse right-hand side. Solution overwrites `rhs`.
`trans` must be `'T'` for transposed solve or `'N'` for forward solve.
"""
function solve!(F::LUFactor, rhs::SparseVector{Float64, Int64}, trans::Char)
    checktrans(trans)
    dim = getdim(F)
    if length(rhs) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim == 0
        return spzeros(dim)
    end
    retcode = ccall((:basiclu_obj_solve_sparse, libbasiclu), Int64,
                    (Ptr{basiclu_object_type},  Int64, Ptr{Int64}, Ptr{Float64}, Cchar),
                    pointer_from_objref(F), nnz(rhs), rowvals(rhs) .- 1, nonzeros(rhs), trans)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    checkretcode("basiclu_obj_solve_sparse", retcode)
    gathersol!(F, rhs)
end

"""
    solve_for_update(F::LUFactor, newcol::SparseVector{Float64, Int64}; getsol::Bool=false) -> SparseVector{Float64, Int64}

Solve forward system in preparation to update the factorization. `newcol` holds the
column to be inserted into the factorized matrix in the next call to
[`update`](@ref). When `getsol = true`, then the solution from the forward
solve with right-hand side `newcol` is returned. Otherwise only the update is
prepared.
"""
function solve_for_update(F::LUFactor, newcol::SparseVector{Float64, Int64}; getsol::Bool=false)
    dim = getdim(F)
    if length(newcol) != dim
        throw(DimensionMismatch("dimension of right-hand side does not match basiclu object"))
    end
    if dim == 0
        return getsol ? spzeros(dim) : nothing
    end
    retcode = ccall((:basiclu_obj_solve_for_update, libbasiclu), Int64,
                    (Ptr{basiclu_object_type}, Int64, Ptr{Int64}, Ptr{Float64}, Cchar, Int64),
                    pointer_from_objref(F), nnz(newcol), rowvals(newcol) .- 1, nonzeros(newcol), 'N', getsol ? 1 : 0)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    if retcode == BASICLU_ERROR_maximum_updates
        throw(ErrorException("maximum number of updates reached"))
    end
    checkretcode("basiclu_obj_solve_for_update", retcode)
    return getsol ? gathersol(F) : nothing
end

"""
    solve_for_update(F::LUFactor, pos::Int64; getsol::Bool=false) -> SparseVector{Float64, Int64}

Solve transposed system in preparation to update the factorization. `pos` holds
the column index of the factorized matrix to be replaced in the next call to
[`update`](@ref). When `getsol = true`, then the solution from the transposed
solve with a unit vector as right-hand side is returned. Otherwise only the
update is prepared.
"""
function solve_for_update(F::LUFactor, pos::Int64; getsol::Bool=false)
    dim = getdim(F)
    if pos < 1 || pos > dim
        throw(DimensionMismatch("column index outside dimension of basiclu object"))
    end
    if dim == 0
        return getsol ? spzeros(dim) : nothing
    end
    retcode = ccall((:basiclu_obj_solve_for_update, libbasiclu), Int64,
                    (Ptr{basiclu_object_type}, Int64, Ptr{Int64}, Ptr{Float64}, Cchar, Int64),
                    pointer_from_objref(F), 0, Ref{Int64}(pos-1), C_NULL, 'T', getsol ? 1 : 0)
    if retcode == BASICLU_ERROR_invalid_call
        throw(ErrorException("no factorization available"))
    end
    if retcode == BASICLU_ERROR_maximum_updates
        throw(ErrorException("maximum number of updates reached"))
    end
    checkretcode("basiclu_obj_solve_for_update", retcode)
    return getsol ? gathersol(F) : nothing
end

"""
    update(F::LUFactor, pivot::Float64) -> Float64

Update the factorization after a column modification. The column position and
the new column must have been set in previous calls to
[`solve_for_update`](@ref).

`pivot` is the pivot element corresponding to the update operation; i.e. when
column `j` of `B` is to be replaced by vector `v`, then `pivot = (B\\v)[j]`. The
absolute difference between `pivot` and a recomputed version can be obtained
with `getinfo(F, :pivotError)`; this is also the return value. A pivot error
larger than 1e-8, say, indicates numerical instability and suggests
refactorization.

An error is thrown when the recomputed pivot element is below the absolute pivot
tolerance. In this case no update is performed and the old factorization remains
valid.
"""
function update(F::LUFactor, pivot::Float64)
    dim = getdim(F)
    if dim != 0
        retcode = ccall((:basiclu_obj_update, libbasiclu), Int64,
                        (Ptr{basiclu_object_type},  Float64),
                        pointer_from_objref(F), pivot)
        if retcode == BASICLU_ERROR_invalid_call
            throw(ErrorException("not prepared for update"))
        end
        checkretcode("basiclu_update", retcode)
    end
    getinfo(F, :pivotError)
end

"""
    nupdate = maxvolume(F::LUFactor, A::SparseMatrixCSC{Float64, Int64}, basis::Vector{Int64}, volumetol::Float64=2.0) -> Int64

Given an initial basis such that `A[:,basis]` is square and nonsingular, make
one pass over the nonbasic columns of `A` and pivot each column into the basis
when it increases the absolute value of the determinant of the basis matrix by
more than a factor `volumetol`. On return `basis` has been updated. Return the
number of basis updates performed.
"""
function maxvolume(F::LUFactor, A::SparseMatrixCSC{Float64, Int64}, basis::Vector{Int64},
                   volumetol::Float64=2.0)
    m, n = size(A)
    if m != getdim(F)
        throw(DimensionMismatch("number of rows of A does not match dimension of basiclu object"))
    end
    if length(basis) != getdim(F)
        throw(DimensionMismatch("basis has incorrect length"))
    end
    isbasic = zeros(Int64, n)
    isbasic[basis] .= 1
    if sum(isbasic) != m
        throw(ErrorException("duplicate index in basis"))
    end
    cbasis = basis .- 1
    Ap = getcolptr(A) .- 1
    Ai = rowvals(A) .- 1
    Ax = nonzeros(A)            # don't need a copy
    p_nupdate = Ref{Int64}(0)
    retcode = ccall((:basiclu_obj_maxvolume, libbasiclu), Int64,
                    (Ptr{basiclu_object_type}, Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Float64},
                     Ptr{Int64}, Ptr{Int64}, Float64, Ptr{Int64}),
                    pointer_from_objref(F), n, Ap, Ai, Ax, cbasis, isbasic, volumetol, p_nupdate)
    basis[:] = cbasis .+ 1
    checkretcode("basiclu_obj_maxvolume", retcode)
    return p_nupdate[]
end

"""
    basis, F = maxvolbasis(A::SparseMatrixCSC{Float64, Int64}; lindeptol::Float64=1e-8,
                           volumetol::Float64=2.0, maxpass::Int64=2, verbose::Bool=true) -> Vector{Int64}, LUFactor

Find a set of column indices for the matrix `AI = [A I]` such that `AI[:,basis]`
is square and nonsingular and the number of slack columns in the basis is
minimum (this is the row rank deficiency of `A`). Return the vector of column
indices of `AI` which form the basis matrix and a `LUFactor` which holds a
factorization of the basis matrix.

Method: Scale the slack columns of `AI` by `lindeptol` and try to find a maximum
volume basis for this matrix by making at most `maxpass` calls to
[`maxvolume`](@ref). If `verbose` is true, then print the number of basis
updates after each call.
"""
function maxvolbasis(A::SparseMatrixCSC{Float64, Int64}; lindeptol::Float64=1e-8,
                     volumetol::Float64=2.0, maxpass::Int64=2,
                     verbose::Bool=true)
    m, n = size(A)
    rowmax = maximum(abs.(A), dims=2)[:]
    rowmax = max.(rowmax, 1.)
    AI = [A spdiagm(0 => lindeptol.*rowmax)]
    basis = collect(n+1:n+m)
    F = LUFactor(m)
    for pass = 1:maxpass
        nupdate = maxvolume(F, AI, basis, volumetol)
        if verbose
            @printf("pass %d: %d updates\n", pass, nupdate)
        end
        if nupdate == 0 break; end
    end
    return basis, F
end

function checktrans(trans::Char)
    if !(trans == 'N' || trans == 'T')
        throw(ArgumentError("trans argument must be 'N' (no transpose) or 'T' (transpose), got '$trans'"))
    end
    trans
end

function checkretcode(funcname::String, retcode::Int64)
    if retcode == BASICLU_ERROR_out_of_memory
        throw(OutOfMemoryError())
    elseif retcode == BASICLU_WARNING_singular_matrix
        msg = @sprintf("%s() encountered singular matrix", funcname)
        throw(ErrorException(msg))
    elseif retcode == BASICLU_ERROR_singular_update
        msg = @sprintf("%s() encountered singular update", funcname)
        throw(ErrorException(msg))
    elseif retcode != BASICLU_OK
        msg = @sprintf("%s() failed with return code %d", funcname, retcode)
        throw(ErrorException(msg))
    end
end

function gathersol(F::LUFactor)
    dim = getdim(F)
    nzind = Vector{Int64}(undef, F.nzlhs)
    nzval = Vector{Float64}(undef, F.nzlhs)
    lhs = SparseVector{Float64, Int64}(dim, nzind, nzval)
    gathersol!(F, lhs)
end

function gathersol!(F::LUFactor, lhs::SparseVector{Float64, Int64})
    dim = getdim(F)
    @assert length(lhs) == dim
    if nnz(lhs) != F.nzlhs
        resize!(lhs.nzind, F.nzlhs)
        resize!(lhs.nzval, F.nzlhs)
    end
    for i = 1:F.nzlhs
        idx = unsafe_load(F.ilhs, i) + 1
        @assert idx >= 1 && idx <= dim
        lhs.nzind[i] = idx
        lhs.nzval[i] = unsafe_load(F.lhs, idx)
    end
    p = sortperm(lhs.nzind)
    lhs.nzind[:] = lhs.nzind[p]
    lhs.nzval[:] = lhs.nzval[p]
    lhs
end

end
