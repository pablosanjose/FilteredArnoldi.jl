module FilteredMatrices

using LinearAlgebra, LinearMaps, UnPack

export delta, order_estimate, Kernels

"""
    delta(H::AbstractMatrix{<:Number}, ε::Number; range, order = 50, kernel = missing)

Construct a `LinearMap` that acts like matrix `δ_K(ε-H)`, where `δ_K` is a Chebyshev
expansion of a Delta function to `order = K`. See `order_estimate` for estimators based on
the desired energy resolution. The `range = (εmin, εmax)` of the spectrum of `H` needs to be
provided for convergence. The `kernel` is a function `k(n, order)`. Built-in kernels are:
    - missing                 : no kernel
    - Kernels.Jackson         : Jackson kernel
    - Kernels.Fejer           : Fejer kernel
    - Kernels.Lanczos(M = 3)  : Lanczos kernel `sinc(n/order)^M` with integer `M`
    - Kernels.Lorentz(λ = 4)  : Lorentz kernel `sinh(λ(1-n/order))/sinh(λ)` with real `λ`
"""
function delta(H::AbstractMatrix{T}, ε::Number; range, order = 50, kernel = (n, order) -> 1) where {T<:Number}
    M, N = size(H)
    M == N || throw(ArgumentError("Only square matrices are supported"))
    bracket = bandbracket(range)
    v = Vector{T}(undef, N)
    v´ = Vector{T}(undef, N)
    args = (H = H, ε = ε, bracket = bracket, order = order, v = v, v´ = v´, kernel = kernel)
    LinearMap{T}((vdst, v0) -> delta_mul!(vdst, v0, args), N, N; ismutating = true, ishermitian = true)
end

bandbracket((εmin, εmax)) = (εmax + εmin)/2, abs(εmax - εmin)

"""
    order_estimate(ε, Δε, range)

Estimate the Chebyshev expansion order to achieve a resolution `Δε` at energy `ε` within a
`range = (εmin, εmax)`. In terms of normalized energies `σ` and `Δσ` in a `(-1,1)` window,
the required order must be greater than 4π*sqrt(1-σ^2)/Δσ. The prefactor `4` is approximate
and can depend on the kernel used and the definition of resolution.
"""
function order_estimate(ε, Δε, range)
    (center, halfwidth) = bandbracket(range)
    σ = (ε - center)/halfwidth
    abs(σ) < 1 || throw(ArgumentError("Energy outside band range"))
    Δσ = Δε/halfwidth
    return ceil(Int, 4π*sqrt(1-σ^2)/Δσ)
end

function delta_mul!(vdst::AbstractVector{T}, v0, args) where {T}
    @unpack H, ε, bracket, order, v, v´, kernel = args

    # We use adjoint of CSC matrix to accelerate multiplication
    H´ = H'
    center, halfwidth = bracket
    α = 2/halfwidth
    β = 2center/halfwidth
    σ = (ε - center)/halfwidth
    ρ = 1/(pi*sqrt(1-σ^2))

    # t and t´ are used to build tn = ρ Tn(σ) recurseively
    t  = ρ
    t´ = σ * t

    # v and v´ are used to build vn = Tn(h)⋅v recursively, where h = (H-center)/halfwidth
    @. v = v0
    fill!(v´, zero(T))
    nextCheby!(v´, H´, v, α, β)

    # vdst accumulates delta(σ-H) * v ≈ ρ (v + 2∑ Tn(σ)Tn(h)⋅v) = t0 v0 + 2 ∑ tn vn
    @. vdst = kernel(0, order) * t * v + kernel(1, order) * 2t´ * v´

    for n in 2:order
        nextCheby!(v, H´, v´, α, β)
        t = nextCheby(t, σ, t´)
        @. vdst += kernel(n, order) * 2t * v
        v´, v = v, v´
        t´, t = t, t´
    end
    return vdst
end

function nextCheby!(v´, H´, v, α, β)
    mul!(v´, H´, v, α, -1)
    @. v´ -= β * v
    return v´
end

nextCheby(t´, σ, t) = 2σ * t - t´

module Kernels

Jackson(n, order) =
    ((order - n + 1) * cos(π * n / (order + 1)) + sin(π * n / (order + 1)) * cot(π / (order + 1))) / (order + 1)

Fejer(n, order) = 1 - n/order

Lorentz(λ::Real = 4) = (n, order) -> sinh(λ*(1-n/order))/sinh(λ)

Lanczos(M::Integer = 3) = (n, order) -> sinc(n/order)^M

end # module

end # module