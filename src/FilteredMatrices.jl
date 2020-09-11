module FilteredMatrices

using LinearAlgebra, LinearMaps, UnPack

export delta, Kernels

function delta(H::AbstractMatrix{T}, ε; range = bandrange(H), order = defaultorder(H), kernel = (n, order) -> 1) where {T<:Number}
    M, N = size(H)
    M == N || throw(ArgumentError("Only square matrices are supported"))
    bracket = bandbracket(range)
    v = Vector{T}(undef, N)
    v´ = Vector{T}(undef, N)
    args = (H = H, ε = ε, bracket = bracket, order = order, v = v, v´ = v´, kernel = kernel)
    LinearMap{T}((vdst, v0) -> delta_mul!(vdst, v0, args), N, N; ismutating = true, ishermitian = true)
end

bandrange(H) = (-1, 1)

bandbracket((εmin, εmax)) = (εmax + εmin)/2, abs(εmax - εmin)

defaultorder(H) = 10 #round(Int, sqrt(size(H,1)))

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

Lorentz(λ = 4) = (n, order) -> sinh(λ*(1-n/order))/sinh(λ)

Lanczos(M = 3) = (n, order) -> sinc(n/order)^M

end # module

end # module