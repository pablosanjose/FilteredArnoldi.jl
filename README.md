# FilteredMatrices [![Build Status](https://github.com/pablosanjose/FilteredMatrices.jl/workflows/CI/badge.svg)](https://github.com/pablosanjose/FilteredMatrices.jl/actions) [![Coverage](https://codecov.io/gh/pablosanjose/FilteredMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pablosanjose/FilteredMatrices.jl)

The FilteredMatrices.jl package aims at computing interior eigenvalues and eigenvectors ("eigenpairs") of large (typically sparse) matrices in a way alternative to more standard shift-invert methods.

To compute a collection of interior eigenpairs of a matrix H close to a given value σ, traditional shift-invert methods use the Arnoldi (or Lanczos) method to compute extremal eigenvalues of `inv(H-σ*I)`. Since the Arnoldi method only requires to compute matrix-vector products, there is no need to compute `inv(H-σ*I)` itself (which would be a *bad* idea for large sparse matrices). Instead, the product `w = inv(H-σ*I) * v` for a vector `v` is obtained by solving the linear equations `(H-σ*I) * w = v`, i.e. `w = v \ (H-σ*I)`. A factorization of `H-σ*I` can be done up-front so that repeated matrix-vector multiplications become very fast.

In this package we take a different approach that requires no factorizations. Instead of applying the Arnoldi algorithm on `inv(H-σ*I)`, we apply it to `δ_K(H-σ)`, where `δ_K(x)` is an approximation to a Delta function using an order-K expansion in orthogonal Chebyshev polynomials. Since `δ_K(ε-σ)` is peaked close to `ε == σ` as `K` increases, exterior eigenvalues of `δ_K(H-σ)` are interior eigenvalues of `H`. Moreover, since the expansion `δ_K(H-σ)` is a polynomial of `H`, computing the matrix-vector product `δ_K(H-σ) * v` requires only a series of `H * w` products.

The package defines the method `delta(H, σ; order = K)` which produces `LinearMap` representing `δ_K(H-σ)`. If `K` is not provided, a reasonable value is computed automatically. This `LinearMap` can then be fed into any iterative extremal eigensolver, such as e.g. KrylovKit.jl or IterativeSolvers.jl, to obtain the interior eigenvalues closest to the maximum of `δ_K(ε-σ)`. The obtained eigenvalues will be `νᵢ = δ_K(εᵢ-σ)`. To recover `εᵢ` use the associated eigenvector `ψᵢ` and `εᵢ = ψᵢ' * H * ψᵢ / dot(ψᵢ, ψᵢ)`.

Since the maximum of `δ_K(ε-σ)` is only at `ε = σ` in the large-`K` limit, the method is not guaranteed to return the eigenpairs in strict order from `σ` unless `K` is large enough.

# Usage example
In this example we compute the zero energy eigenstates of an SSH model using KrylovKit.jl and FilteredMatrices.jl. The system is constructed with Quantica.jl.
```julia
julia> using Quantica, KrylovKit, FilteredMatrices

julia> ssh = LatticePresets.linear() |> hamiltonian(hopping(1)) |> unitcell(1000, modifiers = @hopping!((t,r,dr) -> t + 0.2*mod(r[1],2))) |> unitcell
Hamiltonian{<:Lattice} : Hamiltonian on a 0D Lattice in 1D space
  Bloch harmonics  : 1 (SparseMatrixCSC, sparse)
  Harmonic size    : 1000 × 1000
  Orbitals         : ((:a,),)
  Element type     : scalar (Complex{Float64})
  Onsites          : 0
  Hoppings         : 1998
  Coordination     : 1.998

julia> m = bloch(ssh);

julia> emin, emax = eigsolve(m, 1, :SR)[1][1], eigsolve(m, 1, :LR)[1][1]
(-2.399988200192531, 2.3999882601634264)

julia> o = order_estimate(0, 0.05, (emin, emax))
1207

julia> d = delta(m, 0, (emin, emax), order = o);

julia> states = eigsolve(x -> d*x, size(m,1), 2, :LR)[2];

julia> using Plots; plot(real(states))
```
![image](https://user-images.githubusercontent.com/4310809/92922536-4b7b3b00-f436-11ea-8d13-f474fc89e142.png)
