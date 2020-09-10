# FilteredMatrices [![Build Status](https://github.com/pablosanjose/FilteredMatrices.jl/workflows/CI/badge.svg)](https://github.com/pablosanjose/FilteredMatrices.jl/actions) [![Coverage](https://codecov.io/gh/pablosanjose/FilteredMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pablosanjose/FilteredMatrices.jl)

The FilteredMatrices.jl package aims at computing interior eigenvalues and eigenvectors ("eigenpairs") of large (typically sparse) matrices in a way alternative to more standard shift-invert methods.

To compute a collection of interior eigenpairs of a matrix H close to a given value σ, traditional shift-invert methods use the Arnoldi (or Lanczos) method to compute extremal eigenvalues of `inv(H-σ*I)`. Since the Arnoldi method only requires to compute matrix-vector products, there is no need to compute `inv(H-σ*I)` itself (which would be a *bad* idea for large sparse matrices). Instead, the product `w = inv(H-σ*I) * v` for a vector `v` is obtained by solving the linear equations `(H-σ*I) * w = v`, i.e. `w = v \ (H-σ*I)`. A factorization of `H-σ*I` can be done up-front so that repeated matrix-vector multiplications become very fast. This factorization, however, eventually becomes a problem for very large matrices, as it scales as O(N^3) with standard algorithms.

In this package we take a different approach that requires no factorizations. Instead of applying the Arnoldi algorithm on `inv(H-σ*I)`, we apply it to `δ_K(H-σ)`, where `δ_K(x)` is an approximation to a Delta function using an order-K expansion in orthogonal Chebyshev polynomials. Since `δ_K(ε-σ)` is peaked close to `ε == σ` as `K` increases, exterior eigenvalues of `δ_K(H-σ)` are interior eigenvalues of `H`. Moreover, since the expansion `δ_K(H-σ)` is a polynomial of `H`, computing the matrix-vector product `δ_K(H-σ) * v` requires only a series of `H * w` products.

The package defines the method `delta(H, σ; order = K)` which produces `LinearMap` representing `δ_K(H-σ)`. If `K` is not provided, a reasonable value is computed automatically. This `LinearMap` can then be fed into any iterative extremal eigensolver, such as e.g. KrylovKit.jl or IterativeSolvers.jl, to obtain the interior eigenvalues closest to the maximum of `δ_K(ε-σ)`. The obtained eigenvalues will be `νᵢ = δ_K(εᵢ-σ)`. To recover `εᵢ` use the associated eigenvector `ψᵢ` and `εᵢ = ψᵢ' * H * ψᵢ / dot(ψᵢ, ψᵢ)`.

Since the maximum of `δ_K(ε-σ)` is only at `ε = σ` in the large-`K` limit, the method is not guaranteed to return the eigenpairs in strict order from `σ` unless `K` is large enough.

# Usage example

```
```