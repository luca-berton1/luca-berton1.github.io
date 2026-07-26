---
layout: post
title: "Implementing Zippel's Algorithm: Architecture of the Monic Case"
---

In the previous post, I described the interpolation machinery behind the inductive step of Zippel's algorithm.

The central idea was to combine two reconstruction procedures:

- sparse interpolation, used to reconstruct one image of the GCD from a known monomial structure;
- dense Newton interpolation, used to restore a variable that had previously been evaluated.

In this post, I will explain how these components are organized into a complete implementation of the monic case.

The implementation is part of [SymPy pull request #29824](https://github.com/sympy/sympy/pull/29824). Rather than following the code line by line, I will focus on its general architecture and on the responsibilities of its three main functions:

- `modgcd_multivariate`;
- `sparse_gcd`;
- `zippel_interp`.

Each function works at a different level of the algorithm, and each one performs a different kind of reconstruction.

## The overall architecture

The complete computation starts with two polynomials over the integers:

$$
A,B\in\mathbb{Z}[x_1,\ldots,x_n].
$$

Zippel's algorithm is used after reducing these polynomials modulo a prime number. The computation over a fixed finite field is recursive: one variable is evaluated, a GCD in fewer variables is computed, and the eliminated variable is then reconstructed through interpolation.

The resulting GCD modulo one prime is still not enough to recover the integer coefficients. The same computation is therefore repeated modulo different primes, and the modular images are combined using the Chinese remainder theorem.

The implementation follows this hierarchy:

```text
modgcd_multivariate
    |
    | reduce the integer inputs modulo a prime
    v
sparse_gcd
    |
    | evaluate the last variable and recurse
    v
GCD in one fewer variable
    |
    | use its monomial structure as a skeleton
    v
zippel_interp
    |
    | reconstruct a new GCD image
    v
sparse_gcd
    |
    | restore the evaluated variable with Newton interpolation
    v
GCD modulo the current prime
    |
    | combine several modular images
    v
modgcd_multivariate
    |
    v
GCD over the integers
```

There are therefore three nested reconstruction levels:

1. `zippel_interp` reconstructs one sparse GCD image;
2. `sparse_gcd` reconstructs one evaluated variable;
3. `modgcd_multivariate` reconstructs the integer coefficients.

## `modgcd_multivariate`: reconstruction over the integers

The outermost function is `modgcd_multivariate`.

Its inputs are multivariate polynomials over the integers. Its main responsibility is to move the expensive part of the computation to finite fields and then reconstruct the result over $\mathbb{Z}$.

Before starting the modular loop, the function removes the integer contents of the input polynomials. The GCD of those contents is stored and restored at the end, while the modular algorithm works with primitive polynomials.

The function also identifies primes that should not be used. A prime is unsuitable if reducing modulo it makes an important leading coefficient vanish, because this can change the degrees or the structure of the input polynomials.

For every suitable prime $p$, the inputs are mapped to

$$
\mathbb{F}_p[x_1,\ldots,x_n],
$$

and `sparse_gcd` is called to compute their GCD over the finite field.

Once a modular GCD image has been obtained, it is combined with the images computed modulo the previous primes using the Chinese remainder theorem.

This process continues until the reconstructed polynomial stops changing. Stabilization alone is not considered sufficient: the candidate is accepted only after checking that it divides both original primitive inputs exactly.

The outer function therefore separates two roles:

- modular computation is used to construct a candidate efficiently;
- exact division over the integers is used to verify it.

## `sparse_gcd`: recursion on the number of variables

The function `sparse_gcd` is the recursive core of the implementation.

It works over a fixed finite field:

$$
A,B\in\mathbb{F}_p[x_1,\ldots,x_n].
$$

Its recursion parameter is the number of variables.

When only one variable remains, the problem is univariate and the function calls the existing finite-field GCD routine. This is the base case of the recursion.

When more variables are present, the function first separates the appropriate contents and primitive parts of the inputs. The GCD of the contents is computed independently and multiplied back into the result at the end.

The last variable is then evaluated at a randomly chosen point:

$$
x_n=t_0.
$$

This gives two polynomials in one fewer variable:

$$
A_{t_0},B_{t_0}
\in
\mathbb{F}_p[x_1,\ldots,x_{n-1}].
$$

The function calls itself recursively to compute

$$
G_{t_0}
=
\gcd(A_{t_0},B_{t_0}).
$$

Assuming that the evaluation point preserves the relevant structure, this first image provides the monomial support expected in the full GCD. It is therefore reorganized into a skeletal representation using `skeleton_sorter`.

At this stage, the algorithm knows the expected monomials of the GCD in the first $n-1$ variables, but it does not yet know how their coefficients depend on $x_n$.

To determine this dependence, `sparse_gcd` evaluates $x_n$ at further points:

$$
t_1,t_2,\ldots
$$

For every new point, it calls `zippel_interp` to reconstruct a new GCD image having the same assumed skeleton.

These images provide the values of the unknown coefficient polynomials at the chosen evaluation points. Incremental Newton interpolation is then used to recover those coefficients as polynomials in $x_n$.

Once the Newton interpolations stop changing, the variable $x_n$ has been restored and a complete GCD candidate modulo $p$ can be constructed.

The candidate is then tested by exact division over the finite field. If it divides both inputs, the recursive computation terminates. Otherwise, the algorithm continues with additional evaluation points.

## `zippel_interp`: reconstruction of one sparse image

The function `zippel_interp` works one level below `sparse_gcd`.

At the moment it is called, the outermost variable of the current recursive step has already been fixed at some value. The purpose of `zippel_interp` is therefore not to restore that variable, but to reconstruct one GCD image in the remaining variables.

Suppose its inputs belong to

$$
\mathbb{F}_p[x_1,\ldots,x_k].
$$

The first variable $x_1$ is kept symbolic. The remaining variables

$$
x_2,\ldots,x_k
$$

are evaluated at successive powers of a randomly selected tuple.

After these evaluations, the inputs become univariate polynomials in $x_1$. Their GCDs can therefore be computed using a univariate finite-field algorithm.

The coefficients of these univariate GCDs provide the known values in the sparse interpolation systems described in the previous post.

The assumed monomial support comes from the skeletal GCD. Consequently, the unknowns are only the scalar coefficients attached to those known monomials.

By evaluating the same monomials at successive powers of one tuple, the associated systems acquire a Vandermonde structure. The functions `lag_basis` and `vandermonde_interp` are then used to solve them efficiently.

After all the coefficient groups have been reconstructed, `zippel_interp` returns one complete GCD image for the current value chosen by `sparse_gcd`.

The division of responsibilities is therefore important:

- `zippel_interp` reconstructs one image from univariate GCDs;
- `sparse_gcd` collects several such images and restores the eliminated variable.

## Where normalization enters the algorithm

Several GCDs are computed during the recursion, and a polynomial GCD over a field is determined only up to multiplication by a nonzero scalar.

For interpolation to work, the different images must be represented using a consistent normalization.

The implementation therefore normalizes recursive and interpolated images before combining them through dense Newton interpolation. In particular, the normalization uses information obtained from the GCD of suitable leading coefficients.

This step is needed to make the evaluations correspond to values of the same polynomial coefficient functions.

However, this should not be confused with the additional normalization problem that arises during sparse interpolation in the general non-monic case.

In that case, the univariate GCD images used to build the Vandermonde systems can differ by independent unknown scaling factors. Those factors have to be recovered together with some of the coefficients of the skeletal GCD.

I will return both to the precise normalization used during dense interpolation and to the more difficult non-monic normalization problem in a later post.

## Random evaluation and verification

Zippel's algorithm is probabilistic because it relies on random choices:

- values assigned to the recursively eliminated variable;
- tuples used during sparse interpolation;
- primes used by the outer modular algorithm.

Some choices are unsuitable.

For example, an evaluation can make a leading coefficient vanish, cause two different monomials to become indistinguishable or produce a GCD image whose support does not match the assumed skeleton.

The implementation detects these situations and discards the corresponding attempt.

Randomness is therefore used to find convenient evaluation points, but it does not replace correctness checks.

Before accepting a result:

- `sparse_gcd` verifies divisibility modulo $p$;
- `modgcd_multivariate` verifies divisibility over the integers.

This gives the algorithm a useful separation between candidate construction and final verification.

## Current status

The monic implementation now combines the complete recursive pipeline:

- reduction modulo suitable primes;
- recursion on the number of variables;
- extraction of a skeletal GCD;
- reconstruction of new images through sparse interpolation;
- restoration of evaluated variables through Newton interpolation;
- Chinese remainder reconstruction;
- exact divisibility checks.

The implementation is contained in [pull request #29824](https://github.com/sympy/sympy/pull/29824).

The pull request is still marked as a draft because some comments, documentation and style details need to be improved, but the complete algorithm is implemented and working.

The next post on this pull request will focus on the non-monic case, where the sparse interpolation stage must also determine the unknown scaling factors relating the different univariate GCD images.

## References

1. Suling Yang, *Computing the Greatest Common Divisor of Multivariate Polynomials over Finite Fields*.
2. Keith O. Geddes, Stephen R. Czapor and George Labahn, *Algorithms for Computer Algebra*.
3. [SymPy pull request #29824](https://github.com/sympy/sympy/pull/29824).
