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
A,B\in\mathbb{F}_p[x_1,\ldots,x_n],
$$

and its recursion parameter is the number of variables.

When only one variable remains, the problem is univariate and the function calls the existing finite-field GCD routine. This is the base case of the recursion.

When more variables are present, the function first extracts contents and primitive parts with respect to the first $n-1$ variables. This is different from extracting the scalar content of a polynomial.

Let

$$
X=(x_1,\ldots,x_{n-1}),
$$

and regard the inputs as polynomials in $X$ whose coefficients belong to $\mathbb{F}_p[x_n]$. Thus, we can write

$$
A=\sum_{\alpha}a_{\alpha}(x_n)X^{\alpha},
\qquad
B=\sum_{\alpha}b_{\alpha}(x_n)X^{\alpha}.
$$

The corresponding contents are

$$
a(x_n)=\operatorname{cont}_X(A),
\qquad
b(x_n)=\operatorname{cont}_X(B),
$$

that is, the GCDs in $\mathbb{F}_p[x_n]$ of the coefficient polynomials $a_{\alpha}(x_n)$ and $b_{\alpha}(x_n)$.

The helper `_primitive` separates these contents from the primitive parts:

$$
A=a(x_n)\overline{A},
\qquad
B=b(x_n)\overline{B}.
$$

The function computes

$$
c(x_n)=\gcd(a(x_n),b(x_n))
$$

independently. This factor is restored after the GCD of the primitive parts has been reconstructed.

### Leading coefficients and consistent normalization

The function also computes

$$
g(x_n)
=
\gcd\left(
\operatorname{LC}_X(\overline{A}),
\operatorname{LC}_X(\overline{B})
\right),
$$

where the leading coefficients are taken with respect to the variables in $X$. Consequently, they are polynomials in the remaining variable $x_n$.

Let

$$
H=\gcd(\overline{A},\overline{B})
$$

and denote its leading coefficient with respect to $X$ by

$$
\ell(x_n)=\operatorname{LC}_X(H).
$$

Since $H$ divides both primitive inputs, its leading coefficient must divide the GCD of their leading coefficients:

$$
\ell(x_n)\mid g(x_n).
$$

We can therefore write

$$
g(x_n)=q(x_n)\ell(x_n)
$$

for some polynomial $q(x_n)\in\mathbb{F}_p[x_n]$.

This relation is used to normalize the different GCD images before dense interpolation.

The last variable is first evaluated at a suitable random point:

$$
x_n=t_0.
$$

This produces two polynomials in one fewer variable:

$$
\overline{A}_{t_0},
\overline{B}_{t_0}
\in
\mathbb{F}_p[x_1,\ldots,x_{n-1}].
$$

The function calls itself recursively to compute their GCD. A GCD over a field is determined only up to multiplication by a nonzero scalar, so the result returned by the recursive computation can be written as

$$
\widetilde{H}_{t_0}
=
u_0H_{t_0},
$$

where $u_0\in\mathbb{F}_p^{\times}$ and

$$
H_{t_0}=H(x_1,\ldots,x_{n-1},t_0).
$$

Assuming that the evaluation preserves the relevant degrees, making this image monic gives

$$
\operatorname{monic}(\widetilde{H}_{t_0})
=
\frac{H_{t_0}}{\ell(t_0)}.
$$

The implementation then multiplies it by the evaluated polynomial $g(t_0)$:

$$
g(t_0)\operatorname{monic}(\widetilde{H}_{t_0})
=
\frac{g(t_0)}{\ell(t_0)}H_{t_0}
=
q(t_0)H_{t_0}.
$$

The same operation is performed for every subsequent evaluation point $t_i$. Therefore, the normalized images are not multiplied by unrelated arbitrary factors: they are all evaluations of the same polynomial

$$
\widehat{H}=q(x_n)H.
$$

This is the normalization required for dense interpolation. Without it, every recursive or sparsely reconstructed GCD image could have a different scalar factor, and their coefficients would not be values of a common family of polynomials in $x_n$.

This normalization should not be confused with the additional normalization problem arising inside sparse interpolation in the general non-monic case. I will discuss that separate issue in a later post.

### Constructing the skeletal GCD

After the first recursive image has been normalized, it is passed to `skeleton_sorter`.

This function records two kinds of information.

First, it extracts and organizes the monomial support of the image. This support is used as the assumed skeleton when later copies of the GCD are reconstructed by `zippel_interp`.

Second, it extracts the scalar coefficient attached to each monomial of the first normalized image.

Suppose that the normalized polynomial being reconstructed has the form

$$
\widehat{H}
=
\sum_{\alpha}
C_{\alpha}(x_n)X^{\alpha},
$$

where the monomials $X^{\alpha}$ are given by the skeleton and the coefficients

$$
C_{\alpha}(x_n)\in\mathbb{F}_p[x_n]
$$

are still unknown.

The first normalized image provides the values

$$
C_{\alpha}(t_0).
$$

For each monomial $X^{\alpha}$, `skeleton_sorter` initializes a list containing this first value:

$$
h_{\alpha}=[C_{\alpha}(t_0)].
$$

At this stage, the list is also the initial Newton representation of the unknown coefficient polynomial: an interpolation through a single point is simply the constant polynomial $C_{\alpha}(t_0)$.

The evaluation point $t_0$ is stored separately in the list of Newton interpolation points.

### Reconstructing further GCD images

The function then chooses additional values

$$
t_1,t_2,\ldots
$$

for the last variable.

For each new point $t_i$, the inputs are evaluated again and `zippel_interp` is called. Instead of recomputing the complete recursive GCD from scratch, `zippel_interp` uses the previously determined skeleton and sparse interpolation to reconstruct a new GCD image.

The returned image is again defined only up to a nonzero scalar. It is therefore monicized and multiplied by $g(t_i)$, exactly as for the first recursive image.

After this normalization, the coefficient of each monomial $X^{\alpha}$ is

$$
C_{\alpha}(t_i).
$$

The different images consequently provide a sequence of values

$$
C_{\alpha}(t_0),
C_{\alpha}(t_1),
C_{\alpha}(t_2),
\ldots
$$

for every coefficient polynomial appearing in the skeleton.

### Incremental Newton interpolation

Each coefficient polynomial $C_{\alpha}(x_n)$ is reconstructed independently using incremental Newton interpolation.

After the values at the first $r+1$ points have been processed, its Newton representation has the form

$$
P_r(x_n)
=
v_0
+
v_1(x_n-t_0)
+
v_2(x_n-t_0)(x_n-t_1)
+\cdots+
v_r\prod_{j=0}^{r-1}(x_n-t_j).
$$

For every monomial in the skeleton, the corresponding list initialized by `skeleton_sorter` is progressively extended:

$$
h_{\alpha}
=
[v_0,v_1,\ldots,v_r].
$$

When a new image at $t_{r+1}$ is obtained, its coefficient $C_{\alpha}(t_{r+1})$ is passed to `incremental_newton_interp`, together with the previous evaluation points and Newton coefficients.

The function calculates only the next coefficient $v_{r+1}$, rather than repeating the complete interpolation from the beginning.

Thus, the two main data structures used by `sparse_gcd` have complementary roles:

- the list of evaluation points stores
  $$
  t_0,t_1,\ldots,t_r;
  $$
- each list produced from the first skeletal GCD stores the Newton coefficients of one polynomial $C_{\alpha}(x_n)$.

The implementation treats a zero new Newton coefficient as an indication that the corresponding interpolation has stabilized. Coefficients that have stabilized are skipped during later updates.

Once every coefficient interpolation has stopped changing, the helper `from_newt_to_poly` converts each Newton representation into the usual dense univariate representation in $x_n$.

These reconstructed coefficient polynomials are then attached to their corresponding monomials $X^{\alpha}$, giving the complete candidate

$$
\widehat{H}
=
q(x_n)H.
$$

### Removing the extra content and verifying the result

The normalization factor $q(x_n)$ is common to every coefficient of $\widehat{H}$. It therefore appears as content with respect to the variables $X$.

After the dense interpolation has been completed, `_primitive` is called again. This removes the common factor $q(x_n)$ and recovers the primitive GCD $H$.

Finally, the candidate is tested by exact division over the finite field.

If it divides both primitive inputs, the previously computed content GCD $c(x_n)$ is restored and the result is returned. Otherwise, the interpolation is not accepted and the algorithm continues with further evaluation points.

The role of `sparse_gcd` can therefore be summarized as follows:

1. extract polynomial contents with respect to the last variable;
2. evaluate that variable and recursively compute the first GCD image;
3. normalize all images so that they are evaluations of one common polynomial;
4. use the first image to initialize the skeleton and the Newton representations;
5. obtain further images through `zippel_interp`;
6. update the coefficient polynomials by incremental Newton interpolation;
7. remove the extra normalization content;
8. verify the reconstructed GCD by exact division.
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
