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

and regard the inputs as polynomials in $X$ whose coefficients belong to
$\mathbb{F}_p[x_n]$. Thus, we can write

$$
A=\sum_{\alpha}a_{\alpha}(x_{n})X^{\alpha},
\qquad
B=\sum_{\alpha}b_{\alpha}(x_{n})X^{\alpha}.
$$

The corresponding contents are

$$
a(x_{n})=\operatorname{cont}_{X}(A),
\qquad
b(x_{n})=\operatorname{cont}_{X}(B),
$$

that is, the GCDs in $\mathbb{F}_p[x_n]$

of the coefficient
polynomials $a_{\alpha}(x_{n})$ and $b_{\alpha}(x_{n})$.

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

This is the normalization required for dense interpolation. Every recursive or sparsely reconstructed GCD is multiplied by $q(x_n)$ evaluated in the current evaluation point, and therefore when the last variable is reconstructed, the gcd is inflated by $q(x_n)$. 

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
## `zippel_interp`: reconstructing one sparse image

The function `zippel_interp` operates one level below `sparse_gcd`.

At the moment it is called, the last variable of the current recursive step has already been assigned a fixed value. Its purpose is to reconstruct the corresponding GCD image in the remaining variables:

$$
G\in\mathbb{F}_{p}[x_{1},\ldots,x_{k}].
$$

The function does not reconstruct the eliminated variable. That task is performed later by `sparse_gcd` through dense Newton interpolation.

Instead, `zippel_interp` keeps the first variable $x_{1}$ symbolic, evaluates the other variables at powers of a suitable tuple, computes a collection of univariate GCDs and uses them to recover the coefficients associated with the known skeletal support.

When only $x_{1}$ remains, no sparse interpolation is necessary. The function computes the univariate GCD directly and extracts the coefficients corresponding to the degrees present in the skeleton.

### Choosing a suitable evaluation tuple

Assume that the skeletal GCD has been organized as

$$
G
=
\sum_{j}C_{j}(x_{2},\ldots,x_{k})x_{1}^{j},
$$

where each sparse coefficient has the form

$$
C_{j}
=
\sum_{s=1}^{t_{j}}M_{s,j}W_{s,j}.
$$

The monomials $W_{s,j}$ are known from the skeleton, while the scalar coefficients $M_{s,j}$ must be reconstructed.

The function begins by choosing a random tuple

$$
b=(b_{2},\ldots,b_{k})
\in\mathbb{F}_{p}^{k-1}.
$$

Not every tuple is suitable. Before performing the univariate GCD computations, the implementation checks that the chosen tuple does not immediately destroy important information.

First, it verifies that the leading coefficient of the first input with respect to $x_{1}$ does not vanish at the powers of the tuple that will be used. If this leading coefficient vanished, the degree in $x_{1}$ of the evaluated polynomial could decrease, and the resulting univariate GCD would no longer represent the expected image.

The function also evaluates the monomials belonging to every coefficient group of the skeleton. For each $W_{s,j}$, it computes

$$
\beta_{s,j}=W_{s,j}(b).
$$

Within the same coefficient $C_{j}$, these values must be distinct. If two different skeletal monomials produce the same value,

$$
W_{r,j}(b)=W_{s,j}(b),
$$

the corresponding columns of the interpolation system cannot be distinguished. The tuple is therefore discarded and a new one is selected.

These checks do not prove that every aspect of the evaluation is good, but they eliminate the choices that would immediately make the interpolation systems unusable.

### Building the Vandermonde bases

The mathematical construction of the interpolation systems was described in the previous post. Here, the important point is how the implementation organizes the information required by those systems.

For every power $x_{1}^{j}$ occurring in the skeletal GCD, the function computes the values

$$
\beta_{1,j},\ldots,\beta_{t_{j},j}
$$

of its skeletal monomials at the base tuple.

These values are stored in `all_vand_basis`, a list of lists. Each inner list corresponds to one coefficient $C_{j}$ and contains the Vandermonde values associated with its monomials:

$$
\texttt{all\_vand\_basis}
=
\left[
 [\beta_{1,j_{1}},\ldots,\beta_{t_{j_{1}},j_{1}}],
 \ldots,
 [\beta_{1,j_{m}},\ldots,\beta_{t_{j_{m}},j_{m}}]
\right].
$$

The order of the outer list is the order of the coefficient groups in the skeletal representation. Inside each group, the values follow the order of the corresponding monomials.

Preserving this order is essential: after solving a system, the reconstructed scalar coefficients must be associated with exactly the same monomials from which its Vandermonde basis was constructed.

### Flattening the input polynomials

The next step is to evaluate the input polynomials at successive powers of the same tuple:

$$
b^{r}
=
(b_{2}^{r},\ldots,b_{k}^{r}).
$$

Performing every multivariate evaluation from scratch would repeatedly require computing the same monomial expressions. The current implementation therefore converts the inputs into an evaluation-oriented representation.

Consider the first input written as a polynomial in $x_{1}$:

$$
A
=
\sum_{d}A_{d}(x_{2},\ldots,x_{k})x_{1}^{d}.
$$

Each coefficient $A_{d}$ is a sparse polynomial. For every monomial $W$ appearing in $A_{d}$, the function first computes its value at the base tuple,

$$
\beta=W(b),
$$

and stores the original scalar coefficient under this value.

The flattened representation is therefore a list indexed by the degree $d$ in $x_{1}$. Each entry is a dictionary representing a linear combination of powers of previously computed monomial values.

Conceptually, if

$$
A_{d}
=
c_{1}W_{1}+\cdots+c_{q}W_{q},
$$

the corresponding dictionary stores the pairs

$$
W_{1}(b)\longmapsto c_{1},
\quad\ldots,\quad
W_{q}(b)\longmapsto c_{q}.
$$

If different monomials have the same evaluated value, their scalar coefficients are combined in the same dictionary entry.

The same transformation is applied to $B$.

This representation makes evaluations at powers of the tuple particularly simple. Since

$$
W(b^{r})=W(b)^{r},
$$

the value of $A_{d}$ at $b^{r}$ can be computed as

$$
A_{d}(b^{r})
=
\sum_{\beta}c_{d,\beta}\beta^{r}.
$$

Thus, after the initial flattening, every new evaluation is reduced to exponentiation and a linear combination. The algorithm no longer needs to traverse the exponent vector of every multivariate monomial at every interpolation point.

The implementation starts with the first power of the tuple rather than
with its zeroth power. Mathematically, using the zeroth power would still
give a valid Vandermonde system, whose first row would consist entirely of
ones. However, it would force every interpolation attempt to use the fixed
evaluation point $(1,\ldots,1)$. If this point caused a leading coefficient
to vanish or otherwise produced an unlucky specialization, choosing a new
random base tuple would not avoid the problem. Starting from $b^1=b$ ensures
that every evaluation point depends on the randomly chosen tuple and can
therefore be replaced when it is unsuitable.  

### Computing the right-hand sides

For each required power $b^{r}$, the flattened representations are used to construct two univariate polynomials in $x_{1}$:

$$
A_{r}(x_{1})
=
A(x_{1},b_{2}^{r},\ldots,b_{k}^{r}),
$$

and

$$
B_{r}(x_{1})
=
B(x_{1},b_{2}^{r},\ldots,b_{k}^{r}).
$$

Their univariate GCD is then computed:

$$
H_{r}(x_{1})
=
\gcd(A_{r},B_{r}).
$$

Write this image as

$$
H_{r}(x_{1})
=
\sum_{j}h_{r,j}x_{1}^{j}.
$$

For each degree $j$ appearing in the skeleton, the coefficient $h_{r,j}$ is appended to the collection of values associated with $C_{j}$.

The implementation stores these values in `eval_points`. This is a dictionary whose keys are the degrees $j$ in $x_{1}$ and whose values are lists:

$$
\texttt{eval\_points}[j]
=
[h_{1,j},h_{2,j},\ldots].
$$

Conceptually, this is another list-of-lists structure. Its ordering matches the ordering used for `all_vand_basis`:

- the $i$-th Vandermonde basis belongs to the $i$-th coefficient group in the skeleton;
- the corresponding list of evaluated coefficients contains the right-hand side of that same interpolation system.

The two collections therefore remain aligned throughout the function.

While collecting these values, the implementation also checks the support in the main variable. If a univariate GCD contains a nonzero power of $x_{1}$ that is absent from the skeletal GCD, the assumed skeleton is incompatible with the new image. In that case, the function returns `None`, allowing the surrounding algorithm to discard the current attempt.

### Solving the Vandermonde systems

After enough univariate GCD images have been computed, the function has all the data required for sparse interpolation.

For a fixed coefficient

$$
C_{j}
=
\sum_{s=1}^{t_{j}}M_{s,j}W_{s,j},
$$

the collected equations are

$$
\sum_{s=1}^{t_{j}}
M_{s,j}\beta_{s,j}^{r}
=
h_{r,j},
\qquad
r=1,\ldots,t_{j}.
$$

The values in `all_vand_basis` determine the matrix, while the corresponding values in `eval_points` determine its right-hand side.

For every coefficient group, `lag_basis` first constructs the interpolation basis associated with

$$
\beta_{1,j},\ldots,\beta_{t_{j},j}.
$$

The function `vandermonde_interp` then combines this basis with

$$
h_{1,j},\ldots,h_{t_{j},j}
$$

and solves the transposed Vandermonde system in quadratic time.

Because the powers of the tuple begin at $r=1$, the system can be rewritten as

$$
\sum_{s=1}^{t_{j}}
\left(M_{s,j}\beta_{s,j}\right)
\beta_{s,j}^{r-1}
=
h_{r,j}.
$$

The Vandermonde solver therefore initially reconstructs the quantities

$$
M_{s,j}\beta_{s,j}.
$$

The final multiplication by $\beta_{s,j}^{-1}$ recovers the desired coefficients

$$
M_{s,j}.
$$

The same operation is performed independently for every power of $x_{1}$ in the skeleton. The recovered scalar coefficients are then associated with their original monomials, using the common ordering of the skeleton, the Vandermonde bases and the right-hand sides.

At the end of this process, `zippel_interp` returns a dictionary representing one complete GCD image for the current value of the variable previously fixed by `sparse_gcd`.


## Why does the not monic case need a separate sparse interpolation procedure?

During sparse interpolation, the univariate GCDs used to build the Vandermonde systems can differ by independent unknown scaling factors.
This happens because the unvariate GCDs are always monic (since we have coefficients in a field). 

I will return to this normalization problem in a later post.

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

The next post on this pull request will focus on the non-monic case.

## References

1. Suling Yang, *Computing the Greatest Common Divisor of Multivariate Polynomials over Finite Fields*.
2. Keith O. Geddes, Stephen R. Czapor and George Labahn, *Algorithms for Computer Algebra*.
3. [SymPy pull request #29824](https://github.com/sympy/sympy/pull/29824).
