---
layout: post
title: "Implementing Zippel's Algorithm: The Monic Case"
---

In the previous post, I described the interpolation machinery behind the inductive step of Zippel's algorithm. In particular, I explained how sparse interpolation reconstructs one evaluated image of the GCD, while incremental Newton interpolation restores the variable eliminated by the recursive step.

In this post, I will focus on how these mathematical components are combined into a working implementation for the monic case.

The implementation is part of [SymPy pull request #29824](https://github.com/sympy/sympy/pull/29824). The pull request contains both the monic and non-monic versions of the algorithm, but I will discuss only the monic path here. The normalization problem arising in the non-monic case deserves a separate treatment.

## The three levels of the implementation

The implementation is organized around three main functions:

- `modgcd_multivariate`;
- `sparse_gcd`;
- `zippel_interp`.

They operate at different levels of the algorithm.

`modgcd_multivariate` is the outer modular layer. It starts with polynomials over the integers, reduces them modulo different primes, calls the finite-field GCD algorithm and reconstructs the result over the integers using the Chinese remainder theorem.

`sparse_gcd` is the recursive core of Zippel's algorithm. It works over a fixed finite field and performs induction on the number of variables. It evaluates the last variable, recursively computes a GCD in one fewer variable and then reconstructs the eliminated variable using sparse and dense interpolation.

`zippel_interp` performs one sparse interpolation step. For a fixed value of the variable eliminated by `sparse_gcd`, it uses the known skeletal GCD and a collection of univariate GCD computations to reconstruct a new multivariate GCD image.

The overall flow can be summarized as follows:

```text
Polynomials over ZZ
        |
        | reduction modulo p
        v
modgcd_multivariate
        |
        | GCD over a finite field
        v
sparse_gcd
        |
        | evaluate the last variable
        | and recurse
        v
GCD in one fewer variable
        |
        | use it as the skeleton
        v
zippel_interp
        |
        | compute univariate GCD images
        | and perform sparse interpolation
        v
New GCD image
        |
        | incremental Newton interpolation
        v
GCD modulo p
        |
        | reconstruction over several primes
        v
GCD over ZZ
```

The three functions therefore correspond to three distinct reconstruction levels:

1. reconstruction over the integers from modular images;
2. reconstruction of an eliminated variable;
3. reconstruction of one sparse GCD image from univariate evaluations.

## `modgcd_multivariate`: the outer modular algorithm

The public entry point of the implementation is `modgcd_multivariate`.

Its purpose is to compute the GCD of two polynomials in

$$
\mathbb{Z}[x_1,\ldots,x_n]
$$

by performing most of the work over finite fields.

### Removing the integer content

The first relevant step is to separate the integer content from the primitive parts of the input polynomials:

```python
cf, f = f.primitive()
cg, g = g.primitive()
ch = ring.domain.gcd(cf, cg)

gamma = ring.domain.gcd(f.LC, g.LC)
```

The contents `cf` and `cg` are the GCDs of the integer coefficients of the two polynomials. Their GCD, stored in `ch`, will be restored at the end.

The variable `gamma` is the GCD of the leading coefficients of the primitive input polynomials. It is used to normalize the modular GCD images consistently before Chinese remainder reconstruction.

The function also collects information about leading coefficients with respect to the different variables:

```python
badprimes = ring.domain.one
for i in range(k):
    badprimes *= ring.domain.gcd(
        _swap(f, i).LC,
        _swap(g, i).LC
    )
```

Primes dividing this value are skipped. Reducing modulo such primes could cause relevant leading coefficients to vanish and change the degree information of the inputs.

### Choosing a prime and reducing the inputs

The function then enters a loop over suitable prime numbers:

```python
p = nextprime(p)
while badprimes % p == 0:
    p = nextprime(p)

fp = f.trunc_ground(p)
gp = g.trunc_ground(p)
```

The polynomials `fp` and `gp` are the images of the primitive inputs over the finite field $\mathbb{F}_p$.

The current implementation can dispatch either to the existing Brown modular algorithm or to Zippel's algorithm. In the Zippel branch, it calls:

```python
hp = sparse_gcd(fp, gp, p)
if hp is not None:
    hp = hp.monic()
```

At this level, `sparse_gcd` is responsible for computing the multivariate GCD modulo the chosen prime.

Zippel's algorithm uses many random evaluation points. For this reason, the current implementation starts its prime search from a large value, making it extremely unlikely that the finite field runs out of usable evaluation points.

This is an implementation choice rather than a fundamental requirement of the algorithm, and it may be refined later.

### Reconstructing the integer coefficients

A GCD image modulo one prime is not sufficient to recover the integer result. The function therefore computes images modulo several primes and combines them using the Chinese remainder theorem:

```python
hm = _chinese_remainder_reconstruction_multivariate(
    hp, hlastm, p, m
)
m *= p
```

Here, `hlastm` contains the result reconstructed modulo the product of the previously used primes, while `hp` is the new image modulo `p`.

The reconstruction continues until two consecutive candidates coincide:

```python
if not hm == hlastm:
    hlastm = hm
    continue
```

This equality is used as an early indication that the integer coefficients have stabilized. It is not, however, considered a sufficient correctness check.

### Verifying the reconstructed candidate

Once the reconstruction stabilizes, the candidate is made primitive and tested by exact division:

```python
h = hm.primitive()[1]
fquo, frem = f.div(h)
gquo, grem = g.div(h)

if not frem and not grem:
    ...
    return h, cff, cfg
```

The candidate is accepted only if it divides both original primitive polynomials.

This is an important feature of the modular algorithm: modular computations and random evaluations are used to construct a likely candidate efficiently, but exact division is used before returning the final result.

## `sparse_gcd`: the recursive core

The function `sparse_gcd` computes the GCD over one fixed finite field:

$$
A,B\in\mathbb{F}_p[x_1,\ldots,x_n].
$$

It implements the recursive structure of Zippel's algorithm.

### The univariate base case

The recursion terminates when the polynomial ring contains only one variable:

```python
if ring.ngens == 1:
    G = _gf_gcd(A, B, p)
    return G
```

In this case, no interpolation is needed. The GCD can be computed directly using the existing univariate finite-field algorithm.

Every recursive step reduces the number of variables by one, so the complete multivariate computation is ultimately reduced to several calls to this base case.

### Removing content with respect to the last variable

Before evaluating the last variable, the function separates the primitive parts and contents of the inputs:

```python
a, A = _primitive(A, p)
b, B = _primitive(B, p)

c = _gf_gcd(a, b, p)
c = c.set_ring(ring)

g = _gf_gcd(_LC(A), _LC(B), p)
```

Here, `A` and `B` are viewed as polynomials in the first variables whose coefficients depend on the last variable.

The function `_primitive` extracts the corresponding content and primitive part. The GCD of the contents is stored in `c` and will be multiplied back into the result.

The value `g` is the GCD of the leading coefficients of the primitive parts. It is used to normalize the evaluated GCD images.

### Constructing the first recursive image

The function chooses a random value `t` for the last variable:

```python
t = random.randint(1, p - 1)
gk = g._evaluate({0: t}) % p

if t in M or gk == 0:
    continue
```

The point is rejected if it has already been used or if it makes the normalization factor vanish.

The last variable is then evaluated at `t`:

```python
xk = ring.gens[-1]

A_ = A.evaluate(xk, t).trunc_ground(p)
B_ = B.evaluate(xk, t).trunc_ground(p)
```

The resulting polynomials belong to a ring with one fewer variable.

The function now calls itself recursively:

```python
G = sparse_gcd(A_, B_, p)
```

By the inductive hypothesis, this produces a GCD image in one fewer variable.

If the recursive GCD is equal to one, no interpolation is necessary:

```python
if G == G.ring.one:
    return c
```

Otherwise, this first image provides the assumed monomial structure of the result.

### Normalizing the first image

The recursive GCD is normalized with:

```python
G = (G.monic() * gk).trunc_ground(p)
```

The call to `G.monic()` removes the arbitrary scalar factor of the computed finite-field GCD. Multiplication by `gk` then gives the image a normalization compatible with the leading-coefficient information extracted before the evaluation.

The normalized image is passed to `skeleton_sorter`:

```python
G_s, h, monic, pseudomonic = skeleton_sorter(G)
```

The function groups the monomials according to their degree in the main variable and stores the initial coefficients in a form suitable for interpolation.

At this point:

- `G_s` describes the assumed sparse support;
- `h` stores the first Newton coefficients;
- `monic` records whether the skeletal GCD is monic in the main variable.

The detailed construction of the skeleton and the interpolation systems was discussed in the previous post.

## Why the monic assumption matters

It is useful to clarify what the monic case means here.

The input polynomials themselves do not necessarily have to be monic. The important condition is that the GCD being reconstructed is monic with respect to the chosen main variable.

After evaluating all variables except $x_1$, a univariate GCD algorithm returns a monic representative. In the monic case, these representatives can be normalized consistently and their coefficients can be used directly as interpolation data.

In the non-monic case, this is no longer automatic.

Suppose that the desired evaluated GCD is $H_r(x_1)$. The univariate computation may instead return

$$
\widehat{H}_r(x_1)=m_rH_r(x_1),
$$

where $m_r$ is a nonzero scalar.

The factor $m_r$ can be different for each evaluation point. Therefore, the coefficients collected from different univariate images do not all belong to the same scale.

If they were inserted directly into the Vandermonde systems, the right-hand sides would be inconsistent with one another, and the reconstructed coefficients would generally be incorrect.

The non-monic algorithm must therefore determine these scaling factors as additional unknowns before completing the sparse interpolation. That normalization procedure is implemented in the same pull request, but I will discuss it separately in a future post.

For the remainder of this post, I assume that the monic normalization is available and consistent.

### Computing further images

After obtaining the first skeletal GCD, `sparse_gcd` chooses new values for the last variable:

```python
t = random.randint(1, p - 1)
gk = g._evaluate({0: t}) % p

if t in M or gk == 0:
    continue
```

For every new point, the inputs are evaluated again:

```python
A_ = A.evaluate(xk, t).trunc_ground(p)
B_ = B.evaluate(xk, t).trunc_ground(p)
```

Instead of computing their full multivariate GCD recursively from scratch, the function now uses the assumed skeleton:

```python
G_ = zippel_interp(
    A_, B_, G_s, p, monic, pseudomonic
)
```

The purpose of `zippel_interp` is to reconstruct the new GCD image using sparse interpolation.

If the assumed support is incompatible with the new evaluation or if a usable interpolation system cannot be constructed, the function returns `None` and another evaluation point is chosen.

### Normalizing the new image

`zippel_interp` returns a dictionary containing the reconstructed coefficients rather than a `PolyElement`.

The new image is normalized using the coefficient of the leading monomial:

```python
G_k = pow(G_[lc_mon], -1, p)

for mon in G_:
    G_[mon] = ((G_[mon] * G_k) * gk) % p
```

This is the dictionary-level equivalent of making the image monic and then multiplying it by the normalization factor `gk`.

After this step, the first and subsequent GCD images are expressed using the same normalization and can be combined through Newton interpolation.

### Incrementally restoring the evaluated variable

Each coefficient of the skeletal GCD is treated as a polynomial in the last variable.

For a new evaluation point `t`, `G_` provides one new value for each of these coefficient polynomials.

The function updates their Newton representations independently:

```python
vk = incremental_newton_interp(
    M[:len(h[i])],
    h[i],
    t,
    coeff,
    p
)
```

The list `M` contains the evaluation points already used, while `h[i]` contains the Newton coefficients constructed so far for the `i`-th coefficient of the skeleton.

If the new coefficient `vk` is nonzero, the interpolation has changed:

```python
if vk != 0:
    repeat = True
    h[i].append(vk)
```

The algorithm must therefore compute another GCD image.

If `vk` is zero, the existing interpolation already predicts the new evaluated coefficient correctly:

```python
else:
    skippable.add(i)
```

That coefficient is added to `skippable`, so it does not need to be updated during later iterations.

This is an implementation-level optimization: different coefficients can reach their final degree at different times, and there is no need to continue interpolating coefficients that have already stabilized.

### Reconstructing the polynomial

When every new Newton coefficient is zero, `repeat` remains false. The function then converts each interpolation from Newton form to the usual polynomial representation:

```python
pol = from_newt_to_poly(
    M[:len(h[i])],
    h[i],
    p
)
```

The coefficients of these univariate polynomials are inserted into a dictionary representing the full multivariate candidate:

```python
for j, el in enumerate(pol):
    gcd[mon + (j,)] = el

gcd = ring.from_dict(gcd)
```

The newly restored exponent is appended to the monomial tuple, which reintroduces the last variable into the result.

The candidate is made primitive and converted to a ring over the finite field:

```python
_, gcd = _primitive(gcd, p)
R_mod = ring.clone(domain=FF(p))
```

Finally, it is checked by exact division:

```python
A_mod = A.set_ring(R_mod)
B_mod = B.set_ring(R_mod)
gcd_mod = gcd.set_ring(R_mod)

if A_mod.rem(gcd_mod) == 0 and B_mod.rem(gcd_mod) == 0:
    return (gcd * c).trunc_ground(p)
```

If the candidate divides both inputs, the content `c` is restored and the finite-field GCD is returned.

If the division test fails, the algorithm continues with new evaluation points rather than accepting an incorrect reconstruction.

## `zippel_interp`: reconstructing one GCD image

The role of `zippel_interp` is narrower than that of `sparse_gcd`.

It does not restore the last variable and it does not perform Chinese remainder reconstruction. It works at one fixed value of the outer evaluation variable and reconstructs a single GCD image from univariate GCD computations.

Suppose that its inputs belong to

$$
\mathbb{F}_p[x_1,\ldots,x_k].
$$

The function leaves $x_1$ symbolic and evaluates

$$
x_2,\ldots,x_k
$$

at successive powers of a randomly chosen tuple.

### Extracting basic information

The function starts by identifying the main variable and the degree of the first input in that variable:

```python
ngens = A.ring.ngens
x1 = A.ring.gens[0]

deg_a = A.degree(x1)
lc_A = A.coeff_wrt(x1, deg_a)
```

The leading coefficient `lc_A` is a polynomial in the remaining variables. Evaluation tuples that make it vanish must be rejected, because they could lower the degree of the evaluated input in $x_1$.

The function also determines the number of monomials in each coefficient group of the skeleton:

```python
lengths = list(len(el) for el in G.values())
nt = lengths[-1]
```

Because `skeleton_sorter` orders the groups by size, `nt` is the largest number of unknown coefficients that must be reconstructed in one Vandermonde system.

In the monic case, this determines the number of required univariate GCD images:

```python
if monic:
    z = nt
```

### Choosing a usable tuple

The function randomly chooses a tuple for the variables other than $x_1$:

```python
t = tuple(
    random.randint(-p//2, p//2)
    for _ in range(ngens - 1)
)
```

The tuple must satisfy two main conditions.

First, the leading coefficient of `A` must remain nonzero at all powers of the tuple used during interpolation.

Second, the monomials belonging to the same coefficient group must evaluate to distinct values. Otherwise, the associated Vandermonde matrix would be singular.

The code constructs the corresponding Vandermonde values:

```python
vand_bas = []

for mon in el:
    j = 1
    for i in mon[1:]:
        j = (j * pow(t[i[0]], i[1], p)) % p

    if j in vand_bas:
        is_bad_tuple = True
        break

    vand_bas.append(j)
```

If two monomials produce the same value, the tuple is discarded and a new one is selected.

This is one of the probabilistic parts of the algorithm: random points are used to find a convenient structured system, but unsuitable choices are detected and rejected.

### Preparing repeated evaluations

The inputs must be evaluated at several powers of the same tuple.

The current implementation reorganizes each polynomial according to the degree in $x_1$:

```python
A_flat = list({} for _ in range(deg_a + 1))
B_flat = list({} for _ in range(B.degree(x1) + 1))
```

For every power of $x_1$, the corresponding coefficient polynomial is converted into a dictionary whose keys are the monomial values at the base tuple and whose values are the original scalar coefficients.

Conceptually, a coefficient of the form

$$
c_1M_1+\cdots+c_sM_s
$$

is represented using the values

$$
M_1(t),\ldots,M_s(t).
$$

Evaluation at the $r$-th power of the tuple then becomes

$$
c_1M_1(t)^r+\cdots+c_sM_s(t)^r.
$$

This avoids recomputing every multivariate monomial from scratch for each interpolation point.

The precise low-level representation may change as the sparse polynomial infrastructure is refactored, but its algorithmic purpose is to make repeated evaluations efficient.

### Computing the univariate GCD images

For each required power of the tuple, the function builds two univariate polynomials:

```python
for i in range(z):
    A_ev = []
    B_ev = []
    ...
```

The resulting coefficient lists are passed to the finite-field univariate GCD routine:

```python
G_ev = gf_gcd(
    A_ev[::-1],
    B_ev[::-1],
    p,
    A.ring.domain
)

G_ev.reverse()
```

The coefficients of `G_ev` are then grouped according to their powers of $x_1$:

```python
for deg in G.keys():
    eval_points[deg].append(G_ev[deg])
```

For a fixed degree `deg`, `eval_points[deg]` contains the right-hand side of the Vandermonde system associated with the corresponding coefficient of the skeletal GCD.

The function also checks that the evaluated GCD does not contain unexpected powers of $x_1$:

```python
for j, el in enumerate(G_ev):
    if el != 0 and j not in G.keys():
        return None
```

If a new nonzero term appears outside the assumed support, the skeleton is incompatible with the current image and the interpolation attempt is abandoned.

### Solving the Vandermonde systems

In the monic case, no additional scaling factors have to be recovered. The collected coefficients can therefore be used directly.

For every coefficient group of the skeleton, the function constructs the Lagrange basis:

```python
v = lag_basis(all_vand_basis[i], p)
```

It then solves the transposed Vandermonde system:

```python
c = vandermonde_interp(
    v,
    eval_points[deg][:len(v)],
    p,
    trans=True
)
```

The resulting values are the unknown scalar coefficients of the corresponding sparse coefficient polynomial.

They are inserted into the output dictionary:

```python
for j, mon in enumerate(monoms):
    C[mon[0]] = (
        c[j]
        * pow(all_vand_basis[i][j], -1, p)
    ) % p
```

After all groups have been reconstructed, `zippel_interp` returns one complete GCD image for the current value of the variable handled by the outer `sparse_gcd` call.

## Probabilistic construction and deterministic verification

Several parts of the implementation rely on random choices:

- the primes used by the outer modular algorithm;
- the values assigned to the last variable;
- the tuples used for sparse interpolation.

A random choice can be unsuitable because it causes a leading coefficient to vanish, makes two monomials indistinguishable or produces an image incompatible with the assumed skeleton.

These events do not cause an incorrect result to be returned. They cause the current attempt to be discarded or restarted.

The final candidates are checked through exact polynomial division at two levels:

- `sparse_gcd` checks the candidate modulo $p$;
- `modgcd_multivariate` checks the reconstructed candidate over the integers.

The implementation therefore uses probabilistic methods to construct candidates efficiently, but deterministic divisibility tests before accepting them.

## Current status

The monic path now combines all the required parts of the algorithm:

- recursive reduction of the number of variables;
- extraction of a skeletal GCD;
- sparse interpolation through structured Vandermonde systems;
- incremental Newton interpolation;
- finite-field trial division;
- Chinese remainder reconstruction;
- final exact verification over the integers.

The implementation is included in pull request [#29824](https://github.com/sympy/sympy/pull/29824). The pull request is still marked as a draft because some comments, documentation and style details require cleanup, but the complete monic algorithm is implemented and working.

The next major difficulty is the non-monic case. There, the univariate GCD images cannot be assumed to share the same normalization, so their unknown scaling factors must be reconstructed together with the coefficients of the skeletal GCD. I will discuss that extension in a separate post.

## References

1. Suling Yang, *Computing the Greatest Common Divisor of Multivariate Polynomials over Finite Fields*.
2. Keith O. Geddes, Stephen R. Czapor and George Labahn, *Algorithms for Computer Algebra*.
3. [SymPy pull request #29824](https://github.com/sympy/sympy/pull/29824).
