---
layout: post
title: "The Interpolation Machinery Behind Zippel's Algorithm"
---

Before working on the complete implementation of Zippel's algorithm, I implemented some of the interpolation tools required by it.

These tools were submitted in separate pull requests. This makes them easier to test and review independently, but their purpose becomes clearer when they are considered together: they implement the two interpolation stages that make the inductive step of Zippel's algorithm possible.

The first one is **sparse interpolation**, which reconstructs a polynomial whose monomial structure is already known.

The second one is **dense Newton interpolation**, which restores a variable that had previously been evaluated.

In this post, I will describe how these two forms of interpolation fit into the recursive structure of Zippel's algorithm and how the helper functions I implemented support this process.

## The inductive structure of Zippel's algorithm

Zippel's algorithm computes a multivariate polynomial GCD recursively.

Suppose that we want to compute

$$
G = \gcd(A,B),
$$

where

$$
A,B \in \mathbb{F}_p[x_1,\ldots,x_n].
$$

Here, $\mathbb{F}_p$ is a finite field. In the complete modular algorithm, the original integer polynomials are first reduced modulo a prime $p$, and the results obtained for different primes are later combined using the Chinese remainder theorem.

The recursive idea is to temporarily evaluate the last variable $x_n$ at some value $a_0$. This produces two polynomials in one fewer variable:

$$
A_{a_0}
=
A(x_1,\ldots,x_{n-1},a_0),
$$

$$
B_{a_0}
=
B(x_1,\ldots,x_{n-1},a_0).
$$

By the inductive hypothesis, we assume that we can compute

$$
G_{a_0}
=
\gcd(A_{a_0},B_{a_0})
$$

$$\text{in} \quad \mathbb{F}_p[x_1,\ldots,x_{n-1}]$$.

(The base case of the induction is a univariate gcd, and therefore it is easily taken care of.)

Under suitable choices of the evaluation point, the monomials appearing in $G_{a_0}$ give us the expected structure of the full GCD. This first image is therefore used as a **skeletal GCD**.

The skeletal GCD tells us which monomials should appear, but it does not yet tell us how their coefficients depend on $x_n$.

To recover that information, the algorithm performs two different reconstruction steps:

1. for new values $a_i$ of $x_n$, sparse interpolation reconstructs new images

   $$
   G_{a_i}\in\mathbb{F}_p[x_1,\ldots,x_{n-1}]
   $$

   using the known skeleton;

2. once enough images $G_{a_i}$ have been found, dense interpolation reconstructs their coefficients as polynomials in $x_n$.

Together, these steps turn a GCD in $n-1$ variables into a GCD in $n$ variables.

Remark: $G_{a_i}$ is in general different from the real gcd $G$ evaluated in $a_i$: they differ by multiplication of a scalar. In order to perform dense interpolation succesfully, we need them to always differ by multiplication of the same scalar. Therefore a normalization has to be applied, but for the sake of this post we can omit this problem, and suppose that $G(x_1,\ldots,a_i) = G_{a_i}$. I will treat in more detail this normalization prolem in one of the next posts. 

## Organizing the skeletal GCD

For sparse interpolation, the GCD is considered as a polynomial in the first variable $x_1$, with coefficients in the remaining variables.

We write it as

$$
G
=
C_kx_1^k+\cdots+C_1x_1+C_0,
$$

where each coefficient $C_j$ belongs to

$$
\mathbb{F}_p[x_2,\ldots,x_{n-1}].
$$

The first variable $x_1$ is kept symbolic throughout the sparse interpolation step. It is not evaluated.

Each $C_j$ is itself a sparse polynomial. Since the skeletal GCD tells us which monomials occur, we can write

$$
C_j
=
M_{1,j}W_{1,j}
+\cdots+
M_{t_j,j}W_{t_j,j},
$$

where

$$
W_{s,j}
=
x_2^{a_{2,s,j}}\cdots x_{n-1}^{a_{n-1,s,j}}
$$

are known monomials, while

$$
M_{1,j},\ldots,M_{t_j,j}
$$

are the unknown scalar coefficients that have to be reconstructed.

The helper function `skeleton_sorter` performs this reorganization.

It groups the monomials of the skeletal GCD according to their degree in $x_1$. It also orders these groups by the number of unknown coefficients they contain and creates a compact representation of their nonzero exponents.

The reorganization is mostly preparatory: the important mathematical point is that, after it has been performed, sparse interpolation can treat each coefficient $C_j$ independently.

## From multivariate inputs to univariate GCDs

Let us now fix a new value $a_i$ for $x_n$. We want to reconstruct

$$
G_{a_i}
=
G(x_1,\ldots,x_{n-1},a_i)
$$

using the monomial structure obtained from the skeletal GCD.

To find the unknown coefficients $M_{s,j}$, the variables

$$
x_2,\ldots,x_{n-1}
$$

are evaluated at suitable points, while $x_1$ is left unchanged.

Choose a tuple

$$
b=(b_2,\ldots,b_{n-1}).
$$

The inputs are evaluated at successive powers of this tuple:

$$
b^r
=
(b_2^r,\ldots,b_{n-1}^r).
$$

For every $r$, we therefore compute the univariate GCD

$$
H_r(x_1)
=
\gcd\left(
 A(x_1,b_2^r,\ldots,b_{n-1}^r,a_i),
 B(x_1,b_2^r,\ldots,b_{n-1}^r,a_i)
\right).
$$

Both evaluated inputs are univariate polynomials in $x_1$, so this GCD can be computed with a univariate Euclidean algorithm.

Write the result as

$$
H_r(x_1)
=
h_{r,k}x_1^k+\cdots+h_{r,0}.
$$

For a fixed power $x_1^j$, the value $h_{r,j}$ is the evaluated value of the corresponding skeletal coefficient $C_j$:

$$
h_{r,j}
=
C_j(b_2^r,\ldots,b_{n-1}^r).
$$

Substituting the sparse representation of $C_j$, we obtain

$$
h_{r,j}
=
M_{1,j}W_{1,j}(b^r)
+\cdots+
M_{t_j,j}W_{t_j,j}(b^r).
$$

This is a linear equation in the unknown coefficients

$$
M_{1,j},\ldots,M_{t_j,j}.
$$

By calculating enough univariate GCD images, we obtain enough equations to determine them.

The coefficients of the skeletal GCD are therefore the unknowns of the systems, while the coefficients of the evaluated univariate GCDs form their right-hand sides.

## Why the systems are Vandermonde systems

The evaluation points are chosen as powers of the same tuple precisely because this gives the systems a special structure.

For a monomial

$$
W_{s,j}
=
x_2^{a_2}\cdots x_{n-1}^{a_{n-1}},
$$

we have

$$
W_{s,j}(b^r)
=
(b_2^r)^{a_2}\cdots(b_{n-1}^r)^{a_{n-1}}.
$$

This can be rewritten as

$$
W_{s,j}(b^r)
=
\left(
 b_2^{a_2}\cdots b_{n-1}^{a_{n-1}}
\right)^r
=
W_{s,j}(b)^r.
$$

If we define

$$
\beta_{s,j}=W_{s,j}(b),
$$

then the equations for $C_j$ become

$$
\begin{aligned}
M_{1,j}\beta_{1,j}
+\cdots+
M_{t_j,j}\beta_{t_j,j}
&=h_{1,j},\\
M_{1,j}\beta_{1,j}^2
+\cdots+
M_{t_j,j}\beta_{t_j,j}^2
&=h_{2,j},\\
&\vdots\\
M_{1,j}\beta_{1,j}^{t_j}
+\cdots+
M_{t_j,j}\beta_{t_j,j}^{t_j}
&=h_{t_j,j}.
\end{aligned}
$$

The associated matrix is

$$
\begin{pmatrix}
\beta_{1,j} & \beta_{2,j} & \cdots & \beta_{t_j,j}\\
\beta_{1,j}^2 & \beta_{2,j}^2 & \cdots & \beta_{t_j,j}^2\\
\vdots & \vdots & \ddots & \vdots\\
\beta_{1,j}^{t_j} & \beta_{2,j}^{t_j} & \cdots &
\beta_{t_j,j}^{t_j}
\end{pmatrix}.
$$

Up to the choice of whether the powers begin from zero or one, this is the transposed of a Vandermonde matrix.

This structure is valuable because a Vandermonde system can be solved in quadratic time, rather than by applying a generic linear-system solver.

The mathematical treatment of this interpolation method is described in *Computing the Greatest Common Divisor of Multivariate Polynomials over Finite Fields* by Suling Yang.

## The Vandermonde helper functions

I implemented two helper functions for this stage:

- `lag_basis`;
- `vandermonde_interp`.

The function `lag_basis` constructs the Lagrange basis associated with the values

$$
\beta_{1,j},\ldots,\beta_{t_j,j}.
$$

The basis depends only on these evaluation values. It's for this reason that I decided to write 2 distinct functions: supposing we have a Vandermonde system $Ax = b$, which we want to solve for different values of $b$, we would need to compute the basis only once.
Once it has been computed, it can be used to solve the corresponding interpolation problem.

The function `vandermonde_interp` then combines this basis with the known values

$$
h_{1,j},\ldots,h_{t_j,j}
$$

to recover

$$
M_{1,j},\ldots,M_{t_j,j}.
$$

The same procedure is performed independently for every power of $x_1$ appearing in the skeletal GCD.

After all these systems have been solved, the algorithm has reconstructed one complete image

$$
G_{a_i}
$$

for the chosen value $x_n=a_i$.

In the monic case, the coefficients of the univariate GCD images can be used directly as the right-hand sides of the Vandermonde systems.

In the general non-monic case, the images may differ by unknown scaling factors. Those factors must first be determined through an additional linear system. Once they are known, the same Vandermonde interpolation machinery can be used for the remaining coefficients. I will discuss this normalization problem separately in a later post.

## Restoring the last variable

Sparse interpolation reconstructs a copy of the GCD for one fixed value of $x_n$:

$$
G_{a_i}
=
G(x_1,\ldots,x_{n-1},a_i).
$$

The process is then repeated for different evaluation points

$$
a_0,a_1,a_2,\ldots.
$$

Because all these images share the same skeletal structure, every monomial of the skeleton receives a scalar coefficient at each evaluation point.

Suppose, for example, that one monomial $W$ appears in every image with coefficients

$$
u_0,u_1,u_2,\ldots,
$$

where

$$
u_i=Q(a_i)
$$

for some unknown polynomial $Q(x_n)$.

The coefficient of $W$ in the full GCD is not merely a scalar: it is the polynomial $Q(x_n)$. Recovering $Q$ from the values $Q(a_i)$ is now an ordinary dense interpolation problem in one variable.

This is the second interpolation stage of the inductive step.

Sparse interpolation has determined which monomials are present and has reconstructed their scalar coefficients at a fixed value of $x_n$. Dense interpolation now restores the dependence of those coefficients on $x_n$.

## Incremental Newton interpolation

The dense reconstruction is performed using Newton interpolation.

After the first $k$ evaluation points, an interpolation polynomial can be written in Newton form as

$$
P_k(x_n)
=
v_0
+
v_1(x_n-a_0)
+
v_2(x_n-a_0)(x_n-a_1)
+\cdots+
v_k\prod_{r=0}^{k-1}(x_n-a_r).
$$

When a new GCD image is computed at $a_{k+1}$, it is not necessary to repeat the full interpolation from the beginning. Only the next Newton coefficient has to be calculated.

The function `incremental_newton_interp` performs this update. It receives:

- the previous evaluation points;
- the Newton coefficients already computed;
- the new point $a_{k+1}$;
- the newly reconstructed coefficient value.

It then calculates the next Newton coefficient.

This is well suited to Zippel's algorithm because the GCD images are generated one at a time. Every new image extends the current interpolation rather than replacing it.

When a new Newton coefficient is zero, the current interpolation already predicts the newly obtained value. This provides a termination criterion: once the interpolations have stabilized, the algorithm can attempt to reconstruct the complete GCD and verify it by trial division.

The function `from_newt_to_poly` performs the final conversion from Newton form to the usual expanded polynomial representation.

The Newton interpolation used here is described in *Algorithms for Computer Algebra* by Keith Geddes, Stephen Czapor and George Labahn.

## Completing the inductive step

The two interpolation methods solve different parts of the same recursive problem.

For a fixed value $a_i$ of $x_n$, sparse interpolation performs the following steps:

1. keep $x_1$ as the main variable;
2. evaluate $x_2,\ldots,x_{n-1}$ at powers of a suitable tuple;
3. compute univariate GCDs in $x_1$;
4. extract, for each power $x_1^j$, the corresponding coefficients of those univariate GCDs;
5. use these coefficients as the right-hand sides of Vandermonde systems;
6. solve for the unknown coefficients of the skeletal GCD;
7. obtain the complete image

   $$
   G(x_1,\ldots,x_{n-1},a_i).
   $$

This gives one copy of the GCD associated with one evaluation point of $x_n$.

The process is repeated for several values $a_i$. Dense Newton interpolation then uses the reconstructed coefficients from these different copies to recover each coefficient as a polynomial in $x_n$.

The final result is a polynomial in

$$
x_1,\ldots,x_n.
$$

This completes the inductive passage from $n-1$ variables to $n$ variables.


## References

1. Suling Yang, *Computing the Greatest Common Divisor of Multivariate Polynomials over Finite Fields*.
2. Keith O. Geddes, Stephen R. Czapor and George Labahn, *Algorithms for Computer Algebra*.
