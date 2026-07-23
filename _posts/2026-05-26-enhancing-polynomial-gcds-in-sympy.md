---
layout: post
title: "Enhancing Polynomial GCDs in SymPy: Project Overview"
---

This is the first post of my GSoC blog for SymPy. This and the next few posts are retrospective posts, to keep track with a chronological structure of my GSoC work. This post is meant as a general overview of the project and of the main ideas behind it. In the next posts, I will go into more detail about the concrete pull requests and the implementation work.

This summer I am working on improving polynomial GCD computation in SymPy, with a particular focus on polynomials represented sparsely.

Polynomial GCD computation is a fundamental operation in computer algebra systems. Even when users are not calling `gcd` directly, it can appear internally in many other algorithms. Because of this, improving polynomial GCDs can have an impact beyond a single function: it can affect simplification, solvers, integration, matrix algorithms over polynomial domains, and other parts of SymPy.

The main motivation for this project is that the current choice of algorithms for polynomial GCDs in SymPy can be limited in some important cases, especially for multivariate sparse polynomials. Sparse polynomials are polynomials where the number of nonzero terms is much smaller than the number of possible terms of the same degree. In such cases, algorithms that treat the polynomial as if it were dense can do much more work than necessary.

The project is organized around four main goals.

## 1. Zippel's algorithm

The first major goal is to work on an implementation of Zippel's sparse modular GCD algorithm.

The general idea behind modular GCD algorithms is to avoid doing all computations directly over the integers. Instead, the input polynomials are mapped modulo prime numbers, where computations are usually cheaper. After computing enough modular images, the result over the integers can be reconstructed using the Chinese remainder theorem.

For multivariate polynomials, one can also evaluate variables at suitable points, compute GCDs of smaller polynomials, and then recover the full result by interpolation. Zippel's algorithm improves this process for sparse polynomials by using sparse interpolation: rather than reconstructing all possible terms of a dense polynomial, it tries to recover only the terms that actually appear.

This is especially relevant for the kind of inputs where the GCD has a sparse structure. In those cases, using the sparsity of the polynomial is not just a small optimization, but can be the difference between a practical and an impractical computation.

An important part of this work is handling both the monic and the non-monic case. The non-monic case introduces a normalization problem: the GCDs computed after evaluation are only determined up to multiplication by a nonzero constant, and these scaling factors have to be handled carefully during interpolation. This makes the algorithm more subtle than the monic case.

## 2. Sparse subresultant PRS

The second goal is related to the sparse subresultant polynomial remainder sequence algorithm.

SymPy already has strong polynomial functionality, and previous work had been done toward a sparse subresultant PRS implementation. Part of my project is to revive and refine that work, so that SymPy can have a sparse version of this algorithm available.

This is useful because Zippel's algorithm is not meant to replace every other GCD algorithm (it works only for polynomials with integer coefficients for example, and even in that case, it could not be the best algorithm for certain inputs). Different algorithms can perform better on different kinds of inputs. 
A sparse PRS algorithm gives SymPy another option for polynomial GCD computation, especially in cases where it is better suited than existing alternatives (The PRS algorithm works with coefficients in various rings other than Z).

## 3. Benchmarking

The third goal is benchmarking.

After implementing or improving algorithms, it is not enough to know that they are mathematically correct. We also need to understand when they are fast, when they are slow, and what kinds of inputs they are best suited for.

For this reason, an important part of the project is to create benchmarks for polynomial GCD algorithms. The relevant parameters include the number of variables, the degree of the polynomials, the domain of the coefficients, and the sparsity of the input.

The goal is to gather enough data to make informed decisions about algorithm selection, instead of relying only on theoretical expectations or isolated examples.

## 4. Heuristics and dispatching

The final goal is to use the information from the benchmarks to improve dispatching.

In practice, a user should not have to manually choose the best GCD algorithm for a pair of polynomials. SymPy should inspect the input and select a good strategy automatically.

This requires some preprocessing and analysis of the input polynomials. For example, the algorithm selection may depend on whether the polynomials are univariate or multivariate, dense or sparse, over the integers or over another domain, and on the degrees involved.

The long-term idea is to build a better infrastructure that can analyze the inputs, optionally preprocess them, and dispatch them to the most appropriate GCD algorithm.

## Why this matters

The reason I chose this project is that polynomial GCD computation is a central operation in a computer algebra system. Improving it can make SymPy faster and more reliable in a broad range of situations.

I also like that this project combines several aspects of computer algebra: algebraic theory, algorithm design, and practical integration into a large existing codebase.

In the next posts, I will write about the concrete work I have done so far. I will focus on the pull requests, the design choices, the difficulties I encountered, and the parts of the implementation that changed while discussing the work with the SymPy community.
