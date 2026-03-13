# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision
floating-point arithmetic** to compute roots of ill-conditioned
polynomials. This tool primarily addresses **rounding and
finite-precision issues** present in standard 64-bit floating‑point
implementations by constructing and solving a companion matrix within
the `mpmath` environment.

The solver **cannot eliminate the intrinsic mathematical
ill-conditioning** of certain polynomials, but it can **significantly
reduce the numerical errors caused by rounding**, making the computed
roots far more reliable.

It features **modular code structure**, **residual checks**, and
**robust input validation** to improve usability and reliability.

------------------------------------------------------------------------

## 1. The Mathematical Background

Polynomial root-finding algorithms are sensitive to both **intrinsic
conditioning of the polynomial** and **finite precision arithmetic**.
The mapping from coefficients to roots can be extremely sensitive to
perturbations, a phenomenon illustrated by **Wilkinson's example**.

For a polynomial:

$$
P(x) = \sum_{i=0}^{n} a_i x^i = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0
$$

small changes in coefficients may produce large changes in the roots.
This sensitivity is an inherent mathematical property and **cannot be
fully eliminated by any algorithm**.

However, many practical errors arise from **floating‑point rounding
limitations** in standard numerical implementations. By using
**arbitrary precision arithmetic**, this solver reduces these rounding
effects and therefore improves the stability of the computed roots.

As a classic illustration, the Wilkinson Polynomial ($n=20$) shows that
even extremely small perturbations (e.g., modifying the $x^{19}$
coefficient by $2^{-23}$) can shift real roots into the complex plane.

------------------------------------------------------------------------

## 2. Methodology

Instead of iterative methods such as Newton-Raphson or
Jenkins-Traub---which may struggle with convergence on clustered
roots---this solver utilizes the **Companion Matrix Eigenvalue Method**.

### 1. Normalization

The polynomial is converted to monic form by dividing all coefficients
by the leading coefficient.

### 2. Matrix Construction

A companion matrix $C$ is built such that its characteristic polynomial
is exactly $P(x)$:

$$
C =
\begin{pmatrix}
0 & 0 & \dots & 0 & -a_0/a_n \\
1 & 0 & \dots & 0 & -a_1/a_n \\
0 & 1 & \dots & 0 & -a_2/a_n \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \dots & 1 & -a_{n-1}/a_n
\end{pmatrix}
$$

### 3. High-Precision Eigenvalue Solve

The eigenvalues $\lambda$ of $C$ are computed using `mpmath.eig`.\
These eigenvalues are exactly the **roots of the polynomial**.

By default:

    mp.dps = 100

This provides 100 decimal digits of precision, which significantly
reduces rounding errors compared to standard double precision.

### 4. Residual Verification

Each computed root $r$ is validated by evaluating:

$$
|P(r)|
$$

This residual provides a practical measure of the numerical accuracy of
the computed root.

------------------------------------------------------------------------

## 3. Features

- **Arbitrary Precision:** Default 100 decimal places (configurable).
- **Improved Numerical Stability:** Mitigates rounding errors common
  in double‑precision implementations.
- **Modular Structure:** Separate functions for input, matrix
  construction, eigenvalue computation, evaluation, printing, and
  plotting.
- **Robust Input Handling:** Re-prompts until valid space-separated
  numeric coefficients are entered.
- **Residual Verification:** Displays $|P(r)|$ for each root to assess
  accuracy.
- **Complex Visualization:** Roots are plotted in the complex plane
  with a unit circle reference.

The plotting logic intentionally remains simple to avoid distortions
when visualizing **highly asymmetric root distributions**, such as those
produced by Wilkinson-type polynomials.

------------------------------------------------------------------------

## 4. Usage

Run the script and enter the coefficients separated by spaces.

Example input:

    Coefficients (space separated): 1 0 -4

The program outputs:

- Polynomial degree
- Polynomial equation string
- Computed roots
- Residual values $|P(r)|$
- A complex-plane plot of the roots with a unit circle reference

------------------------------------------------------------------------

## 5. Test Cases

### 1. Wilkinson's Polynomial (Roots 1 through 20)

    1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500 11310276995381 -135585182899530 1307535010540395 -10142299865511450 63030812099294896 -311333643161390640 1206647803780373360 -3599979517947607200 8037811822645051776 -12870931245150988800 13803759753640704000 -8752948036761600000 2432902008176640000

### 2. Multiple Clustered Roots ((x-1)\^10)

    1 -10 45 -120 210 -252 210 -120 45 -10 1

### 3. Large Dynamic Range

    1 0 0 0 0 100000000000000000000

------------------------------------------------------------------------

## Key Idea

The solver does **not attempt to solve the fundamental conditioning
problem of polynomial root finding**. Instead, it focuses on **reducing
numerical rounding errors** using high‑precision arithmetic, which in
practice leads to far more reliable root approximations for difficult
polynomials.
