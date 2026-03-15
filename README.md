# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision
floating-point arithmetic** to compute roots of ill-conditioned
polynomials. This tool primarily addresses **rounding and
finite-precision issues** present in standard 64-bit floating-point
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

Polynomial root-finding algorithms are sensitive to both the **intrinsic
conditioning of the polynomial** and the **finite precision of
arithmetic**. The mapping from coefficients to roots can be extremely
sensitive to perturbations, a phenomenon illustrated by **Wilkinson's
example**.

For a polynomial

$$
P(x) = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0
$$

small changes in coefficients may produce large changes in the roots.
This sensitivity is an inherent mathematical property and **cannot be
fully eliminated by any algorithm**.

However, many practical errors arise from **floating-point rounding
limitations** in standard numerical implementations. By using
**arbitrary precision arithmetic**, this solver reduces these rounding
effects and therefore improves the stability of the computed roots.

As a classic illustration, the Wilkinson Polynomial ($n = 20$) shows
that even extremely small perturbations (e.g., modifying the $x^{19}$
coefficient by $2^{-23}$) can shift real roots into the complex plane.

### Multiple Roots and Sensitivity

Polynomials containing **multiple roots** are particularly
ill-conditioned. Near a root of multiplicity $m$, small perturbations in
coefficients can produce root shifts approximately proportional to

$$
|\Delta x| \sim |\Delta a|^{1/m}
$$

This means even extremely small coefficient perturbations may split a
multiple root into a cluster of nearby roots.

High-precision arithmetic helps ensure that such behavior reflects the
**true mathematical sensitivity of the polynomial**, rather than
artificial numerical noise.

------------------------------------------------------------------------

## 2. Methodology

Instead of relying on iterative root-polishing techniques such as
Newton-Raphson or specialized polynomial solvers like Jenkins--Traub,
this solver computes all roots simultaneously using the **Companion
Matrix Eigenvalue Method**.

### 1. Normalization

The polynomial is converted to monic form by dividing all coefficients
by the leading coefficient.

### 2. Matrix Construction

A companion matrix $C$ is constructed such that its **characteristic
polynomial equals the normalized polynomial**:

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

The eigenvalues $\lambda$ of $C$ are computed using `mpmath.eig`.

These eigenvalues are exactly the **roots of the polynomial**.

By default

    mp.dps = 100

This provides **100 decimal digits of working precision**, significantly
reducing rounding errors compared to standard double precision.

Increasing `mp.dps` increases the precision used during both the
eigenvalue computation and polynomial evaluation.

### 4. Residual Verification

Each computed root $r$ is validated by evaluating

$$
|P(r)|
$$

This residual measures **how well the computed value satisfies the
polynomial equation**.

Residual magnitudes often scale with the working precision (e.g. roughly
$10^{-dps}$), but the exact value may vary depending on the structure of
the polynomial and numerical conditioning.

A small residual indicates that $P(r) \approx 0$, but **does not always
guarantee that the root itself is highly accurate**, particularly for
ill-conditioned polynomials or multiple roots.

------------------------------------------------------------------------

## 3. Features

- **Arbitrary Precision:** Default 100 decimal places (configurable)
- **Improved Numerical Stability:** Mitigates rounding errors common
    in double-precision implementations
- **Modular Structure:** Separate functions for input, matrix
    construction, eigenvalue computation, evaluation, printing, and
    plotting
- **Robust Input Handling:** Re-prompts until valid space-separated
    numeric coefficients are entered
- **Residual Verification:** Displays $|P(r)|$ for each root
- **Complex Visualization:** Roots are plotted in the complex plane
    with a unit circle reference

The plotting logic intentionally remains simple to avoid distortions
when visualizing **highly asymmetric root distributions**, such as those
produced by Wilkinson-type polynomials.

------------------------------------------------------------------------

## 4. Usage

Run the script and enter the coefficients separated by spaces.

Example input

    Coefficients (space separated): 1 0 -4

The program outputs

- Polynomial degree
- Polynomial equation string
- Computed roots
- Residual values $|P(r)|$
- A complex-plane plot of the roots with a unit circle reference

------------------------------------------------------------------------

## 5. Test Cases

### 1. Wilkinson's Polynomial (Roots 1 through 20)

    1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500 11310276995381 -135585182899530 1307535010540395 -10142299865511450 63030812099294896 -311333643161390640 1206647803780373360 -3599979517947607200 8037811822645051776 -12870931245150988800 13803759753640704000 -8752948036761600000 2432902008176640000

### 2. Multiple Clustered Roots $((x-1)^10)$

    1 -10 45 -120 210 -252 210 -120 45 -10 1

### 3. Large Dynamic Range

    1 0 0 0 0 100000000000000000000

------------------------------------------------------------------------

## 6. Interpreting Residuals and Precision

Residual values provide useful diagnostic information but should be
interpreted carefully.

In arbitrary-precision arithmetic with `mp.dps = D`, rounding errors
typically occur at roughly

$$
10^{-D}
$$

Therefore residual magnitudes often appear near this scale.

However several factors influence the exact value:

- Polynomial structure
- Root multiplicity
- Conditioning of the polynomial
- Cancellation during polynomial evaluation

For well-conditioned polynomials, small residuals usually indicate that
the computed roots are accurate.

For ill-conditioned polynomials---especially those with **multiple
roots**---tiny residuals may still coexist with noticeable root errors.
This occurs because the polynomial becomes extremely flat near such
roots.

In other words:

- **Residuals measure equation satisfaction**
- **They do not directly measure root accuracy**

Consequently residuals should be interpreted as **diagnostic indicators
rather than absolute guarantees of correctness**.

------------------------------------------------------------------------

## Key Idea

The solver does **not attempt to solve the fundamental conditioning
problem of polynomial root finding**.

Instead, it focuses on **reducing numerical rounding errors** using
high-precision arithmetic, which in practice leads to far more reliable
root approximations for difficult polynomials.

------------------------------------------------------------------------

## Further Reading

Medium article:

<https://medium.com/@ratwolf/the-floating-point-catastrophe-9e795d46cfb1>
