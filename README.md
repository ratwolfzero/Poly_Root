# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision
floating-point arithmetic** to compute roots of ill-conditioned
polynomials.
This tool primarily addresses **rounding and
finite-precision issues** present in standard 64-bit floating-point
implementations by constructing and solving a companion matrix within
the `mpmath` environment.
The solver **cannot eliminate the intrinsic mathematical
ill-conditioning** of certain polynomials, but it can **significantly
reduce the numerical errors caused by rounding**, making the computed
roots far more reliable.
It features **modular code structure**, **residual checks**, and
**robust input validation** to improve usability and reliability

---

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
that even extremely small perturbations (for example modifying the
$x^{19}$ coefficient by $2^{-23}$) can shift real roots into the complex
plane.

### Multiple Roots and Sensitivity

Polynomials containing **multiple roots** are particularly
ill-conditioned. Near a root of multiplicity $m$, small perturbations in
coefficients can produce root shifts approximately proportional to

$$
|\Delta x| \sim |\Delta a|^{1/m}
$$

This means even extremely small coefficient perturbations may split a
multiple root into a cluster of nearby roots.
High‑precision arithmetic helps ensure that such behavior reflects the
**true mathematical sensitivity of the polynomial**, rather than
artificial numerical noise

---

## 2. Methodology

Instead of relying on iterative root‑polishing techniques such as
Newton--Raphson or specialized polynomial solvers like Jenkins--Traub,
this solver computes all roots simultaneously using the **Companion
Matrix Eigenvalue Method**.

### Normalization

The polynomial is converted to monic form by dividing all coefficients
by the leading coefficient.

### Matrix Construction

A companion matrix $C$ is constructed such that its **characteristic
polynomial corresponds to the normalized polynomial**:

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

The eigenvalues of this matrix correspond to the **roots of the
polynomial**.

### High‑Precision Eigenvalue Solve

The eigenvalues of the companion matrix are computed using `mpmath.eig`,
which operates with arbitrary precision.
By default
    mp.dps = 100
This provides **100 decimal digits of working precision**, which
significantly reduces rounding errors compared to standard double
precision.
Increasing `mp.dps` increases the precision used during both the
eigenvalue computation and polynomial evaluation.
**Important (non-artificial approach):** The solver uses the **raw,
unmodified companion matrix** with **no artificial improvements**
whatsoever — no balancing, no scaling, no root polishing, no deflation,
no Wilkinson's shift, and no cluster handling. All numerical
difficulties are reported transparently via diagnostics (see section 3).

### Residual Verification

Each computed root $r$ is validated by computing the **relative residual**

$$
\frac{|P(r)|}{\sum |a_i| \cdot |r|^{n-i}}
$$

This measures **how well the computed value satisfies the polynomial
equation relative to its magnitude**. Residuals are displayed for every
root

---

## 3. Diagnostics (the "non-artificial but diagnostic" approach)

The solver **never hides or fixes** ill-conditioning. Instead it prints
clear diagnostics exactly when they occur:

- **Coefficient scaling** — ratio max/min > 10⁶ (NOTICE) or > 10¹² (WARNING)
- **Companion matrix conditioning** — 1-norm condition number > 10⁶ (NOTICE), > 10¹² or singular (WARNING)
- **Clustered roots** — groups closer than 10⁻⁸ × |root| (NOTICE)
- **Large residuals** — absolute or relative residuals exceeding the working-precision expectation (WARNING)
This philosophy guarantees you see the **true behaviour** of the raw
companion-matrix method on ill-conditioned polynomials.

---

## 4. Features

- **Arbitrary Precision:** Default 100 decimal places (configurable)
- **Improved Numerical Stability:** Mitigates rounding errors common
    in double‑precision implementations
- **Modular Structure:** Separate functions for input, matrix
    construction, eigenvalue computation, evaluation, printing, and
    plotting
- **Robust Input Handling:** Re‑prompts until valid space‑separated
    numeric coefficients are entered
- **Transparent Diagnostics:** Prints every warning/notice exactly as
    the raw method encounters it
- **Relative Residual Verification:** Displays relative residual for each root
- **Combined Visualization:** Roots plotted in the complex plane
    **plus** the polynomial curve on the real line (with automatic
    symlog scaling for huge dynamic ranges)

---

## 5. Usage

Run the script and enter the coefficients separated by spaces.
Example input
    Coefficients (space separated): 1 0 -4
The program outputs

- All diagnostics (scaling, conditioning, clusters, residuals)
- Polynomial degree
- Polynomial equation string
- Computed roots with relative residuals
- A combined plot (complex plane + real-domain curve) --> use the "l" and "y" key to toggle between linear and symlog scale

---

## 6. Interpreting Residuals

The displayed value is the **relative residual**. With

    mp.dps = 100

typical values lie far below the warning threshold of roughly

$$
10^{-45}
$$

A very small relative residual indicates that the computed root satisfies
the polynomial equation to near machine precision. However, this **does
not guarantee high root accuracy** in cases where the polynomial is
intrinsically ill-conditioned, such as for multiple roots or highly
sensitive coefficient-to-root mappings.

Diagnostics (such as conditioning warnings or singular companion
matrices) describe the numerical difficulty of the **companion matrix
eigenvalue problem**, not directly the accuracy of the computed roots.

Therefore, it is entirely possible to observe warnings while still
obtaining highly accurate roots, provided the relative residuals remain
very small. Conversely, small residuals should always be interpreted in
the context of the reported diagnostics, especially for ill-conditioned
problems.

---

## 7. Test Cases & Stress Tests

### 1. Classic Wilkinson Polynomial (distinct roots 1–20)

    1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500 11310276995381 -135585182899530 1307535010540395 -10142299865511450 63030812099294896 -311333643161390640 1206647803780373360 -3599979517947607200 8037811822645051776 -12870931245150988800 13803759753640704000 -8752948036761600000 2432902008176640000

    **Expected diagnostics:** extreme coefficient scaling + ill-conditioned matrix.

### 2. Wilkinson-style with multiplicity 4 (roots 1–5 each multiplicity 4)

    1 -60 1690 -29700 365071 -3334800 23477380 -130374600 579693631 -2082967740 6077950570 -14418383700 27741179521 -43026101880 53234263960 -51699564000 38463221776 -21114635520 8041766400 -1893888000 207360000

    **Expected diagnostics:** scaling, ill-conditioned matrix, clustered roots (exactly as in the example run).

### 3. Extreme high-multiplicity real root

$(x-1)^{10}$

    1 -10 45 -120 210 -252 210 -120 45 -10 1

Remark: Compare default `dps=100`with `dps=200`!

$(x-1)^{20}$

    1 -20 190 -1140 4845 -15504 38760 -77520 125970 -167960 184756 -167960 125970 -77520 38760 -15504 4845 -1140 190 -20 1

Remark: Compare default `dps=100`with `dps=400`!

### 4. High-multiplicity complex roots: $(x^2 + 1)^5$

    1 0 5 0 10 0 10 0 5 0 1

    **Expected:** clustered roots notice (five identical roots at $i$ and at $-i$).

### 5. Explicit root at zero (singular companion matrix)

    1 -5 0

    **Expected:** "NOTICE: Polynomial has root at x = 0 → companion matrix is singular."

### 6. Variant without leading coefficient (negative leading)

    -20 190 -1140 4845 -15504 38760 -77520 125970 -167960 184756 -167960 125970 -77520 38760 -15504 4845 -1140 190 -20 1

### 7. Large dynamic range / extreme scaling

    1 0 0 0 0 100000000000000000000

    1 0 -1000101 0 101000100 0 -100000000

    524288 0 -2621440 0 +5570560 0 -6553600 0 +4845120 0 -2325120 0 +715400 0 -133760 0 +14400 0 -800 0 1

    ------------------------------------------------------------------------

## Key Idea

The solver does **not attempt to solve the fundamental conditioning
problem of polynomial root finding**.
Instead, it focuses on **reducing numerical rounding errors** using
high‑precision arithmetic **while deliberately applying zero artificial
stabilizations**. The diagnostics let you see the raw, unfiltered
behaviour of the companion-matrix method.

---

## Further Reading

Medium article:
<https://medium.com/@ratwolf/the-floating-point-catastrophe-9e795d46cfb1>
