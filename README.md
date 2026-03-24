# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision
floating-point arithmetic** to compute roots of ill-conditioned
polynomials.

This tool primarily addresses **rounding and finite-precision issues**
present in standard 64-bit floating-point implementations by constructing
and solving a companion matrix within the `mpmath` environment.

The solver **cannot eliminate intrinsic mathematical ill-conditioning**,
but it can **significantly reduce numerical rounding errors**, making the
computed roots more reliable.

---

## 1. Mathematical Background

Polynomial root-finding is sensitive to both:

- intrinsic conditioning
- finite precision arithmetic

Given:

$$
P(x) = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0
$$

Small coefficient perturbations can produce large root changes
(**Wilkinson's phenomenon**).

### Multiple Roots

Near a root of multiplicity $m$:

$$
|\Delta x| \sim |\Delta a|^{1/m}
$$

---

## 2. Methodology

### Normalization

Polynomial is converted to monic form.

### Exact Zero-Root Deflation

Trailing zero coefficients are removed **purely algebraically**.
Zero roots are reinserted exactly.

### Companion Matrix

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

Eigenvalues correspond to polynomial roots.

### High-Precision Eigenvalue Solve

- Uses `mpmath.eig`
- Default: `mp.dps = 100`

### Linear Case

Degree-1 polynomials are solved directly (no matrix construction).

### Display Normalization (cosmetic only)

Tiny imaginary parts (below a relative tolerance) are set to zero for readability.
This primarily removes numerical noise from theoretically real roots, but may also suppress very small imaginary components in near-real roots.

Example: High-Multiplicity Root Behaviour
For polynomials with high-multiplicity roots, such as:

$(x-1)^{10}$

all roots are theoretically real and identical.
However, in numerical computation:
the multiple root typically splits into a cluster of nearby roots
some roots may acquire small imaginary parts
others may remain real after display normalization
This behavior is not a numerical bug, but a consequence of:
intrinsic ill-conditioning of multiple roots
sensitivity described by:

$$
|\Delta x| \sim |\Delta a|^{1/m}
$$

Effect of Precision
Increasing precision:
reduces random rounding noise
shrinks the cluster
but does not eliminate the root splitting
This provides a direct visualization of the instability of multiple roots.

    --- Robust Companion Matrix Polynomial Solver ---
Coefficients (space separated): 1 -10 45 -120 210 -252 210 -120 45 -10 1.
NOTICE: Moderately ill-conditioned matrix (cond ≈ 1.065e+6)  
NOTICE: Clustered roots detected → high sensitivity likely.

Polynomial Degree: 10
Equation:
    x^10 - 10x^9 + 45x^8 - 120x^7 + 210x^6 - 252x^5 + 210x^4 - 120x^3 + 45x^2 - 10x + 1 = 0

---

    Root #  Real       Imaginary   |z|     Rel.Residual

    1        1.0            +0.0j   1.0     6.13212e-102
    2        1.0            +0.0j   1.0     2.23937e-102
    3        1.0            +0.0j   1.0     1.31851e-102
    4        1.0            +0.0j   1.0     4.97406e-102
    5        1.0     -1.0029e-10j   1.0      6.9069e-102
    6        1.0    +1.39469e-10j   1.0     6.31195e-102
    7        1.0    -1.48179e-10j   1.0     6.14457e-102
    8        1.0    +1.48179e-10j   1.0     6.12143e-102
    9        1.0    -1.39469e-10j   1.0     6.64146e-102
    10       1.0     +1.0029e-10j   1.0     6.01245e-102

**Important:**  
No artificial numerical stabilization is applied:

- no balancing
- no scaling
- no polishing
- no deflation (except exact zero roots)

---

## 3. Diagnostics

The solver reports issues without fixing them.

### Coefficient Scaling

- NOTICE: ratio > 1e6
- WARNING: ratio > 1e12

### Matrix Conditioning (approximate indicator)

- NOTICE: cond > 1e6
- WARNING: cond > 1e12

### Cluster Detection

Roots closer than:

$$
10^{-8} \cdot \max(|r|, 1)
$$

### Residual Warnings

Expected scale:

$$
10^{-\frac{\text{dps}}{2} + 5}
$$

---

## 4. Residual Definition

$$
\frac{|P(r)|}{\sum |a_i| \cdot |r|^{n-i}}
$$

Small residuals indicate numerical consistency,
but not necessarily accurate roots for ill-conditioned problems.

---

## 5. Features

- Arbitrary precision (configurable)
- Exact zero-root handling
- Companion matrix eigenvalue method
- Transparent diagnostics
- Relative residual verification
- Visualization (complex plane + real domain)

---

## 6. Usage

Example input:

   1 0 -4 --> x^2 - 4 = 0

   1 1 4  --> x^2 + x + 4 = 0

   1 1 0  --> x^2 + x = 0

Outputs:

- diagnostics
- polynomial equation
- roots with residuals
- plots

Controls:

Toggle y-sxale of real domain

- "l" linear
- "y" symlog

---

## 7. Interpretation

- Small residual → numerically consistent root
- Warnings → structural/numerical difficulty
- Clusters → high sensitivity
- Conditioning → approximate matrix indicator (not eigenvalue sensitivity)

---

## 8. Test Cases

- Wilkinson polynomial
- Multiple roots: (x-1)^10, (x-1)^20
- Complex multiplicity: (x^2 + 1)^5
- Zero-root cases
- Extreme scaling cases

### 1. Classic Wilkinson Polynomial (distinct roots 1–20)

    1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500 11310276995381 -135585182899530 1307535010540395 -10142299865511450 63030812099294896 -311333643161390640 1206647803780373360 -3599979517947607200 8037811822645051776 -12870931245150988800 13803759753640704000 -8752948036761600000 2432902008176640000

Expected diagnostics: extreme coefficient scaling + ill-conditioned matrix.

### 2. Wilkinson-style with multiplicity 4 (roots 1–5 each multiplicity 4)

    1 -60 1690 -29700 365071 -3334800 23477380 -130374600 579693631 -2082967740 6077950570 -14418383700 27741179521 -43026101880 53234263960 -51699564000 38463221776 -21114635520 8041766400 -1893888000 207360000

Expected diagnostics: scaling, ill-conditioned matrix, clustered roots (exactly as in the example run).

### 3. Extreme high-multiplicity real root

$(x-1)^{10}$

    1 -10 45 -120 210 -252 210 -120 45 -10 1

Remark: Compare default `dps=100`with `dps=200`!

$(x-1)^{20}$

    1 -20 190 -1140 4845 -15504 38760 -77520 125970 -167960 184756 -167960 125970 -77520 38760 -15504 4845 -1140 190 -20 1

Remark: Compare default `dps=100`with `dps=400`!

### 4. High-multiplicity complex roots: $(x^2 + 1)^5$

    1 0 5 0 10 0 10 0 5 0 1

Expected: clustered roots notice (five identical roots at $i$ and at $-i$).

### 5. Explicit root at zero

    1 -5 0

Expected:  
`NOTICE: Polynomial has 1 root(s) at x = 0`  
The zero root appears with **Rel.Residual = 0** and there is **no**
large-residual warning and **no** singular-matrix message.

### 6. Variant without leading coefficient (negative leading)

    -20 190 -1140 4845 -15504 38760 -77520 125970 -167960 184756 -167960 125970 -77520 38760 -15504 4845 -1140 190 -20 1

### 7. Large dynamic range / extreme scaling

    1 0 0 0 0 100000000000000000000
---
    1 0 -1000101 0 101000100 0 -100000000
---
    524288 0 -2621440 0 +5570560 0 -6553600 0 +4845120 0 -2325120 0 +715400 0 -133760 0 +14400 0 -800 0 1

---

## Key Idea

The solver does **not solve the conditioning problem**.

It instead provides:

- high precision arithmetic
- raw companion matrix method
- full numerical transparency

---

## Further Reading

<https://medium.com/@ratwolf/the-floating-point-catastrophe-9e795d46cfb1>
