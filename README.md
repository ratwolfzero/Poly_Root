# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision floating-point arithmetic** to compute roots of ill-conditioned polynomials. This tool bypasses the precision limitations of standard 64-bit float implementations by constructing a companion matrix within the `mpmath` environment.

## 1. The Mathematical Problem

Standard root-finding algorithms are notoriously sensitive to small perturbations in coefficients, a phenomenon known as **Wilkinson's Catastrophe**. For a polynomial:

$$P(x) = \sum_{i=0}^{n} a_i x^i = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0$$

The mapping from coefficients $\vec{a}$ to roots $\vec{r}$ can be extremely ill-conditioned. In the case of the Wilkinson Polynomial ($n=20$), a change in the $x^{19}$ coefficient by as little as $2^{-23}$ can shift real roots into the complex plane.

## 2. The Methodology

Instead of iterative Newton-Raphson or Jenkins-Traub methods which may struggle with convergence on clustered roots, this solver utilizes the **Companion Matrix Eigenvalue Method**:

1. **Normalization:** The polynomial is converted to its monic form by dividing all coefficients by $a_n$.
2. **Matrix Construction:** A companion matrix $C$ is built such that its characteristic polynomial is exactly $P(x)$:
   $$C = \begin{pmatrix} 0 & 0 & \dots & 0 & -a_0/a_n \\ 1 & 0 & \dots & 0 & -a_1/a_n \\ 0 & 1 & \dots & 0 & -a_2/a_n \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \dots & 1 & -a_{n-1}/a_n \end{pmatrix}$$
3. **High-Precision Eigensolve:** We calculate the eigenvalues $\lambda$ of $C$ using `mpmath.eig`, which are precisely the roots of $P(x)$. By setting `mp.dps = 50`, we maintain enough guard bits to prevent precision loss during the QR/Schur decomposition.

## 3. Features

* **Arbitrary Precision:** Default 50 decimal places (configurable).
* **Robust Logic:** Trims leading zeros and handles linear/edge cases.
* **Complex Visualization:** Automatic plotting of roots against the unit circle in the complex plane.

## Test

### 1. Wilkinson's Polynomial (Roots 1 through 20)

1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500 11310276995381 -135585182899530 1307535010540395 -10142299865511450 63030812099294896 -311333643161390640 1206647803780373360 -3599979517947607200 8037811822645051776 -12870931245150988800 13803759753640704000 -8752948036761600000 2432902008176640000

### 2. Multiple Clustered Roots ((x-1)^10 = 0)

1 -10 45 -120 210 -252 210 -120 45 -10 1

### 3. Large Dynamic Range

1 0 0 0 0 100000000000000000000
