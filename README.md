# High-Precision Companion Matrix Root Finder

A robust polynomial solver leveraging **arbitrary-precision floating-point arithmetic** to compute roots of ill-conditioned polynomials. This tool bypasses the precision limitations of standard 64-bit float implementations by constructing a companion matrix within the `mpmath` environment.

It features **modular code structure**, **residual checks**, and **robust input validation** to improve usability and reliability.

---

## 1. The Mathematical Problem

Standard root-finding algorithms are notoriously sensitive to small perturbations in coefficients, a phenomenon known as **Wilkinson's Catastrophe**. For a polynomial:

$$
P(x) = \sum_{i=0}^{n} a_i x^i = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0
$$

The mapping from coefficients $\vec{a}$ to roots $\vec{r}$ can be extremely ill-conditioned. In the case of the Wilkinson Polynomial ($n=20$), a change in the $x^{19}$ coefficient by as little as $2^{-23}$ can shift real roots into the complex plane.

---

## 2. Methodology

Instead of iterative Newton-Raphson or Jenkins-Traub methods which may struggle with convergence on clustered roots, this solver utilizes the **Companion Matrix Eigenvalue Method**:

1. **Normalization:** The polynomial is converted to its monic form by dividing all coefficients by the leading coefficient.  
2. **Matrix Construction:** A companion matrix $C$ is built such that its characteristic polynomial is exactly $P(x)$:

$$
C = \begin{pmatrix}
0 & 0 & \dots & 0 & -a_0/a_n \\
1 & 0 & \dots & 0 & -a_1/a_n \\
0 & 1 & \dots & 0 & -a_2/a_n \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \dots & 1 & -a_{n-1}/a_n
\end{pmatrix}
$$

1. **High-Precision Eigensolve:** Compute the eigenvalues $\lambda$ of $C$ using `mpmath.eig`, which are the roots of $P(x)$. By default, `mp.dps = 100` ensures enough precision for ill-conditioned polynomials.  
2. **Residual Check:** Each computed root $r$ is validated by evaluating $|P(r)|$, providing a measure of numerical accuracy.

---

## 3. Features

* **Arbitrary Precision:** Default 100 decimal places (configurable).  
* **Modular & Structured:** Functions separate input, companion matrix construction, eigenvalue computation, evaluation, printing, and plotting.  
* **Robust Input Handling:** Loops until valid space-separated numeric coefficients are entered.  
* **Residual Verification:** Prints $|P(r)|$ for each root to confirm numerical correctness.  
* **Complex Visualization:** Roots are plotted in the complex plane with a unit circle reference. Original plotting logic preserved to avoid distortion for highly asymmetric root distributions (e.g., Wilkinson polynomial).  

---

## 4. Usage

Run the script and enter the coefficients separated by spaces. Example:
Coefficients (space separated): 1 0 -4

Output includes:

* Polynomial degree  
* Polynomial equation string  
* List of roots with residuals $|P(r)|$  
* Complex-plane plot of roots with unit circle  

---

## 5. Test Cases

### 1. Wilkinson's Polynomial (Roots 1 through 20)

1 -210 20615 -1256850 53327946 -1672280820 40171771630 -756111184500
11310276995381 -135585182899530 1307535010540395 -10142299865511450
63030812099294896 -311333643161390640 1206647803780373360
-3599979517947607200 8037811822645051776 -12870931245150988800
13803759753640704000 -8752948036761600000 2432902008176640000

### 2. Multiple Clustered Roots ((x-1)^10)

1 -10 45 -120 210 -252 210 -120 45 -10 1

### 3. Large Dynamic Range

1 0 0 0 0 100000000000000000000
