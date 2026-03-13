import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, mpf, matrix, eig
def solve_and_plot(dps=50): # Allow customizable precision
   mp.dps = dps # Set decimal precision for robustness against ill-conditioning

   print("--- Robust Companion Matrix Polynomial Solver ---")
   print("Enter coefficients (e.g., '1 0 -4' for x^2 - 4 = 0)")

   try:
       user_input = input("Coefficients: ").strip()
       if not user_input:
           return

       # Convert to list of mpf for high precision
       raw_coeffs = [mpf(x) for x in user_input.split()]
       # Trim leading zeros (highest degree non-zero)
       i = 0
       while i < len(raw_coeffs) - 1 and raw_coeffs[i] == 0:
           i += 1
       coeffs = raw_coeffs[i:]

       # Validation: Need at least two coefficients
       if len(coeffs) < 2:
           print("Error: Need at least two coefficients for a valid polynomial.")
           return

       # Build equation string with only non-zero terms
       n = len(coeffs) - 1 # degree
       terms = []
       for i in range(n + 1):
           c = coeffs[i]
           if c == 0:
               continue
           deg = n - i
           # Sign handling
           if not terms:
               sign = '' if c > 0 else '-'
           else:
               sign = ' + ' if c > 0 else ' - '
           abs_c = abs(c)
           # Coefficient string
           if deg == 0 or abs_c != 1:
               coeff_str = mp.nstr(abs_c, 4)
           else:
               coeff_str = ''
           # Variable part
           if deg > 1:
               var_str = f'x^{deg}'
           elif deg == 1:
               var_str = 'x'                                                                   
           else:
               var_str = ''
           term = f"{sign}{coeff_str}{var_str}"
           terms.append(term)
       equation = ''.join(terms) + ' = 0'

       # 1. Normalize to monic
       leading = coeffs[0]
       if leading == 0:
           print("Error: Leading coefficient cannot be zero.")
           return
       a = [c / leading for c in coeffs[1:]] # Coefficients for companion (without leading 1)
       n = len(a)

       # 2. Build the Companion Matrix
       if n == 1:
           # Linear: x + a0 = 0 => root = -a0
           roots = [-a[0]]
           roots_np = np.array(roots, dtype=complex) # For plotting
       else:
           # Create companion matrix                   1 -10 45 -120 210 -252 210 -120 45 -10 1
           C = matrix(n, n)
           for row in range(1, n):
               C[row, row-1] = mpf(1)
           for row in range(n):
               C[row, n-1] = -a[n-1 - row] # Reversed for standard companion

           # Compute eigenvalues (roots) with high precision
           eigenvalues = eig(C, left=False, right=False)

           # Convert to Python complex for plotting
           roots = [complex(ev.real, ev.imag) for ev in eigenvalues]
           roots_np = np.array(roots)
													   1 -10 45 -120 210 -252 210 -120 45 -10 1
       # Sort roots by real part, then imag (groups conjugates and orders reals)
       roots_np = np.sort_complex(roots_np)

       # 3. Display Results
       print(f"\nPolynomial Degree: {n}")
       print(f"Equation: {equation}")
       print("-" * 30)
       for i, r in enumerate(roots_np, 1):
           print(f"Root {i}: {r.real:10.4f} {'+' if r.imag >= 0 else '-'} {abs(r.imag):.4f}j")

       # 4. Plotting
       plt.figure(figsize=(8, 8))
       plt.axhline(0, color='black', lw=1)
       plt.axvline(0, color='black', lw=1)

       # Plot roots
       plt.scatter(roots_np.real, roots_np.imag, color='red', marker='.', s=100, label='Roots', zorder=5)

       # Unit Circle for reference
       t = np.linspace(0, 2*np.pi, 100)
       plt.plot(np.cos(t), np.sin(t), ls='--', color='gray', alpha=0.5, label='Unit Circle')
       plt.grid(True, linestyle=':', alpha=0.6)
       plt.axis('equal')
       plt.title(f"Roots in Complex Plane\n{equation}")
       plt.xlabel("Real")
       plt.ylabel("Imaginary")
       plt.legend()
       plt.show()

   except ValueError:
       print("Error: Please enter only numbers separated by spaces.")
   except Exception as e:
       print(f"An unexpected error occurred: {e}")
if __name__ == "__main__":
   solve_and_plot()
