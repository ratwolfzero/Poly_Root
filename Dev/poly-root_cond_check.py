from mpmath import mp, mpf, mpc, matrix, eig
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Ensures TkAgg backend works well on macOS


def parse_coefficients(text):
    """Parse coefficient string into a valid coefficient list."""
    parts = text.strip().split()
    if not parts:
        raise ValueError("Please enter at least two numbers.")
    try:
        coeffs = [mpf(x) for x in parts]
    except ValueError:
        raise ValueError("Invalid input. Enter numbers separated by spaces.")
    # Trim leading zeros
    i = 0
    while i < len(coeffs) - 1 and coeffs[i] == 0:
        i += 1
    coeffs = coeffs[i:]
    if len(coeffs) < 2:
        raise ValueError("Need at least two coefficients.")
    return coeffs


def get_coefficients():
    """Prompt user until valid coefficient input is provided."""
    while True:
        try:
            text = input("Coefficients (space separated): ")
            return parse_coefficients(text)
        except ValueError as e:
            print(e)


def polynomial_string(coeffs):
    """Create a readable polynomial string."""
    def fmt(num):
        if num == int(num):
            return str(int(num))
        return mp.nstr(num, 4)
    n = len(coeffs) - 1
    terms = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        deg = n - i
        if not terms:
            sign = "" if c > 0 else "-"
        else:
            sign = " + " if c > 0 else " - "
        abs_c = abs(c)
        if deg == 0:
            coeff_str = fmt(abs_c)
        elif abs_c == 1:
            coeff_str = ""
        else:
            coeff_str = fmt(abs_c)
        var_str = f"x^{deg}" if deg > 1 else ("x" if deg == 1 else "")
        terms.append(f"{sign}{coeff_str}{var_str}")
    return "".join(terms) + " = 0"


def build_companion(coeffs):
    """Construct companion matrix."""
    leading = coeffs[0]
    if leading == 0:
        print("WARNING: Leading coefficient is zero - companion matrix construction failed.")
        raise ZeroDivisionError("Leading coefficient is zero.")
    a = [c / leading for c in coeffs[1:]]
    n = len(a)
    if n == 1:
        return None, [-a[0]]
    C = matrix(n)
    for row in range(1, n):
        C[row, row - 1] = mpf(1)
    for row in range(n):
        C[row, n - 1] = -a[n - 1 - row]
    return C, None


def compute_roots(coeffs):
    """Compute roots using companion matrix."""
    C, linear_root = build_companion(coeffs)
    if linear_root is not None:
        roots_mp = [mpc(linear_root[0])]
    else:
        try:
            eigenvalues = eig(C, left=False, right=False)
            roots_mp = [mpc(ev) for ev in eigenvalues]
        except Exception as e:
            print(f"WARNING: Companion matrix eigenvalue computation failed: {type(e).__name__} - {e}")
            print("Possible cases: invalid matrix dimensions, convergence failure, defective matrix, or numerical instability.")
            raise
    # Critical condition check (companion matrix accuracy)
    residuals = [abs(poly_eval(coeffs, r)) for r in roots_mp]
    if residuals:
        max_res = max(residuals)
        if max_res > mpf(1):
            print("WARNING: Companion matrix eigenvalue computation produced inaccurate roots.")
            print(f"   Max |p(r)| ≈ {mp.nstr(max_res, 6)} (expected ~0)")
            print("   Possible cases: low dps vs coefficient magnitude/degree (numerical instability),")
            print("   ill-conditioned polynomial, defective matrix, or rounding errors in input.")
    # Sort roots deterministically
    roots_mp.sort(key=lambda z: (mp.re(z), mp.im(z)))
    return roots_mp


def poly_eval(coeffs, x):
    """Evaluate polynomial using Horner's method."""
    p = mpf(0)
    for c in coeffs:
        p = p * x + c
    return p


def print_roots(coeffs, roots_mp):
    """Display roots and high-precision residuals."""
    print("-" * 40)
    for i, r in enumerate(roots_mp, 1):
        residual = abs(poly_eval(coeffs, r))
        r_real = float(mp.re(r))
        r_imag = float(mp.im(r))
        print(
            f"Root {i:2d}: {r_real:10.6f} "
            f"{'+' if r_imag >= 0 else '-'} {abs(r_imag):.6f}j"
            f" |p(r)| ≈ {mp.nstr(residual, 6)}"
        )

									 
def mpc_to_numpy_complex(roots_mp):
    """Convert a list of mpmath mpc numbers to a NumPy complex array."""
    return np.array(
        [complex(float(mp.re(z)), float(mp.im(z))) for z in roots_mp],
        dtype=complex
    )                                                                                    

									 
def plot_roots(roots_mp, equation, coeffs):
    """Complex plane plot with a button to open polynomial curve."""
    roots_np = mpc_to_numpy_complex(roots_mp)
    plt.figure(figsize=(8, 8))
    plt.axhline(0, color="black", lw=1)
    plt.axvline(0, color="black", lw=1)                                
    plt.scatter(roots_np.real, roots_np.imag, color="red",
                marker=".", s=100, label="Roots", zorder=5)

    t = np.linspace(0, 2 * np.pi, 200)
    plt.plot(np.cos(t), np.sin(t), ls="--", color="gray",
             alpha=0.5, label="Unit Circle")

    plt.grid(True, linestyle=":", alpha=0.6)
    plt.axis("equal")
    plt.title(f"Roots in Complex Plane\n{equation}", fontsize=10, pad=20)
    plt.xlabel("Real")
    plt.ylabel("Imaginary")
    plt.legend(loc='best')

    # Toggle button (Top-Right placement)
    button_ax = plt.axes([0.75, 0.02, 0.15, 0.03])
    button = Button(button_ax, 'Poly Curve',
                    color='lightgoldenrodyellow', hovercolor='gold')

    def open_poly(event):
        plot_polynomial_curve(coeffs, roots_mp, equation)

    button.on_clicked(open_poly)                                    
    plt.show()


def plot_polynomial_curve(coeffs, roots_mp, equation):
    """Smarter, automatically scaled polynomial curve."""
    fig_name = "PolyCurveWindow"

    # Avoid opening duplicates                                      
    if plt.fignum_exists(fig_name):
        fig = plt.figure(fig_name)
        try:
            fig.canvas.manager.window.lift()
        except:
            pass
        return

    coeffs_float = np.array([float(c) for c in coeffs], dtype=float)
    real_parts = [float(mp.re(r)) for r in roots_mp]

    # --- Smarter X-axis scaling ---
    spread = max(real_parts) - min(real_parts)                      
    x_center = sum(real_parts) / len(real_parts)
    x_pad = 1.5 * max(spread, 1.0)  # at least 1.0 unit padding
    x_min = x_center - x_pad
    x_max = x_center + x_pad
    x_vals = np.linspace(x_min, x_max, 2000)
    y_vals = np.polyval(coeffs_float, x_vals)

    # --- Smarter Y-axis scaling (95th percentile) ---
    y_abs_sorted = np.sort(np.abs(y_vals))
    y_limit = max(y_abs_sorted[int(0.95*len(y_abs_sorted))], 1e-9)

    # --- Plot ---
    plt.figure(num=fig_name, figsize=(10, 6))                       
    plt.plot(x_vals, y_vals, color="blue", lw=2, label="p(x)")
    plt.axhline(0, color="black", lw=1)
    plt.axvline(0, color="black", lw=0.5, ls="--")                     

    # Mark real roots
    tol = 1e-7
    real_roots_to_mark = [r for i, r in enumerate(
        real_parts) if abs(float(mp.im(roots_mp[i]))) < tol]
    if real_roots_to_mark:
        plt.plot(real_roots_to_mark, [0]*len(real_roots_to_mark),
                 "ro", markersize=8, label="Real Roots")

    plt.grid(True, linestyle=":", alpha=0.7)
    plt.title(f"Polynomial Curve p(x)\n{equation}", fontsize=11)
    plt.xlim(x_min, x_max)
    plt.ylim(-y_limit, y_limit)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show(block=False)


def solve_and_plot(dps=20):

    mp.dps = dps
    print("\n--- Robust Companion Matrix Polynomial Solver ---")
    coeffs = get_coefficients()                                     
    equation = polynomial_string(coeffs)
    roots_mp = compute_roots(coeffs)
    print(f"\nPolynomial Degree: {len(coeffs) - 1}")
    print(f"Equation: {equation}")
    print_roots(coeffs, roots_mp)

    # Complex-plane plot with built-in button
    plot_roots(roots_mp, equation, coeffs)


if __name__ == "__main__":
    solve_and_plot()
