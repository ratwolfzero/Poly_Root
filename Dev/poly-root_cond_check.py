import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, matrix, eig


def parse_coefficients(text):
    """Parse coefficient string into a valid coefficient list."""
    parts = text.strip().split()
    if not parts:
        raise ValueError("Please enter at least two numbers.")
    try:
        coeffs = [mpf(x) for x in parts]
    except ValueError:
        raise ValueError(
            "Invalid input. Please enter numbers separated by spaces.")
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
    """Compute roots using companion matrix.
    Safely handles singular matrices (e.g. root at 0) and returns friendly error only when needed."""
    C, linear_root = build_companion(coeffs)
    if linear_root is not None:
        return [mpc(linear_root[0])], None

    # === SAFE condition number (mp.cond can raise on exact singularity) ===
    cond = None
    try:
        cond = mp.cond(C)
    except ZeroDivisionError as e:
        if "numerically singular" in str(e).lower():
            cond = mp.inf          # treat exact singularity as "infinite condition"
        else:
            raise

    # === Now compute eigenvalues (eig is more tolerant than cond) ===
    try:
        eigenvalues = eig(C, left=False, right=False)
        roots_mp = [mpc(ev) for ev in eigenvalues]
        roots_mp.sort(key=lambda z: (mp.re(z), mp.im(z)))
        return roots_mp, cond

    except ZeroDivisionError as e:
        if "numerically singular" in str(e).lower():
            raise RuntimeError(
                f"Companion matrix is NUMERICALLY SINGULAR (cond = ∞)\n\n"
                "This happens when:\n"
                "   • Constant term is exactly (or extremely close to) zero → multiple root at 0\n"
                "   • Polynomial is severely ill-conditioned\n"
                "   • Coefficients lead to a defective (Jordan) matrix\n\n"
                "SOLUTIONS:\n"
                "   1. Increase precision dramatically: solve_and_plot(dps=500) or dps=1000\n"
                "   2. If constant term = 0, factor out x manually and solve the reduced polynomial\n"
                "   3. Check coefficients — tiny changes can move roots a lot"
            ) from e
        raise


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
    """Convert mpmath mpc to NumPy complex."""
    return np.array(
        [complex(float(mp.re(z)), float(mp.im(z))) for z in roots_mp],
        dtype=complex
    )


def plot_roots(roots_mp, equation):
    """Plot roots in complex plane."""
    roots_np = mpc_to_numpy_complex(roots_mp)
    plt.figure(figsize=(8, 8))
    plt.axhline(0, color="black", lw=1)
    plt.axvline(0, color="black", lw=1)
    plt.scatter(roots_np.real, roots_np.imag, color="red", marker=".", s=100,
                label="Roots", zorder=5)
    t = np.linspace(0, 2 * np.pi, 200)
    plt.plot(np.cos(t), np.sin(t), ls="--", color="gray", alpha=0.5,
             label="Unit Circle")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.axis("equal")
    plt.title(f"Roots in Complex Plane\n{equation}", fontsize=10)
    plt.xlabel("Real")
    plt.ylabel("Imaginary")
    plt.legend()
    plt.show()


def solve_and_plot(dps=16):

    mp.dps = dps

    print("\n--- Robust Companion Matrix Polynomial Solver ---")
    print("   (safe condition number + full singularity protection)")

    coeffs = get_coefficients()
    equation = polynomial_string(coeffs)

    roots_mp, cond = compute_roots(coeffs)

    print(f"\nPolynomial Degree: {len(coeffs) - 1}")
    print(f"Equation: {equation}")

    # === Display condition number + warnings ===
    if cond is not None:
        if mp.isinf(cond):
            print(f"\nCompanion matrix condition number: ∞")
            print(
                "   → Matrix is singular (e.g. root(s) exactly at 0 or multiple roots).")
            print("Matrix inversion failed → extremely ill-conditioned (not necessarily singular)")
        elif cond > mpf('1e12'):
            print(f"\nCompanion matrix condition number: {mp.nstr(cond, 6)}")
            print("   → Polynomial is ill-conditioned. Roots may be very sensitive.")
            print(
                "   → Consider increasing dps (e.g. solve_and_plot(dps=500)) for better accuracy.")

    print_roots(coeffs, roots_mp)
    plot_roots(roots_mp, equation)


if __name__ == "__main__":
    solve_and_plot()
