import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp, mpf, mpc, matrix, eig


def get_coefficients():
    """Prompt user until valid coefficient input is provided."""
    while True:
        user_input = input("Coefficients (space separated): ").strip()

        if not user_input:
            print("Please enter at least two numbers.")
            continue

        try:
            coeffs = [mpf(x) for x in user_input.split()]
        except ValueError:
            print("Invalid input. Please enter numbers separated by spaces.")
            continue

        # Trim leading zeros
        i = 0
        while i < len(coeffs) - 1 and coeffs[i] == 0:
            i += 1
        coeffs = coeffs[i:]

        if len(coeffs) < 2:
            print("Need at least two coefficients.")
            continue

        return coeffs


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
    """Compute roots using companion matrix."""
    C, linear_root = build_companion(coeffs)

    if linear_root is not None:
        roots_mp = [mpc(linear_root[0])]
    else:
        eigenvalues = eig(C, left=False, right=False)
        roots_mp = [mpc(ev) for ev in eigenvalues]

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
            f"   |p(r)| ≈ {mp.nstr(residual, 6)}"
        )


def mpc_to_numpy_complex(roots_mp):
    """Convert a list of mpmath mpc numbers to a NumPy complex array."""
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

    plt.scatter(
        roots_np.real,
        roots_np.imag,
        color="red",
        marker=".",
        s=100,
        label="Roots",
        zorder=5,
    )

    t = np.linspace(0, 2 * np.pi, 200)

    plt.plot(
        np.cos(t),
        np.sin(t),
        ls="--",
        color="gray",
        alpha=0.5,
        label="Unit Circle",
    )

    plt.grid(True, linestyle=":", alpha=0.6)

    plt.axis("equal")

    plt.title(f"Roots in Complex Plane\n{equation}", fontsize=10)
    plt.xlabel("Real")
    plt.ylabel("Imaginary")
    plt.legend()

    plt.show()


def solve_and_plot(dps=100):
    mp.dps = dps

    print("\n--- Robust Companion Matrix Polynomial Solver ---")

    coeffs = get_coefficients()

    equation = polynomial_string(coeffs)

    roots_mp = compute_roots(coeffs)

    print(f"\nPolynomial Degree: {len(coeffs) - 1}")
    print(f"Equation: {equation}")

    print_roots(coeffs, roots_mp)

    plot_roots(roots_mp, equation)


if __name__ == "__main__":
    solve_and_plot()
