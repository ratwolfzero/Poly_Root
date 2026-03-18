from mpmath import mp, mpf, mpc, matrix, eig
from matplotlib.widgets import Button
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')


# ----------------------------- Input Parsing ----------------------------- #
def parse_coefficients(text):
    parts = text.strip().split()
    if len(parts) < 2:
        raise ValueError("Please enter at least two numbers.")
    try:
        coeffs = [mpf(x) for x in parts]
    except ValueError:
        raise ValueError("Invalid input. Enter numbers separated by spaces.")

    while len(coeffs) > 1 and coeffs[0] == 0:
        coeffs.pop(0)

    if len(coeffs) < 2:
        raise ValueError(
            "Invalid polynomial: leading zeros reduce this to a constant (no roots).")
    return coeffs


def get_coefficients():
    while True:
        try:
            return parse_coefficients(input("Coefficients (space separated): "))
        except ValueError as e:
            print(e)


# ----------------------------- Diagnostics ----------------------------- #
def coefficient_diagnostics(coeffs):
    mags = [abs(c) for c in coeffs if c != 0]
    if not mags:
        return

    ratio = max(mags) / min(mags)

    if ratio > 1e12:
        print(
            f"WARNING: Extreme coefficient scaling (ratio ≈ {mp.nstr(ratio, 4)})")
    elif ratio > 1e6:
        print(
            f"NOTICE: Large coefficient scaling (ratio ≈ {mp.nstr(ratio, 4)})")


def companion_diagnostics(C):
    try:
        cond_number = mp.norm(C, 1) * mp.norm(C**-1, 1)
    except Exception:
        cond_number = mp.inf

    if cond_number == mp.inf:
        print("WARNING: Companion matrix is singular or near-singular.")
    elif cond_number > 1e12:
        print(
            f"WARNING: Ill-conditioned matrix (cond ≈ {mp.nstr(cond_number, 4)})")
    elif cond_number > 1e6:
        print(
            f"NOTICE: Moderately ill-conditioned matrix (cond ≈ {mp.nstr(cond_number, 4)})")

    return cond_number


def detect_clusters(roots, tol=1e-4):
    clusters = []
    used = [False] * len(roots)

    for i, r1 in enumerate(roots):
        if used[i]:
            continue
        cluster = [r1]
        used[i] = True

        for j, r2 in enumerate(roots):
            if not used[j] and abs(r1 - r2) < tol:
                cluster.append(r2)
                used[j] = True

        if len(cluster) > 1:
            clusters.append(cluster)

    return clusters


# ----------------------------- Polynomial Helpers ----------------------------- #
def polynomial_string(coeffs):
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
        sign = "" if not terms and c > 0 else (" + " if c > 0 else " - ")
        abs_c = abs(c)
        coeff_str = "" if deg > 0 and abs_c == 1 else fmt(abs_c)
        var_str = f"x^{deg}" if deg > 1 else ("x" if deg == 1 else "")
        terms.append(f"{sign}{coeff_str}{var_str}")
    return "".join(terms) + " = 0"


def poly_eval(coeffs, x):
    p = mpf(0)
    for c in coeffs:
        p = p * x + c
    return p


def relative_residual(coeffs, r):
    numerator = abs(poly_eval(coeffs, r))
    denom = mpf(0)
    n = len(coeffs) - 1
    for i, c in enumerate(coeffs):
        denom += abs(c) * abs(r) ** (n - i)
    return numerator / denom if denom != 0 else numerator


# ----------------------------- Companion Matrix ----------------------------- #
def build_companion(coeffs):
    leading = coeffs[0]
    if leading == 0:
        raise ZeroDivisionError("Leading coefficient cannot be zero.")

    if coeffs[-1] == 0:
        print("NOTICE: Polynomial has root at x = 0 → companion matrix is singular.")

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


# ----------------------------- Root Computation ----------------------------- #
def compute_roots(coeffs):
    coefficient_diagnostics(coeffs)

    C, linear_root = build_companion(coeffs)

    if linear_root is not None:
        roots_mp = [mpc(linear_root[0])]
    else:
        companion_diagnostics(C)

        try:
            eigenvalues = eig(C, left=False, right=False)
            roots_mp = [mpc(ev) for ev in eigenvalues]
        except Exception as e:
            print("ERROR: Eigenvalue computation failed.")
            print("   Likely due to extreme ill-conditioning.")
            print(f"   {type(e).__name__}: {e}")
            return []

    # Residual checks
    abs_residuals = [abs(poly_eval(coeffs, r)) for r in roots_mp]
    rel_residuals = [relative_residual(coeffs, r) for r in roots_mp]

    max_abs = max(abs_residuals) if abs_residuals else 0
    max_rel = max(rel_residuals) if rel_residuals else 0

    if max_abs > mpf(1):
        print(
            f"WARNING: Large absolute residuals (max ≈ {mp.nstr(max_abs, 6)})")

    if max_rel > mpf(1e-6):
        print(
            f"WARNING: Large relative residuals (max ≈ {mp.nstr(max_rel, 4)})")

    # Cluster detection (replaces old multiplicity heuristic)
    clusters = detect_clusters(roots_mp)
    if clusters:
        print("NOTICE: Clustered roots detected → high sensitivity likely.")

    roots_mp.sort(key=lambda z: (mp.re(z), mp.im(z)))
    return roots_mp


# ----------------------------- Display ----------------------------- #
def print_roots(coeffs, roots_mp):
    print("-" * 40)
    for i, r in enumerate(roots_mp, 1):
        residual = abs(poly_eval(coeffs, r))
        r_real, r_imag = float(mp.re(r)), float(mp.im(r))
        print(
            f"Root {i:2d}: {r_real:10.6f} "
            f"{'+' if r_imag >= 0 else '-'} {abs(r_imag):.6f}j"
            f" |p(r)| ≈ {mp.nstr(residual, 6)}"
        )


def mpc_to_numpy_complex(roots_mp):
    return np.array([complex(float(mp.re(z)), float(mp.im(z))) for z in roots_mp])


# ----------------------------- Plotting ----------------------------- #
def plot_roots(roots_mp, equation, coeffs):
    if not roots_mp:
        return

    roots_np = mpc_to_numpy_complex(roots_mp)

    plt.figure(figsize=(8, 8))
    plt.axhline(0, lw=1)
    plt.axvline(0, lw=1)
    plt.scatter(roots_np.real, roots_np.imag, color="red",
                s=30, label="Roots", zorder=5)

    t = np.linspace(0, 2 * np.pi, 200)
    plt.plot(np.cos(t), np.sin(t), ls="--", alpha=0.5, label="Unit Circle")

    plt.grid(True, linestyle=":", alpha=0.6)
    plt.axis("equal")
    plt.title(f"Roots in Complex Plane\n{equation}", fontsize=10, pad=20)
    plt.xlabel("Real")
    plt.ylabel("Imaginary")
    plt.legend(loc='best')

    button_ax = plt.axes([0.75, 0.02, 0.15, 0.03])
    button = Button(button_ax, 'Poly Curve')
    button.on_clicked(lambda event: plot_polynomial_curve(
        coeffs, roots_mp, equation))

    plt.show()


def plot_polynomial_curve(coeffs, roots_mp, equation):
    if not roots_mp:
        return

    coeffs_float = np.array([float(c) for c in coeffs])
    real_parts = [float(mp.re(r)) for r in roots_mp]

    spread = max(real_parts) - min(real_parts)
    x_center = sum(real_parts) / len(real_parts)
    x_pad = 1.5 * max(spread, 1.0)

    x_vals = np.linspace(x_center - x_pad, x_center + x_pad, 2000)
    y_vals = np.polyval(coeffs_float, x_vals)

    y_limit = max(np.sort(np.abs(y_vals))[int(0.95 * len(y_vals))], 1e-9)

    # Reuse window
    fig = plt.figure("Polynomial Curve", figsize=(10, 6))
    fig.clf()

    # Plot the polynomial curve
    plt.plot(x_vals, y_vals, lw=2, label="p(x)")

    # Axes
    plt.axhline(0, lw=1)
    plt.axvline(0, lw=1)

    # Plot real roots as discrete points
    tol = 1e-7
    real_roots = [
        float(mp.re(r)) for r in roots_mp if abs(float(mp.im(r))) < tol
    ]

    if real_roots:
        plt.scatter(
            real_roots,
            [0] * len(real_roots),
            color="red",
            s=60,
            label="Real Roots",
            zorder=5
        )

    # Styling
    plt.grid(True, linestyle=":", alpha=0.7)
    plt.title(f"Polynomial Curve p(x)\n{equation}")
    plt.ylim(-y_limit, y_limit)
    plt.xlabel("x")
    plt.ylabel("p(x)")
    plt.legend(loc='best')
    plt.tight_layout()

    plt.show(block=False)


# ----------------------------- Main ----------------------------- #
def solve_and_plot(dps=100):
    mp.dps = dps
    print("\n--- Robust Companion Matrix Polynomial Solver ---")

    coeffs = get_coefficients()
    equation = polynomial_string(coeffs)

    roots_mp = compute_roots(coeffs)

    print(f"\nPolynomial Degree: {len(coeffs) - 1}")
    print(f"Equation: {equation}")

    if roots_mp:
        print_roots(coeffs, roots_mp)

    plot_roots(roots_mp, equation, coeffs)


if __name__ == "__main__":
    solve_and_plot()
