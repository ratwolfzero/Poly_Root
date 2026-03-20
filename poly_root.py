from mpmath import mp, mpf, mpc, matrix, eig
import matplotlib.pyplot as plt
import numpy as np


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
        abs_c = abs(c)
        coeff_str = "" if deg > 0 and abs_c == 1 else fmt(abs_c)
        var_str = f"x^{deg}" if deg > 1 else ("x" if deg == 1 else "")
        term_str = f"{coeff_str}{var_str}"

        if not terms:  # first term
            sign = "-" if c < 0 else ""
            terms.append(f"{sign}{term_str}")
        else:
            sign = " + " if c > 0 else " - "
            terms.append(f"{sign}{term_str}")

    if not terms:
        return "0 = 0"
    return "".join(terms) + " = 0"


def poly_eval(coeffs, x):
    p = mpc(0)
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
            print("Likely due to extreme ill-conditioning.")
            print(f"{type(e).__name__}: {e}")
            return []

    abs_residuals = [abs(poly_eval(coeffs, r)) for r in roots_mp]
    rel_residuals = [relative_residual(coeffs, r) for r in roots_mp]

    if max(abs_residuals, default=0) > 1:
        print(f"WARNING: Large absolute residuals")

    if max(rel_residuals, default=0) > mpf(1e-6):
        print(f"WARNING: Large relative residuals")

    if detect_clusters(roots_mp):
        print("NOTICE: Clustered roots detected → high sensitivity likely.")

    roots_mp.sort(key=lambda z: (mp.re(z), mp.im(z)))
    return roots_mp


# ----------------------------- Display ----------------------------- #
def print_roots(coeffs, roots_mp):
    # Header - Imag column widened slightly to account for the 'j'
    print("-" * 75)
    print(f"{'Root #':<6} {'Real':>14} {'Imaginary':>16} {'|z|':>10} {'Residual':>14}")
    print("-" * 75)

    for i, r in enumerate(roots_mp, 1):
        residual = abs(poly_eval(coeffs, r))
        r_real = float(mp.re(r))
        r_imag = float(mp.im(r))
        mag = float(abs(r))

        # Format strings:
        # :>14.6f  -> Right-aligned, 14 chars wide, 6 decimal places
        # :+15.6f  -> Right-aligned, 15 chars wide, 6 decimals, ALWAYS shows +/-
        print(
            f"{i:<6d} "
            f"{r_real:14.9f} "
            # The '+' forces the sign and keeps it in the column
            f"{r_imag:+15.9f}j "
            f"{mag:10.6f} "
            f"{mp.nstr(residual, 9):>14}"
        )


# ----------------------------- Combined Plot ----------------------------- #
def plot_combined(coeffs, roots_mp, equation):
    if not roots_mp:
        return

    roots_np = np.array([
        complex(float(mp.re(z)), float(mp.im(z)))
        for z in roots_mp
    ])

    # Increase the top margin slightly to make room for the Suptitle
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(16, 8),
        gridspec_kw={'width_ratios': [1, 1.5]}
    )

    fig.suptitle(f"Polynomial Equation: {equation}",
                 fontsize=10, fontweight='bold', wrap=True)
    fig.canvas.manager.set_window_title(f"Solver Output")

    # ---- Complex Plane ----
    ax1.axhline(0, color='black', lw=1, alpha=0.5)
    ax1.axvline(0, color='black', lw=1, alpha=0.5)
    ax1.scatter(roots_np.real, roots_np.imag, color="red",
                s=30, zorder=5, label="Roots")

    t = np.linspace(0, 2*np.pi, 200)
    ax1.plot(np.cos(t), np.sin(t), ls="--", alpha=0.5,
             color="gray", label="Unit Circle")

    ax1.set_title("Roots in Complex Plane")
    ax1.set_xlabel("Real")
    ax1.set_ylabel("Imaginary")
    ax1.legend(loc="upper right")
    ax1.grid(True, linestyle=":", alpha=0.6)

    # Square scaling for the complex plane
    limit = max(np.max(np.abs(roots_np.real)) if len(roots_np) > 0 else 1,
                np.max(np.abs(roots_np.imag)) if len(roots_np) > 0 else 1, 1.1)
    ax1.set_xlim(-limit*1.1, limit*1.1)
    ax1.set_ylim(-limit*1.1, limit*1.1)
    ax1.set_aspect('equal')

    # ---- Polynomial Curve ----
    real_parts = [float(mp.re(r)) for r in roots_mp]
    x_min, x_max = min(real_parts), max(real_parts)
    margin = max((x_max - x_min) * 0.5, 2.0)
    x_vals = np.linspace(x_min - margin, x_max + margin, 1000)

    y_vals = [float(mp.re(poly_eval(coeffs, mpf(x)))) for x in x_vals]

    ax2.plot(x_vals, y_vals, color='tab:blue', lw=2, label="p(x) [Real]")
    ax2.axhline(0, color='black', lw=1)

    # Highlight real roots on the curve
    tol = 1e-9
    real_roots = [r.real for r in roots_np if abs(r.imag) < tol]
    if real_roots:
        ax2.scatter(real_roots, [0]*len(real_roots),
                    color="blue", s=40, zorder=5, label="Real Roots")

    ax2.set_title("Polynomial Curve (Real Domain)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("p(x)")
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend()

    # Dynamic Y-scaling to keep the curve visible
    # Ignore extreme outliers for scaling
    y_limit = np.percentile(np.abs(y_vals), 95)
    if y_limit > 1e6:
        ax2.set_yscale('symlog', linthresh=100)
    else:
        ax2.set_ylim(-y_limit*2, y_limit*2)

    # Adjust layout to make room for suptitle
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


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

    plot_combined(coeffs, roots_mp, equation)


if __name__ == "__main__":
    solve_and_plot()
