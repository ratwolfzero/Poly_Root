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


def companion_diagnostics(C, structural_singular=False):
    # Skip meaningless condition number if singular by construction
    if structural_singular:
        print("INFO: Companion matrix is singular by construction (numerically safe).")
        return None

    try:
        cond_number = mp.norm(C, 1) * mp.norm(C**-1, 1)
    except Exception:
        cond_number = mp.inf

    if cond_number == mp.inf:
        print("WARNING: Companion matrix is singular or numerically unstable.")
    elif cond_number > 1e12:
        print(f"WARNING: Ill-conditioned matrix (cond ≈ {mp.nstr(cond_number, 4)})")
    elif cond_number > 1e6:
        print(f"NOTICE: Moderately ill-conditioned matrix (cond ≈ {mp.nstr(cond_number, 4)})")

    return cond_number


def detect_clusters(roots, tol=mpf('1e-8')):
    """Relative-tolerance cluster detection (scaled to each root's magnitude)."""
    clusters = []
    used = [False] * len(roots)
    for i, r1 in enumerate(roots):
        if used[i]:
            continue
        cluster = [r1]
        used[i] = True
        mag = max(abs(r1), mpf(1))
        for j, r2 in enumerate(roots):
            if not used[j] and abs(r1 - r2) < tol * mag:
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

    has_zero_root = (coeffs[-1] == 0)

    if has_zero_root:
        print("NOTICE: Polynomial has root at x = 0")

    a = [c / leading for c in coeffs[1:]]
    n = len(a)

    if n == 1:
        return None, [-a[0]], has_zero_root

    C = matrix(n)
    for row in range(1, n):
        C[row, row - 1] = mpf(1)

    for row in range(n):
        C[row, n - 1] = -a[n - 1 - row]

    return C, None, has_zero_root

# ----------------------------- Root Computation ----------------------------- #


def compute_roots(coeffs):
    coefficient_diagnostics(coeffs)

    C, linear_root, has_zero_root = build_companion(coeffs)

    if linear_root is not None:
        roots_mp = [mpc(linear_root[0])]
    else:
        companion_diagnostics(C, structural_singular=has_zero_root)
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
        print("WARNING: Large absolute residuals")

    expected_res = mpf(10) ** (-mp.dps // 2 + 5)
    if max(rel_residuals, default=0) > expected_res:
        print(f"WARNING: Large relative residuals (expected ~1e{-mp.dps//2})")

    if detect_clusters(roots_mp):
        print("NOTICE: Clustered roots detected → high sensitivity likely.")

    roots_mp.sort(key=lambda z: (
        mp.re(z),
        abs(z),
        mp.atan2(mp.im(z), mp.re(z))
    ))

    return roots_mp

# ----------------------------- Display ----------------------------- #


def root_tolerance(roots):
    """Compute consistent tolerance based on root magnitudes."""
    if not roots:
        return mpf('1e-10')
    max_mag = max(abs(r) for r in roots)
    return mpf('1e-10') * max(mpf(1), max_mag)


def print_roots(coeffs, roots_mp):
    if not roots_mp:
        print("No roots to display.")
        return

    # --- Human-friendly precision ---
    digits = min(6, max(10, int(0.4 * mp.dps)))

    tol = root_tolerance(roots_mp)

    # ---------- Normalize roots (snap tiny imaginary parts) ----------
    normalized = []
    for r in roots_mp:
        imag = mp.im(r)
        if abs(imag) < tol:
            imag = mpf(0)
        normalized.append(mpc(mp.re(r), imag))

    # ---------- Separate real and complex ----------
    real_roots = [r for r in normalized if mp.im(r) == 0]
    complex_roots = [r for r in normalized if mp.im(r) != 0]

    # ---------- Group conjugate pairs ----------
    used = [False] * len(complex_roots)
    pairs = []
    leftovers = []

    for i, r1 in enumerate(complex_roots):
        if used[i]:
            continue
        conj = mpc(mp.re(r1), -mp.im(r1))
        found = False

        for j, r2 in enumerate(complex_roots):
            if i != j and not used[j]:
                if abs(r2 - conj) < tol:
                    pairs.append((r1, r2))
                    used[i] = used[j] = True
                    found = True
                    break

        if not found:
            leftovers.append(r1)
            used[i] = True

    # ---------- Formatting ----------
    def format_root(r):
        rel_res = relative_residual(coeffs, r)
        real_str = mp.nstr(mp.re(r), digits)
        imag_val = mp.im(r)
        imag_str = mp.nstr(imag_val, digits)
        if imag_val >= 0:
            imag_str = "+" + imag_str
        mag_str = mp.nstr(abs(r), digits)
        res_str = mp.nstr(rel_res, 6)
        return real_str, imag_str, mag_str, res_str

    rows = []

    # --- Real roots first ---
    for r in sorted(real_roots, key=lambda z: mp.re(z)):
        rows.append(("real", format_root(r)))

    # --- Complex conjugate pairs ---
    for r1, r2 in pairs:
        rows.append(("pair", format_root(r1)))
        rows.append(("pair", format_root(r2)))

    # --- Any leftovers (rare, non-conjugate numerical artifacts) ---
    for r in leftovers:
        rows.append(("complex", format_root(r)))

    # ---------- Column widths ----------
    w_real = max(len(r[1][0]) for r in rows) + 2
    w_imag = max(len(r[1][1]) for r in rows) + 3
    w_mag = max(len(r[1][2]) for r in rows) + 2

    total_width = 10 + w_real + w_imag + w_mag + 20

    # ---------- Header ----------
    print("-" * total_width)
    print(
        f"{'Root #':<6} "
        f"{'Real':>{w_real}} "
        f"{'Imaginary':>{w_imag}} "
        f"{'|z|':>{w_mag}} "
        f"{'Rel.Residual':>16}"
    )
    print("-" * total_width)

    # ---------- Rows ----------
    idx = 1
    for kind, (real_str, imag_str, mag_str, res_str) in rows:
        print(
            f"{idx:<6d} "
            f"{real_str:>{w_real}} "
            f"{imag_str:>{w_imag}}j "
            f"{mag_str:>{w_mag}} "
            f"{res_str:>16}"
        )
        idx += 1
# ----------------------------- Combined Plot ----------------------------- #

def plot_complex_plane(ax, roots_mp):
    roots_np = np.array([
        complex(float(mp.re(z)), float(mp.im(z)))
        for z in roots_mp
    ])

    ax.axhline(0, lw=1)
    ax.axvline(0, lw=1)
    ax.scatter(roots_np.real, roots_np.imag, color="red",
                s=10, zorder=5, label="Roots")
    t = np.linspace(0, 2*np.pi, 200)
    ax.plot(np.cos(t), np.sin(t), ls="--", alpha=0.5,
             color="gray", label="Unit Circle")

    max_modulus = np.max(np.abs(roots_np)) if roots_np.size > 0 else 0
    view_radius = 1.1 * max(max_modulus, 1.0)

    ax.set_xlim(-view_radius, view_radius)
    ax.set_ylim(-view_radius, view_radius)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_title("Roots in Complex Plane")
    ax.set_xlabel("Real")
    ax.set_ylabel("Imaginary")
    ax.legend(loc="best")
    ax.grid(True, linestyle=":", alpha=0.6)


def plot_polynomial_curve(ax, coeffs, roots_mp):
    # First, compute real roots (needed for insertion and markers)
    if roots_mp:
        tol = root_tolerance(roots_mp)
        real_roots = [
            float(mp.re(r)) for r in roots_mp if abs(mp.im(r)) < tol
        ]
    else:
        real_roots = []

    # Build x values including real roots
    real_parts = [float(mp.re(r)) for r in roots_mp]
    spread = max(real_parts) - min(real_parts) if real_parts else 1.0
    x_center = sum(real_parts) / len(real_parts) if real_parts else 0.0
    x_pad = 1.5 * max(spread, 1.0)

    # Initial sampling (use 2000 for speed; increase if needed)
    x_vals_initial = np.linspace(x_center - x_pad, x_center + x_pad, 2000)

    # Insert each real root (if not already present)
    x_list = x_vals_initial.tolist()
    for root in real_roots:
        if not any(abs(root - x) < 1e-12 for x in x_list):
            x_list.append(root)
    x_list.sort()
    x_vals = np.array(x_list)

    # Evaluate polynomial at each x (high precision)
    y_vals = np.array([float(mp.re(poly_eval(coeffs, mpf(xx))))
                      for xx in x_vals])

    # Prevent overflow
    max_abs = np.max(np.abs(y_vals)) if len(y_vals) > 0 else 1.0
    if np.isinf(max_abs) or max_abs > 1e300:
        print("NOTICE: Polynomial values exceed float range → clipping plot")
        y_vals = np.clip(y_vals, -1e300, 1e300)
        max_abs = 1e300

    # Plot curve and markers
    ax.plot(x_vals, y_vals, lw=1, label=f"p(x)")
    ax.axhline(0, lw=1)
    ax.axvline(0, lw=1)

    if real_roots:
        ax.scatter(real_roots, [0]*len(real_roots),
                    color="blue", s=10, zorder=5, label="Real Roots")

    # Initial scale setting
    if max_abs > 1e6:
        linthresh = 1.0
        ax.set_yscale('symlog', linthresh=linthresh, linscale=1.0)
        print(f"NOTICE: symlog y-scale activated (linthresh = {linthresh})")
    else:
        ax.set_ylim(-max_abs * 1.05, max_abs * 1.05)

    ax.set_title("Polynomial Curve and Roots in Real Domain")
    ax.set_xlabel("x")
    ax.set_ylabel("$f(x)$")
    ax.legend(loc="best")
    ax.grid(True, linestyle=":", alpha=0.7)

    return max_abs


def add_interactive_scale(fig, ax, max_abs):
    def on_key(event):
        if event.key == 'l':
            ax.set_yscale('linear')
            ax.set_ylim(-max_abs * 1.05, max_abs * 1.05)
            fig.canvas.draw_idle()
            print("Switched to linear scale")
        elif event.key == 'y':
            ax.set_yscale('symlog', linthresh=1.0, linscale=1.0)
            fig.canvas.draw_idle()
            print("Switched to symlog scale")
    fig.canvas.mpl_connect('key_press_event', on_key)


def plot_combined(coeffs, roots_mp, equation):
    if not roots_mp:
        return

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(24, 7),
        gridspec_kw={'width_ratios': [1, 1.5]}
    )
    fig.canvas.manager.set_window_title(
        f"Complex Plane and Polynomial Curve")
    fig.suptitle(f"Polynomial: {equation}", wrap=True)

    plot_complex_plane(ax1, roots_mp)
    max_abs = plot_polynomial_curve(ax2, coeffs, roots_mp)

    plt.subplots_adjust(top=0.85)

    add_interactive_scale(fig, ax2, max_abs)

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
