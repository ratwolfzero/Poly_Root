from math import prod

# -----------------------------
# INPUT: roots + multiplicities
# Format: (root, multiplicity)
# -----------------------------
roots_with_mult = [
    (-1, 5),   # (x - 1)^2
    (0, 33),   # (x - 2)^3
    (1, 5)    # (x - 5)^1
]

# Optional: leading coefficient
leading_coeff = 1


# -----------------------------
# POLYNOMIAL MULTIPLICATION
# -----------------------------
def multiply_poly(p, q):
    """Multiply two polynomials (lists of coefficients)."""
    res = [0] * (len(p) + len(q) - 1)
    for i in range(len(p)):
        for j in range(len(q)):
            res[i + j] += p[i] * q[j]
    return res


# -----------------------------
# BUILD POLYNOMIAL
# -----------------------------
def build_polynomial(roots_with_mult, leading_coeff=1):
    poly = [leading_coeff]  # represents constant

    for r, m in roots_with_mult:
        factor = [1]  # will become (x - r)^m
        base = [1, -r]  # (x - r)

        for _ in range(m):
            factor = multiply_poly(factor, base)

        poly = multiply_poly(poly, factor)

    return poly


# -----------------------------
# PRETTY PRINT
# -----------------------------
def format_polynomial(coeffs):
    """Convert coefficient list to readable polynomial string."""
    terms = []
    degree = len(coeffs) - 1

    for i, c in enumerate(coeffs):
        power = degree - i
        if c == 0:
            continue

        # sign
        sign = "+" if c > 0 else "-"
        c_abs = abs(c)

        if power == 0:
            term = f"{c_abs}"
        elif power == 1:
            term = f"{'' if c_abs == 1 else c_abs}x"
        else:
            term = f"{'' if c_abs == 1 else c_abs}x^{power}"

        terms.append((sign, term))

    # build string
    if not terms:
        return "0"

    first_sign, first_term = terms[0]
    result = ("" if first_sign == "+" else "-") + first_term

    for sign, term in terms[1:]:
        result += f" {sign} {term}"

    return result


# -----------------------------
# RUN
# -----------------------------
coeffs = build_polynomial(roots_with_mult, leading_coeff)

print("Polynomial:")
print("f(x) =", format_polynomial(coeffs))

print("\nCoefficients (highest degree → constant):")
print(" ".join(str(c) for c in coeffs))