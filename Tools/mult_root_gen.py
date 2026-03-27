from sympy import symbols, expand, prod

x = symbols('x')

n = 10        # number of roots
m = 3          # multiplicity of each root

roots = list(range(1, n + 1))   # [1,2,3,...,n]

# build the polynomial
f = prod((x - r)**m for r in roots)

# expand and extract coefficients
f_expanded = expand(f)
coeffs = f_expanded.as_poly(x).all_coeffs()

# print like Sage
print(" ".join(str(c) for c in coeffs))
