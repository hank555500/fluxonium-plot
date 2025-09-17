import sympy as sp
from sympy import pprint
from sympy import simplify
from sympy import factor


x, xp, y, z, zp, t, h = sp.symbols('x xp y z zp t h')

H = sp.Matrix([
    [0,  x,   0],
    [xp, y,   z],
    [0,  zp,  0]
])
i = sp.I
U = sp.exp(i * H * t / h)
simpler_U = U.applyfunc(simplify)

print(H)
pprint(simpler_U[0, 0])