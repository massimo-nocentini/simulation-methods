
from functools import partial

import sympy

from sympy import symbols, Function, IndexedBase
from sympy.abc import a, b, c, d, z, r # latin symbols
from sympy.abc import alpha, beta, pi # greek symbols

from commons import *
from sequences import *
from matrix_functions import *

def function(name, z):
    F = Function(name)
    return F(z)

def make(eigendata, Phi_polynomials, subs_eigenvals=True, *, f):

    data, eigenvals, multiplicities = eigendata.rhs # unpacking to use `eigenvals` in `subs`

    g = Hermite_interpolation_polynomial(f, eigendata, Phi_polynomials)
    g = g.subs(eigenvals, simultaneous=True) if subs_eigenvals else g
    with lift_to_matrix_function(g) as G:
        return f, g, G

power = partial(make, f=define(function('P', z), z**r))
inverse = partial(make, f=define(function('I', z), 1/z))
square_root = partial(make, f=define(function('R', z), sympy.sqrt(z)))
exp = partial(make, f=define(function('E', z), sympy.exp(alpha*z)))
log = partial(make, f=define(function('L', z), sympy.log(z)))
sin = partial(make, f=define(function('S', z), sympy.sin(z)))
cos = partial(make, f=define(function('C', z), sympy.cos(z)))
mobius = partial(make, f=define(function('M', z), (a*z+b)/(c*z+d)))
normal = partial(make, f=define(function('N', z), sympy.exp(- z**2/2)/sympy.sqrt(2*pi))
