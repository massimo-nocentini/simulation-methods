
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
    if subs_eigenvals: 
        g = define(g.lhs, g.rhs.subs(eigenvals, simultaneous=True), ctor=MFEq)

    return f, g

power = partial(make, f=define(function('P', z), z**r, ctor=FEq))
inverse = partial(make, f=define(function('I', z), 1/z, ctor=FEq))
square_root = partial(make, f=define(function('R', z), sympy.sqrt(z), ctor=FEq))
exp = partial(make, f=define(function('E', z), sympy.exp(alpha*z), ctor=FEq))
log = partial(make, f=define(function('L', z), sympy.log(z), ctor=FEq))
sin = partial(make, f=define(function('S', z), sympy.sin(z), ctor=FEq))
cos = partial(make, f=define(function('C', z), sympy.cos(z), ctor=FEq))
mobius = partial(make, f=define(function('M', z), (a*z+b)/(c*z+d), ctor=FEq))
normal = partial(make, f=define(function('N', z), sympy.exp(- z**2/2)/sympy.sqrt(2*pi), ctor=FEq))
