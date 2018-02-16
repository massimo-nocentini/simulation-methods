
from contextlib import contextmanager

from sympy import *
from sympy.abc import n, i, N, x, lamda, phi, z, j, r, k, a

from commons import *
from sequences import *

#{{{ generic indexed bases according Brugnano notation
phi_abstract_coefficient = IndexedBase(r'\phi')
lamda_indexed = IndexedBase(r'\lambda')
mul_indexed = IndexedBase(r'm')
#}}}


@contextmanager
def lift_to_matrix_function(f_def):
    
    def lift(matrix, post=lambda c: c.simplify()):
       
        z, *rest = f_def.lhs.args
        gpp, *_ = poly_from_expr(f_def.rhs, z)

        I = eye(matrix.rows)
        Z = f_def.rhs.coeff(z, gpp.degree())*I

        for d in range(gpp.degree()-1, -1, -1): 
            Z = Z*matrix + f_def.rhs.coeff(z, d)*I

        return Z.applyfunc(post) if callable(post) else Z

    yield lift

def Phi_poly_ctor(deg, z=Symbol('z'), i=Symbol('i'), j=Symbol('j')):
    k = symbols('k')
    Phi = Function(r'\Phi')
    terms = Sum(phi_abstract_coefficient[i, j, deg-k] * z**k, (k, 0, deg)).doit()
    return Eq(Phi(z, i, j), terms)

def eigen_data(matrix):
    data = {}
    eigenvals = {}
    multiplicities = {}
    for i, (eigen_value, multiplicity) in enumerate(matrix.eigenvals().items(), start=1):
        lamda, mul = lamda_indexed[i], mul_indexed[i]
        data[i] = lamda, mul
        eigenvals[lamda] = eigen_value
        multiplicities[mul] = multiplicity

    return data, eigenvals, multiplicities

def Phi_poly_define(Phi_poly, eigendata):
    
    data, eigenvals, multiplicities = eigendata
    z, *_ = Phi_poly.lhs.args
    
    def maker(i, j, verbose_return=False):

        Phi_def = Function(r'\Phi_{{ {}, {} }}'.format(i, j))

        with lift_to_Lambda(Phi_poly) as Phi: 
            Phi_ij = Phi(z, i, j)

        eqs = set()
        for l, (lamda_l, m_l) in data.items():

            derivatives = {r:Eq(delta_r(z), Phi_ij.diff(z, r-1))
                           for r in range(1, multiplicities[m_l]+1)
                           for delta_r in [Function(r'\delta_{{{}}}'.format(r))]}

            for r, der in derivatives.items():
                with lift_to_Lambda(der) as der_fn:
                    eqs.add(Eq(der_fn(lamda_l), KroneckerDelta(l, i)*KroneckerDelta(r, j)))

        respect_to = [phi_abstract_coefficient[i,j,k] 
                      for k in range(poly(Phi_ij, z).degree()+1)]

        sols = solve(eqs, respect_to)

        lhs = Phi_def(z)
        rhs = Phi_ij.subs(sols, simultaneous=True).collect(z)
        baked_poly = Eq(lhs, rhs)
        return (baked_poly, Eq(lhs, Phi_ij), eqs, sols) if verbose_return else baked_poly
    
    return maker
        
def component_polynomials(Phi_poly, eigendata):
    
    make = Phi_poly_define(Phi_poly, eigendata)
    data, eigenvals, multiplicities = eigendata
    polynomials = {(i, j): make(i, j)
                   for i, (lamda, mul) in data.items()
                   for j in range(1, multiplicities[mul]+1)}

    return polynomials

def component_polynomials_riordan(degree):
    z = symbols('z')
    return {(1, j): define( let=Function(r'\Phi_{{ {}, {} }}'.format(1, j))(z), 
                            be=Sum((-1)**(j-1-k)*z**k/(factorial(k)*factorial(j-1-k)), (k, 0, j-1)).doit()) 
                for j in range(1, degree+1)}


def g_poly(f_eq, eigendata, Phi_polys, matrix_form=False):
    
    data, eigenvals, multiplicities = eigendata
    delta, g = Function(r'\delta'), Function('g')
    Z = IndexedBase('Z')
    z, *rest = f_eq.lhs.args
    
    g_poly = Integer(0)
    
    for i, (lamda, mul) in data.items():
        for j in range(1, multiplicities[mul]+1):
            with lift_to_Lambda(f_eq) as f_fn,\
                 lift_to_Lambda(Eq(delta(z), f_fn(z).diff(z, j-1))) as der_fn,\
                 lift_to_Lambda(Phi_polys[i, j]) as Phi:
                    derivative_term = der_fn(lamda)
                    Phi_term = Z[i,j] if matrix_form else Phi(z)
                    g_poly += derivative_term*Phi_term
    
    return Eq(g(z), g_poly.expand().collect(z))


def component_matrices(matrix, Phi_polys, Z=IndexedBase('Z')):
    
    Z_matrices = {}
    
    for (i, j), cp in Phi_polys.items():
        with lift_to_matrix_function(cp) as cp_fn:
            Z_matrices[i, j] = Eq(Z[i, j], cp_fn(matrix), evaluate=False)
    
    return Z_matrices


def M_space(cmatrices, x=IndexedBase('x')):

    def M_S(v):

        M_space = {}

        for (i, j), Z_eq in cmatrices.items():

            m = Z_eq.rhs.rows

            if i not in M_space: 
                M_space[i] = {}

            Z_i1 = cmatrices[i,1].rhs 
            Z_i2 = cmatrices[i,2].rhs if (i,2) in cmatrices else zeros(m, m)

            x_ij = Z_i2**(j-1) * Z_i1 * v

            M_space[i][j] = Eq(x[i,j], x_ij, evaluate=False)

        return M_space

    return M_S






