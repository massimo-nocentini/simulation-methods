
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

def generalized_eigenvectors_matrices(M_space, X=IndexedBase('X')):

    Xs = {}

    for i, M_i in M_space.items():
        m = M_i[1].rhs.rows # at least, an eigenvalue should exist so take the first
        X_i = Matrix(m, len(M_i), lambda n, k: M_i[k+1].rhs[n, 0])
        Xs[i] = Eq(X[i], X_i, evaluate=False)

    return Xs


def Jordan_normalform(eigendata, matrices, syms=symbols('X J')):

    data, eigenvals, multiplicities = eigendata
    A, M_space, blocks = matrices
    X, J = syms

    X_cols = sum(len(M_i) for i, M_i in M_space.items())
    Xmatrix = zeros(A.cols, X_cols)

    J_cols = sum(len(J_i) for i, J_i in blocks.items())
    Jmatrix = zeros(X_cols, J_cols)

    row, col = 0, 0
    for i, (λ, m_λ) in data.items():

        M_i = M_space[i]
        J_i = blocks[i]

        assert M_i.keys() == J_i.keys()

        for j in M_i.keys():
            Xmatrix[:, col] = M_i[j].rhs
            Jmatrix[row:row+multiplicities[m_λ], col] = J_i[j]
            col += 1

        row += multiplicities[m_λ]

    Xeq = Eq(X, Xmatrix, evaluate=False)
    Jeq = Eq(J, Jmatrix, evaluate=False)

    return Xeq, Jeq


def Jordan_blocks(eigendata, J=IndexedBase('J')):

    data, eigenvals, multiplicities = eigendata

    blocks = { i: { j: Matrix(m, 1, lambda n,k: λ if n+1 == j else 
                                                1 if n == j else 
                                                0) for j in range(1, m+1) }
               for i, (λ, m_λ) in data.items()
               for m in [multiplicities[m_λ]] # binding
             }

    return blocks


def generalized_eigenvectors_relations(eigendata):

    data, eigenvals, multiplicities = eigendata

    def GER(A, M_space, post=simplify, check=True):

        eqs = {}
        for i, (λ, m_λ) in data.items():

            eig, mul = eigenvals[λ], multiplicities[m_λ]
            Jordan_chain = M_space[i]

            eqs[i] = {}
            for j in range(1, mul):
                x_ij, x_i_succj = Jordan_chain[j].rhs, Jordan_chain[j+1].rhs
                x_ij_eig = x_ij.applyfunc(lambda x: x * λ)
                eq = Eq(A*x_ij, x_ij_eig + x_i_succj)
                eqs[i][j] = Eq(eq.lhs.applyfunc(post), eq.rhs.applyfunc(post), evaluate=False)

            x_i_mul = Jordan_chain[mul].rhs
            x_i_mul_eig = x_i_mul.applyfunc(lambda x: x * λ)
            eq = Eq(A*x_i_mul, x_i_mul_eig)
            eqs[i][mul] = Eq(eq.lhs.applyfunc(post), eq.rhs.applyfunc(post), evaluate=False)

        if check:
            def S(i):
                return i.subs(eigenvals).simplify()

            assert all(lhs.applyfunc(S) == rhs.applyfunc(S)
                       for i, M_i in eqs.items() 
                       for j, eq in M_i.items()
                       for lhs, rhs in [(eq.lhs, eq.rhs)])

        return eqs

    return GER

