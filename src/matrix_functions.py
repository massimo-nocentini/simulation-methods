
'''
    Start of :code:`matrix_functions` module documentation

    >>> from sympy import *

    >>> from commons import *
    >>> from sequences import *
    >>> from matrix_functions import *

'''

from contextlib import contextmanager
from functools import lru_cache

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
def lift_to_matrix_function(g_def):
    
    def lift(matrix_def, post=lambda c: c.simplify()):
       
        matrix = matrix_def.rhs
        z, *_ = g_def.lhs.args
        gpp, *_ = poly_from_expr(g_def.rhs, z)
        coeffs = dict(gpp.terms())

        I = eye(matrix.rows)
        Z = coeffs.get((gpp.degree(),), 0)*I # pay attention to tuple `(gpp.degree(),)`

        for d in range(gpp.degree()-1, -1, -1): 
            Z = Z*matrix + coeffs.get((d,), 0)*I

        lhs = g_def.lhs.func(matrix_def.lhs)
        rhs = Z.applyfunc(post) if callable(post) else Z
        return Eq(lhs, rhs, evaluate=False)

    yield lift

def spectrum(matrix_def):

    matrix_name, matrix = matrix_def.lhs, matrix_def.rhs

    data = {}
    eigenvals = {}
    multiplicities = {}
    EVs = matrix.eigenvals() # in this call lies the complexity
    for i, (eigen_value, multiplicity) in enumerate(EVs.items(), start=1):
        lamda, mul = lamda_indexed[i], mul_indexed[i]
        data[i] = lamda, mul
        eigenvals[lamda] = eigen_value
        multiplicities[mul] = multiplicity

    eigendata = data, eigenvals, multiplicities

    return Eq(Function(r'\sigma').__call__(matrix_name), eigendata, evaluate=False)


def Phi_poly_ctor(deg):
    '''
    Ctor.

    >>> p = Phi_poly_ctor(deg=3)
    >>> p
    Eq(\Phi(z, i, j), z**3*\phi[i, j, 0] + z**2*\phi[i, j, 1] + z*\phi[i, j, 2] + \phi[i, j, 3])

    .. testcode::
        :hide:

        save_latex_repr(p, './source/latex-snippets/Phi_poly_ctor-0.rst')

    .. include:: latex-snippets/Phi_poly_ctor-0.rst

    '''

    Phi = Function(r'\Phi')
    z, i, j, k = symbols('z, i, j, k')

    terms = Sum(phi_abstract_coefficient[i, j, deg-k] * z**k, (k, 0, deg)).doit()

    return Eq(Phi(z, i, j), terms)


def component_polynomials(eigendata, early_eigenvals_subs=False, verbose_return=False):

    data, eigenvals, multiplicities = eigendata.rhs

    deg = sum(multiplicities.values()) - 1

    Phi_poly = Phi_poly_ctor(deg)
    z, *_ = Phi_poly.lhs.args

    polynomials = {}

    for i, (lamda, mul) in data.items():

        for j in range(1, multiplicities[mul] + 1):

            with lift_to_Lambda(Phi_poly) as Phi: 
                Phi_ij = Phi(z, i, j)

                @lru_cache(maxsize=None)
                def Phi_ij_derivative(j):
                    return Phi_ij_derivative(j-1).diff(z) if j else Phi_ij

            eqs = set()
            for l, (lamda_l, mul_l) in data.items():

                λ_l, m_l = eigenvals[lamda_l] , multiplicities[mul_l]

                derivatives = {r: Eq(partial_r(z), Phi_ij_derivative(r-1))
                               for r in range(1, m_l + 1)
                               for partial_r in [Function(r'\partial^{{({})}}'.format(r))]}

                for r, der_eq in derivatives.items():
                    with lift_to_Lambda(der_eq) as der_fn:
                        lhs = der_fn(λ_l if early_eigenvals_subs else lamda_l)
                        rhs = KroneckerDelta(l, i) * KroneckerDelta(r, j)
                        eqs.add(Eq(lhs, rhs))

            respect_to = [phi_abstract_coefficient[i,j,k]
                          for k in range(poly(Phi_ij, z).degree()+1)]

            sols = solve(eqs, respect_to)

            lhs = Function(r'\Phi_{{ {}, {} }}'.format(i, j)).__call__(z)
            rhs = Phi_ij.subs(sols, simultaneous=True).collect(z)
            baked_poly = Eq(lhs, rhs)
            formal_poly = Eq(lhs, Phi_ij)

            polynomials[i,j] = (baked_poly, formal_poly, eqs, sols) if verbose_return else baked_poly

    return polynomials

def component_polynomials_riordan(degree):
    '''
    Quick generation of component polynomials to be used as bases in functions applications to Riordan arrays.
    

    >>> m = 8
    >>> R_cal, d = IndexedBase(r'\mathcal{R}'), IndexedBase('d') # helpers bases
    >>> R = define(R_cal[m], Matrix(m, m, riordan_matrix_by_recurrence(m, lambda n, k: {(n, k):1 if n == k else d[n, k]})))
    >>> eigendata = spectrum(R)
    >>> assert (component_polynomials(eigendata, early_eigenvals_subs=True) == 
    ...         component_polynomials_riordan(degree=m-1))

    '''
    z = symbols('z')
    return {(1, j): define( let=Function(r'\Phi_{{ {}, {} }}'.format(1, j))(z), 
                            be=Sum((-1)**(j-1-k)*z**k/(factorial(k)*factorial(j-1-k)), (k, 0, j-1)).doit()) 
                for j in range(1, degree+2)}


def Hermite_interpolation_polynomial(f_eq, eigendata, Phi_polys, matrix_form=False):

    data, eigenvals, multiplicities = eigendata.rhs
    delta = Function(r'\delta')
    g = Function('{}_{{ {} }}'.format(f_eq.lhs.func, sum(multiplicities.values())))
    Z = IndexedBase('Z')
    z, *_ = f_eq.lhs.args

    g_poly = Integer(0)
    
    with lift_to_Lambda(f_eq) as f_fn:

        @lru_cache(maxsize=None)
        def F_derivative(j):
            return F_derivative(j-1).diff(z) if j else f_fn(z)  

        for i, (lamda, mul) in data.items():
            for j in range(1, multiplicities[mul]+1):
                 with lift_to_Lambda(Eq(delta(z), F_derivative(j-1))) as der_fn,\
                      lift_to_Lambda(Phi_polys[i, j]) as Phi:
                        derivative_term = der_fn(lamda)
                        Phi_term = Z[i,j] if matrix_form else Phi(z)
                        g_poly += derivative_term*Phi_term

    return Eq(g(z), g_poly.expand().collect(z))


def component_matrices( matrix_eq, 
                        Phi_polys, 
                        Zctor=lambda lhs: IndexedBase(r'Z^{{\left[ {} \right]}}'.format(lhs)), 
                        post=factor):
    
    matrix = matrix_eq.rhs
    Z = Zctor(matrix_eq.lhs)
    Z_matrices = {}
    
    for (i, j), cp in Phi_polys.items():
        with lift_to_matrix_function(cp) as cp_fn:
            Z_ij = cp_fn(matrix).applyfunc(post)
            Z_matrices[i, j] = Eq(Z[i, j], Z_ij , evaluate=False)
    
    return Z_matrices


def M_space(cmatrices, x=IndexedBase(r'\boldsymbol{x}')):

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

    data, eigenvals, multiplicities = eigendata.rhs
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

    data, eigenvals, multiplicities = eigendata.rhs

    blocks = { i: { j: Matrix(m, 1, lambda n,k: λ if n+1 == j else 
                                                1 if n == j else 
                                                0) for j in range(1, m+1) }
               for i, (λ, m_λ) in data.items()
               for m in [multiplicities[m_λ]] # binding
             }

    return blocks


def generalized_eigenvectors_relations(eigendata):

    data, eigenvals, multiplicities = eigendata.rhs

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
            eq = Eq(A*x_i_mul, x_i_mul_eig, evaluate=False)
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

def split_X_matrix(X, v, factor=False, normalization=one):

    split = []
    factorized = []
    for i in range(X.rows):
        indicator_i = {v[j]:v[i] if j == i else 0 for j in range(X.rows)}
        X_i = Matrix(X.rows, X.cols, lambda n, k: X[n,k].subs(indicator_i))
        X_i_factorized = Mul(Matrix(X.rows, X.cols, lambda n, k: X_i[n,k]/(v[i]*normalization(k))),
                             Matrix(X.rows, X.cols, lambda n, k: v[i]*normalization(k) if n==k else 0),
                             evaluate=False)

        assert X_i == X_i_factorized.args[0]*X_i_factorized.args[1]

        split.append(X_i)
        factorized.append(X_i_factorized)

    assert X == sum(split, zeros(X.rows, X.cols))

    return factorized if factor else split 


