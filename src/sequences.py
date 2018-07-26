
from functools import lru_cache
from itertools import product
from collections import namedtuple

from sympy import Symbol, Integer, poly, Eq, expand, zeros, symbols, Matrix, factorial, IndexedBase, solve, simplify

from commons import *

nature = namedtuple('nature', ['is_ordinary', 'is_exponential'])

# the following are indexed bases relative to commonly known RAs
mathcal_R, mathcal_P, mathcal_C, mathcal_S = (
    IndexedBase(r'\mathcal{R}'), # stands for the abstract RA
    IndexedBase(r'\mathcal{P}'), # stands for the Pascal triangle
    IndexedBase(r'\mathcal{C}'), # stands for the Catalan triangle
    IndexedBase(r'\mathcal{S}'), # stands for the Stirling triangle
)

gencoeff_d = IndexedBase('d') # stands for a generic coeff in a matrix


def riordan_matrix_exponential(RA):
    return lambda i,j: factorial(i)*RA(i,j)/factorial(j)

def riordan_matrix_by_convolution(dim, d, h):

    t = symbols('t')                                                # Local symbol to denote a fresh formal variable `t`.

    with lift_to_Lambda(d, return_eq=True) as D:                    # Lift equations `d` and `h` to become callables
      with lift_to_Lambda(h, return_eq=True) as H:                  # objects, returning equations as well.
        d_eq, h_eq = D(t), H(t)                                     # Evaluate both of them to get terms in the var `t`.

    @lru_cache(maxsize=None)
    def column(j):                                                  # columns are memoized for the sake of efficiency
        if not j: return d_eq
        lhs = column(j-1).lhs * h_eq.lhs
        rhs = column(j-1).rhs * h_eq.rhs
        return Eq(lhs, rhs)

    @lru_cache(maxsize=None)
    def C(j):
        return column(j).rhs.series(t, n=dim).removeO()

    return lambda i, j: C(j).coeff(t, i).expand()          # return a lambda to be plugged in a `Matrix` ctor.

def riordan_matrix_by_AZ_sequences(dim, seqs, init={(0,0):1},
                                   ctor=zeros, post=expand,
                                   lattice=None):

    if not lattice: lattice = [(n, k) for n in range(1, dim)
                                      for k in range(n+1)]

    Zseq, Aseq = seqs
    t = symbols('t')

    with lift_to_Lambda(Zseq) as Z, lift_to_Lambda(Aseq) as A:
        Z_series = Z(t).series(t, n=dim).removeO()
        A_series = A(t).series(t, n=dim).removeO()

    R = ctor(dim)

    for cell, i in init.items(): R[cell] = i

    for n, k in lattice:

        if k:
            v = sum(R[n-1, j] * A_series.coeff(t, n=i)
                    for i, j in enumerate(range(k-1, dim)))
        else:
            v = sum(R[n-1, j] * Z_series.coeff(t, n=j)
                    for j in range(dim))

        R[n, k] = v

    if callable(post): R = R.applyfunc(post)

    return lambda n, k: R[n, k]


def riordan_matrix_by_recurrence(dim, rec, init={(0,0):1},
                                 ctor=zeros, post=expand,
                                 lattice=None):

    if not lattice: lattice = [(n, k) for n in range(1, dim)
                                      for k in range(n+1)]

    R = ctor(dim)

    for cell, i in init.items(): R[cell] = i

    for cell in lattice:
        for comb_cell, v in rec(*cell).items():

            try:
                comb = 1 if cell == comb_cell else R[comb_cell]
                combined = v * comb
            except IndexError:
                combined = 0

            R[cell] += combined

    if callable(post): R = R.applyfunc(post)

    return lambda n, k: R[n, k]

def diagonal_matrix(seq, default=0):

    def D(n, k):
        return seq[n] if n == k else default

    return D

def identity_matrix():

    def I(n, k):
        return 1 if n == k else 0

    return I

def frobenius_matrix(seq):

    m = len(seq)
    def F(n, k):
        return 1 if n+1 == k else -seq[-1-k] if n == m-1 else 0

    return F

def unit_vector(i, offset=-1):

    def U(n, k):
        return 1 if n == i+offset else 0

    return U

def Asequence(M):

    A = zeros(M.rows-1, M.cols)
    t = symbols('t')
    current_row_poly = M[0,0]

    for i in range(A.rows):

        next_row_poly = sum(M[i+1,j]*t**(i+1-j) for j in range(i+2))
        div_poly = (next_row_poly/current_row_poly).series(t, n=A.cols).removeO()

        for j in range(A.cols):
            A[i, j] = div_poly.coeff(t, n=j)

        current_row_poly = next_row_poly

    return A

def one(i):
    return 1

def columns_symmetry(M):
    return lambda i, j: M[i, i-j]

def rows_shift_matrix(by):
    return lambda i, j: 1 if i + by == j else 0

def diagonal_func_matrix(f):
    return lambda i, j: f(i) if i == j else 0

def production_matrix(M, exp=False):

    f = factorial if exp else one
    U = Matrix(M.rows, M.cols, rows_shift_matrix(by=1))
    F = Matrix(M.rows, M.cols, diagonal_func_matrix(f))
    F_inv = F**(-1)
    V = F_inv * U * F
    O = F_inv * M * F
    O_inv = O**(-1)
    PM = F * O_inv * V * O * F_inv
    PM = F_inv * PM * F if exp else PM
    PM = PM[:-1, :-1]
    return PM.applyfunc(simplify)


def inspect(M):

    P = production_matrix(M, exp=False)
    C = production_matrix(M, exp=True)

    is_ord = all(P[:1-i, 1] == P[i-1:, i] for i in range(2, P.cols))

    diagonals = { d: [C[j+d,j] for j in range(1,C.rows-d)]
                  for d in range(C.rows-3) }

    is_exp = all(map(is_arithmetic_progression, diagonals.values()))

    return nature(is_ord, is_exp)


def is_arithmetic_progression(prog):

    steps = len(prog)-1
    for _ in range(steps):
        prog = [(b-a).simplify() for a, b in zip(prog, prog[1:])]

    assert len(prog) == 1

    return prog.pop() == 0


def compositional_inverse(h_eq, y=symbols('y'), check=True):

    spec, body = h_eq.lhs, h_eq.rhs
    t, = spec.args

    sols = solve(Eq(y, body), t)
    for sol in sols:
        L = Lambda(y, sol)
        if L(0) == 0:

            if check: assert L(body).simplify() == t

            h_bar = Function(r'\bar{{ {} }}'.format(str(spec.func)))
            eq = Eq(h_bar(y), sol.factor())
            return eq

    raise ValueError

def group_inverse(d_eq, h_eq, post=identity, check=True):

    t, y = symbols('t y')
    g_fn, f_fn = Function('g'), Function('f')

    with lift_to_Lambda(d_eq, return_eq=False) as D:                    # Lift equations `d` and `h` to become callables
      with lift_to_Lambda(h_eq, return_eq=True) as H:                  # objects, returning equations as well.
        f_eq = compositional_inverse(H(t), y, check)
        with lift_to_Lambda(f_eq, return_eq=False) as F:                  # objects, returning equations as well.
            F_t = F(t)
            g = post(1/D(F_t))
            f = post(F_t)

    g_eq, f_eq = Eq(g_fn(t), g.simplify()), Eq(f_fn(t), f.simplify())

    if check:
        with lift_to_Lambda(g_eq, return_eq=False) as G:                    # Lift equations `d` and `h` to become callables
          with lift_to_Lambda(f_eq, return_eq=False) as F:                  # objects, returning equations as well.
            H_rhs = H(t).rhs
            assert (D(t)*G(H_rhs)).simplify() == 1
            assert F(H_rhs).simplify() == t

    return g_eq, f_eq

