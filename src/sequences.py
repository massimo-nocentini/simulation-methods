
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

    t = symbols('t')                                                # Local symbol to denote a formal variable `t`.

    with lift_to_Lambda(d, return_eq=True) as D:                    # Lift equations `d` and `h` to become callables
      with lift_to_Lambda(h, return_eq=True) as H:                  # objects, returning equations as well, in order to
        d_eq, h_eq = D(t), H(t)                                     # let both of them depend on variable `t`.

    @lru_cache(maxsize=None)
    def column(j):                                                  # Columns are memoized for the sake of efficiency.
        if not j: return d_eq                                       # Base case.
        lhs = column(j-1).lhs * h_eq.lhs                            # Otherwise, use already computed column to build
        rhs = column(j-1).rhs * h_eq.rhs                            # the current `lhs` and `rhs`.
        return Eq(lhs, rhs)

    @lru_cache(maxsize=None)
    def C(j):                                                       # Local function that performs Taylor expansion
        return column(j).rhs.series(t, n=dim).removeO()             # of columns, which are symbolic terms up to now.

    return lambda i, j: C(j).coeff(t, i).expand()                   # Return a lambda to be plugged in a `Matrix` ctor.

def riordan_matrix_by_AZ_sequences(dim, seqs, init={(0,0):1},
                                   ctor=zeros, post=expand,
                                   lattice=None):

    if not lattice: lattice = [(n, k) for n in range(1, dim)        # `lattice` plays the same role as in function
                                      for k in range(n+1)]          # `riordan_matrix_by_recurrence`, see that.

    Zseq, Aseq = seqs                                               # Destruct `seqs` into objs denoting Z and A seqs.
    t = symbols('t')                                                # Local formal variable `t`.

    with lift_to_Lambda(Zseq) as Z, lift_to_Lambda(Aseq) as A:      # Promote both sequences as callable objects to
        Z_series = Z(t).series(t, n=dim).removeO()                  # evaluate both of them in `t` and then perform
        A_series = A(t).series(t, n=dim).removeO()                  # two Taylor expansions with respect to `t`.

    R = ctor(dim)                                                   # Init the matrix.

    for cell, i in init.items(): R[cell] = i                        # Set buondary values for combinations.

    for n, k in lattice:                                            # Visit each cell according to the order `lattice`

        if k:                                                       # If the cell lies not on the first column,
            v = sum(R[n-1, j] * A_series.coeff(t, n=i)              # combine using `A_series` coefficients
                    for i, j in enumerate(range(k-1, dim)))         # up to the cell in position `(k-1, dim-1)`.
        else:                                                       # Otherwise, the cell lies on the first column,
            v = sum(R[n-1, j] * Z_series.coeff(t, n=j)              # combine using `Z_series` coefficients
                    for j in range(dim))                            # up to the cell in position `(k-1, dim-1)`.

        R[n, k] = v                                                 # Store the combination.

    if callable(post): R = R.applyfunc(post)                        # If some post processing is desired, do it.

    return lambda n, k: R[n, k]                                     # Return a lambda that uses the computed `R`.


def riordan_matrix_by_recurrence(dim, rec, init={(0,0):1},
                                 ctor=zeros, post=expand,
                                 lattice=None):

    if not lattice: lattice = [(n, k) for n in range(1, dim)        # `lattice` denotes the order in which coeffs in
                                      for k in range(n+1)]          # the array will be computed; it can be plugged in.

    R = ctor(dim)                                                   # Initial array as base for recursive construction.

    for cell, i in init.items(): R[cell] = i                        # Set boundary values for the recurrence `rec`.

    for cell in lattice:                                            # Visit cells as they appear in the order `lattice`;
        for comb_cell, v in rec(*cell).items():                     # then, get dependencies with respect the each cell.

            try:
                comb = 1 if cell == comb_cell else R[comb_cell]     # If it is possible to access the dependee cell,
                combined = v * comb                                 # then perform the combination using the given `v`.
            except IndexError:
                combined = 0                                        # Otherwise, take `0` because combination fails.

            R[cell] += combined                                     # Finally, accumulate the current combination.

    if callable(post): R = R.applyfunc(post)                        # If some post processing is desired, do it.

    return lambda n, k: R[n, k]                                     # Return a lambda that uses the computed `R`.

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

    P = production_matrix(M, exp=False)                             # Ordinary production matrix.
    C = production_matrix(M, exp=True)                              # "Exponential" production matrix.

    is_ord = all(P[:1-i, 1] == P[i-1:, i]                           # A RA is ordinary if columns of its PM
                 for i in range(2, P.cols))                         # are shifted but equals.

    diagonals = {d: [C[j+d,j] for j in range(1,C.rows-d)]           # Fetch matrix diagonals for the exponential case.
                 for d in range(C.rows-3) }

    is_exp = all(map(is_arithmetic_progression,                     # A Riordan array is exponential if `diagonals`
                     diagonals.values()))                           # are arithmetic progressions, all of them.

    return nature(is_ord, is_exp)


def is_arithmetic_progression(prog):

    steps = len(prog)-1
    for _ in range(steps):                                          # Reduce the list by consecutive differences,
        prog = [(b-a).simplify() for a, b in zip(prog, prog[1:])]   # till a single object survives.

    assert len(prog) == 1                                           # Consistency check of the last comment.

    return prog.pop() == 0                                          # Finally, the reduced object has to vanish.


def compositional_inverse(h_eq, y=symbols('y'), check=True):

    spec, body = h_eq.lhs, h_eq.rhs                                 # Destructuring function $h$ denoted by `h_eq`,
    t, = spec.args                                                  # let `t` be the formal var of function $h$.
                                                                    # $\bar{h}(h(t)) = t \leftrightarrow \left[ \bar{h}(y) = t\, | \, y = h(t) \right]$ allows us
    sols = solve(Eq(y, body), t)                                    # to solve $y = h(t)$ with respect to $t$ because
    for sol in sols:                                                # $h(t)$ is known, denoted by `body`. For each
        L = Lambda(y, sol)                                          # solution `sol`, which depends on $y$, we build
        if L(0) == 0:                                               # a callable object `L`. If it vanishes in $0$ and
                                                                    # $\bar{h}(h(t))=t$, then it is a compositional
            if check: assert L(body).simplify() == t                # inverse of $h$.

            h_bar = Function(r'\bar{{ {} }}'.format(str(spec.func)))# Prepare the name for function $\bar{h}$ and
            eq = Eq(h_bar(y), sol.factor())                         # build the corresponding equation that defines
            return eq                                               # $\bar{h}(y)$, compositional inverse of $h(t)$.

    raise ValueError                                                # If the above code fails to return, raise an error.

def group_inverse(d_eq, h_eq, post=identity, check=True):

    t, y = symbols('t y')                                           # Formal symbols for symbolic functions $f$ and $g$.
    g_fn, f_fn = Function('g'), Function('f')

    with lift_to_Lambda(d_eq, return_eq=False) as D:                # Promote equations `d` and `h` to become callables
      with lift_to_Lambda(h_eq, return_eq=True) as H:               # objects, returning equations as well.
        f_eq = compositional_inverse(H(t), y, check)                # Let function $f$ be the compositional inverse of
        with lift_to_Lambda(f_eq, return_eq=False) as F:            # function $h$, then promote it as callable and
            F_t = F(t)                                              # (i) evaluate it at $t$;
            g = post(1/D(F_t))                                      # (ii) build and refine $g$, the first function of
            f = post(F_t)                                           # the new Ra; (iii) refine $f$, the second function.

    g_eq = Eq(g_fn(t), g.simplify())                                # Build corresponding equation objects with
    f_eq = Eq(f_fn(t), f.simplify())                                # simplified expressions.

    if check:                                                       # If it is required to certify the computed result,
        with lift_to_Lambda(g_eq, return_eq=False) as G:            # then promote the just built equations `g` and `f`
          with lift_to_Lambda(f_eq, return_eq=False) as F:          # to perform group operation in order to check that
            H_rhs = H(t).rhs                                        # it yields the identity element, which has
            assert (D(t)*G(H_rhs)).simplify() == 1                  # (i) $1$ as first component and
            assert F(H_rhs).simplify() == t                         # (ii) $t$ as second component.

    return g_eq, f_eq

