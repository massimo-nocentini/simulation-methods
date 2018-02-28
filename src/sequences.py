
from functools import lru_cache
from itertools import product
from collections import namedtuple

from sympy import Symbol, Integer, poly, Eq, expand, zeros, symbols, Matrix, factorial, IndexedBase, solve, simplify

from commons import *

nature = namedtuple('nature', ['is_ordinary', 'is_exponential'])

def convolution(sequences, t):
    
    if not sequences: return Integer(0)
    
    l = max(map(lambda s: poly(s.rhs, t).degree(), sequences))+1
    
    def convolve(a, b):
        
        new = sum(sum(a.coeff(t, j)*b.coeff(t, i-j) for j in range(i+1)) * t**i 
                  for i in range(l))
        p = poly(new, t).args[0]
        return sum(p.coeff(t, i)*t**i for i in range(l))
    
    conv_eq, *rest = sequences
    conv_def, conv = conv_eq.lhs, conv_eq.rhs
    while rest:
        another_eq, *rest = rest
        b_def, b = another_eq.lhs, another_eq.rhs
        conv_def *= b_def
        conv = convolve(conv, b)
    
    return Eq(conv_def, conv)


def riordan_matrix_by_convolution(dim, d, h):

    t = symbols('t')

    with lift_to_Lambda(d, return_eq=True) as D, \
         lift_to_Lambda(h, return_eq=True) as H:
        d_eq, h_eq = D(t), H(t)
        d_eq = Eq(d_eq.lhs, d_eq.rhs.series(t, n=dim).removeO())
        h_eq = Eq(h_eq.lhs, h_eq.rhs.series(t, n=dim).removeO())

    @lru_cache(maxsize=None)
    def column(j):
        return convolution([column(j-1), h_eq], t) if j else d_eq
        #return convolution([d] + [h]*j, t)

    return lambda i, j: column(j).rhs.coeff(t, i).expand()

def riordan_matrix_by_AZ_sequences(dim, seqs, init={(0,0):1}, ctor=zeros, post=expand, lattice=None):
    
    if not lattice: lattice = [(n, k) for n in range(1, dim) for k in range(n+1)]

    Zseq, Aseq = seqs
    t = symbols('t')

    with lift_to_Lambda(Zseq) as Z, lift_to_Lambda(Aseq) as A:
        Z_series = Z(t).series(t, n=dim).removeO()
        A_series = A(t).series(t, n=dim).removeO()
    
    R = ctor(dim)
    
    for cell, i in init.items(): R[cell] = i
    
    for n, k in lattice:

        if k:
            v = sum(R[n-1, j] * A_series.coeff(t, n=i) for i, j in enumerate(range(k-1, dim)))
        else:
            v = sum(R[n-1, j] * Z_series.coeff(t, n=j) for j in range(dim))
        
        R[n, k] = v

    
    if callable(post): R = R.applyfunc(post)

    return lambda n, k: R[n, k]


def riordan_matrix_by_recurrence(dim, rec, init={(0,0):1}, ctor=zeros, post=expand, lattice=None):
    
    if not lattice: lattice = [(n, k) for n in range(1, dim) for k in range(n+1)]
    
    R = ctor(dim)
    
    for cell, i in init.items(): R[cell] = i
    
    for cell in lattice:
        for comb_cell, v in rec(*cell).items():
            
            try:
                combined = v * (1 if cell == comb_cell else R[comb_cell])
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

def production_matrix(M, exp=False):
    """
    Returns the production matrix of the given RA `M`.

    Implementation according to Barry's book ``RA: a Primer'', page 215.

    """

    U = Matrix(M.rows, M.cols, rows_shift_matrix(by=1))
    F = Matrix(M.rows, M.cols, diagonal_func_matrix(f=factorial if exp else one))
    F_inv = F**(-1)
    V = F_inv * U * F
    O = F_inv * M * F
    O_inv = O**(-1)
    PM = F * O_inv * V * O * F_inv
    PM = F_inv * PM * F if exp else PM
    PM = PM[:-1, :-1]
    return PM.applyfunc(simplify)

def rows_shift_matrix(by):
    return lambda i, j: 1 if i + by == j else 0

def diagonal_func_matrix(f):
    return lambda i, j: f(i) if i == j else 0

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

    return prog[0] == 0


def compositional_inverse(h_eq, y=symbols('y'), check=True):

    spec, body = h_eq.lhs, h_eq.rhs
    t, = spec.args

    sols = solve(Eq(y, body), t)
    for sol in sols:
        L = Lambda(y, sol)
        if L(0) == 0:

            # partialmethod(L(body), post) # TODO
            if check: assert L(body).simplify() == t

            h_bar = Function(r'\bar{{ {} }}'.format(str(spec.func)))
            eq = Eq(h_bar(y), sol.factor())
            return eq

    raise ValueError




