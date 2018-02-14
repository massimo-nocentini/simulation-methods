
from functools import lru_cache
from itertools import product

from sympy import Symbol, Integer, poly, Eq, expand, zeros, symbols, Matrix, factorial, IndexedBase, solve

def convolution(sequences, variable=Symbol('t'), op=max):
    
    if not sequences: return Integer(0)
    
    t = variable
    l = op(map(lambda s: poly(s.rhs, t).degree(), sequences))+1
    
    def convolve(a, b):
        
        new = sum(sum(a.coeff(t, j)*b.coeff(t, i-j) for j in range(i+1)) * t**i for i in range(l))
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


def riordan_matrix_by_convolution(d, h, t):

    @lru_cache(maxsize=None)
    def column(j):
        return convolution([column(j-1), h], t) if j else d
        #return convolution([d] + [h]*j, t)

    return lambda i, j: column(j).rhs.coeff(t, i).expand()



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

def frobenius_matrix(seq):

    m = len(seq)
    def F(n, k):
        return 1 if n+1 == k else -seq[-1-k] if n == m-1 else 0

    return F

def unit_vector(i, offset=-1):

    def U(n, k):
        return 1 if n == i+offset else 0
     
    return U

def Asequence_(M):

    A = zeros(M.rows-1, M.cols-1)
    t = symbols('t')
    funz1 = M[0,0]

    for i in range(A.rows):
        current_row = M[i+1,:]
        funz2 = sum(current_row[j]*t**(i-j) for j in range(i+2))
        funz = (funz2/funz1).series(t, n=A.cols).removeO()
        for j in range(A.cols):
            A[i, j]=funz.coeff(t, n=j)
        funz1=funz2

    return A#[1:,:-1]

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
    return PM[:-1, :-1]

def rows_shift_matrix(by):
    return lambda i, j: 1 if i + by == j else 0

def diagonal_func_matrix(f):
    return lambda i, j: f(i) if i == j else 0

def is_ordinary_RA(M, show_witness=False):

    PM = production_matrix(M)

    is_ord = True
    witness = []
    for i in range(2, PM.cols):
        if not (PM[:1-i, 1] == PM[i-1:, i]):
            is_ord = False
            witness.append(PM[:, i]) # for the witness we return the entire column

    return (is_ord, witness) if show_witness else is_ord


def is_exponential_RA(M, show_witness=False):

    C = production_matrix(M, exp=True) 

    diagonals = { d: [C[j+d,j] for j in range(1,C.rows-d)] 
                  for d in range(C.rows-2) }

    k = IndexedBase('k')
    sols = {}
    for d, l in diagonals.items():
        eqs = []
        unknowns = []
        for i in range(len(l)-1):
            a, b = l[i], l[i+1]
            k_d = k[d]
            eq = Eq(k_d, a-b)
            unknowns.append(k_d)
            eqs.append(eq)
        sols[d] = solve(eqs, unknowns)

    is_exp = all(len(offsets) == 1 for d, offsets in sols.items())

    return (is_exp, diagonals, sols) if show_witness else is_exp

