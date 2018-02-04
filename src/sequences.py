
from functools import lru_cache

from sympy import Symbol, Integer, poly, Eq, expand, zeros, symbols, Matrix, factorial

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

def production_matrix(M, exp=False):
    U = Matrix(M.rows, M.cols, lambda i, j: 1 if i+1 == j else 0)
    pm = M**(-1) * U * M
    pm = pm[:-1, :-1]
    if exp:
        for n in range(pm.rows):
            for k in range(pm.cols):
                pm[n, k] = pm[n, k] / (factorial(n)/factorial(k))
    return pm



