
from functools import lru_cache

from sympy import Symbol, Integer, poly, Eq, expand, zeros

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
