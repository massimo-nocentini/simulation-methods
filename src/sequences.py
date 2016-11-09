
from functools import lru_cache

from sympy import Symbol, Integer, poly, Eq

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



