
from functools import lru_cache, reduce
import itertools
import operator, math

from sympy import *

from commons import * # our own module with basic utilities

w = symbols('w0 w1', commutative=True) # to denote False and True, respectively.
u = symbols('u0:51', commutative=True)
o = symbols('o0:51', commutative=False)
emptybox = Symbol(r'␣', commutative=True)

class ty:
    
    def __init__(self, *types):
        try:
            iterable, = types
            self.types = list(iterable)
        except (ValueError, TypeError):
            self.types = types # the list of types that I depend on.
    
    def tyvars(self):
        vs = map(operator.methodcaller('tyvars'), self.types)
        return foldr(operator.or_, lambda: set(), vs)
    
    def label(self):
        raise NotImplemented # I'm an abstract type, nameless indeed.
    
    def gf_lhs(self):
        L = Function(self.label())
        return L(*self.tyvars())
    
    def gf(self):
        return [define(self.gf_lhs(), rhs, ctor=FEq) 
                for rhs in self.gf_rhs(*self.types)]
        
    def gf_rhs(self, *types):
        return self.definition(*types).gf_rhs(*types)
    
    def definition(self, *types):
        raise NotImplemented
        
    def gfs_space(self):
        return itertools.product(*map(lambda ty: ty.gf(), self.types))
    
    def __or__(self, other):
        return du(self, other)
    
    def __mul__(self, other):
        return cp(self, other)
    
    def __rpow__(self, base):
        
        if base == 2: 
            return powerset(self)
        elif base == -2: 
            return ipowerset(self) * lst(self)
        else:
            raise ValueError
    
    def __invert__(self):
        return cycle(self)

    def __getitem__(self, key):
        if not isinstance(key,  Basic):
            raise TypeError
        return self * tyvar(key)

class cp(ty):
        
    def gf_rhs(self, *types):
        return [foldr(lambda gf, acc: gf.rhs * acc, 
                      lambda: Integer(1), 
                      gfs)
                for gfs in self.gfs_space()]
        
    def label(self):
        return r'\times'

class du(ty):
    
    def label(self):
        return r'\cup'
            
    def gf_rhs(self, *types):
        return [foldr(lambda gf, acc: gf.rhs + acc,
                      lambda: Integer(0), 
                      gfs) 
                for gfs in self.gfs_space()]

class tyvar(ty):
    
    def label(self):
        return r'\mathcal{V}'
        
    def gf_rhs(self, sym):
        return [sym]
    
    def tyvars(self):
        sym, = self.types
        args = sym.args
        syms = filter(lambda a: a.is_symbol, args) if args else [sym]
        return set(syms)

class maybe(ty):
    
    def definition(self, alpha):
        return tyvar(emptybox) | alpha
    
    def label(self):
        return r'\mathcal{M}'


class rec(ty):
    
    def me(self):
        return tyvar(self.gf_lhs())
    
    def gf(self):
        eqs = super().gf()
        return [define(eq.lhs, sol, ctor=FEq) 
                for eq in eqs
                for sol in solve(define(eq.lhs, eq.rhs), [eq.lhs])]
    
class lst(rec):
         
    def definition(self, alpha):
        return cp() | (alpha * self.me())
    
    def label(self):
        return r'\mathcal{L}'

class nnlst(rec):
         
    def definition(self, alpha):
        return alpha | (alpha * self.me())
    
    def label(self):
        return r'\mathcal{L}_{+}'

class bin_tree(rec):
    
    def definition(self, alpha):
        return cp() | (alpha * self.me() * self.me())
    
    def label(self):
        return r'\mathcal{B}'

class nnbin_tree(rec):
    
    def definition(self, alpha):
        return alpha | (alpha * self.me() * self.me())
    
    def label(self):
        return r'\mathcal{B}_{+}'

def occupancy(eq, syms, objects='unlike', boxes='unlike'):

    bullet = Symbol(r'\bullet', commutative=False)
    circ = Symbol(r'\circ', commutative=True)
    
    osyms = [Symbol(s.name, commutative=False) for s in syms]

    def S(expr, assocs):
        return expr.subs(assocs, simultaneous=True)

    oemptybox = Symbol(emptybox.name, commutative=False)
    rhs = S(eq.rhs, {emptybox: oemptybox})

    if (objects, boxes) == ('unlike', 'unlike'):
        gf = S(rhs, dict(zip(syms, osyms)))
    elif (objects, boxes) == ('like', 'unlike'):
        gf = S(rhs, dict(zip(syms, itertools.repeat(bullet))))
    elif (objects, boxes) == ('unlike', 'like'):
        gf = rhs #S(rhs, {oemptybox: Integer(1)})
    elif (objects, boxes) == ('like', 'like'):
        gf = S(S(rhs, dict(zip(syms, itertools.repeat(circ)))), 
                 {oemptybox: Integer(1)})
    else:
        raise ValueError('Unknown configuration')

    f = Function('gf')
    return define(f(*gf.free_symbols), gf, ctor=FEq)

# ______________________________________________________________________________
# Basic concrete types
# ______________________________________________________________________________

truth, falsehood = tyvar(w[1]), tyvar(w[0])
boolean = truth | falsehood

