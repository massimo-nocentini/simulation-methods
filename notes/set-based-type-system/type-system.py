
from functools import lru_cache, reduce
import itertools
import operator, math

from sympy import *

from commons import * # our own module with basic utilities

class ty:
    
    def __init__(self, *types):
        self.types = types # the list of types that I depend on.
    
    def tyvars(self):
        vs = map(operator.methodcaller('tyvars'), self.types)
        return reduce(operator.or_, vs, set())
    
    def label(self):
        raise NotImplemented # I'm an abstract type, nameless indeed.
    
    def gf_lhs(self):
        L = Function(self.label())
        return L(*self.tyvars())
    
    def gf(self):
        return [define(self.gf_lhs(), rhs.simplify(), ctor=FEq) 
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

class cp(ty):
        
    def gf_rhs(self, *types):
        return [foldr(lambda gf, acc: gf.rhs * acc, 
                      lambda: Integer(1), gfs)
                for gfs in self.gfs_space()]
        
    def label(self):
        return r'\times'

class du(ty):
    
    def label(self):
        return r'\cup'
            
    def gf_rhs(self, *types):
        return [reduce(lambda acc, gf: gf.rhs + acc, gfs, Integer(0)) 
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
        return cp() | alpha
    
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
