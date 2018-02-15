
from contextlib import contextmanager

from sympy import Eq, Lambda

def define(let, be, **kwds):
    return Eq(let, be, **kwds)

@contextmanager
def lift_to_Lambda(eq, return_eq=False):
    yield Lambda(eq.lhs.args, eq if return_eq else eq.rhs)

