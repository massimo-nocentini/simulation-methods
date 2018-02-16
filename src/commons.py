
from contextlib import contextmanager

from sympy import Eq, Lambda, Function, Indexed

def define(let, be, **kwds):
    return Eq(let, be, **kwds)

@contextmanager
def lift_to_Lambda(eq, return_eq=False, lhs_handler=lambda args: []):
    lhs = eq.lhs
    args = (lhs.args[1:] if isinstance(lhs, Indexed) else 
            lhs.args if isinstance(lhs, Function) else 
            lhs_handler(lhs))
    yield Lambda(args, eq if return_eq else eq.rhs)

