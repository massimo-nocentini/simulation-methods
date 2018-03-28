from contextlib import contextmanager, redirect_stdout

from sympy import Eq, Lambda, Function, Indexed, latex

def define(let, be, **kwds):
    return Eq(let, be, **kwds)

@contextmanager
def lift_to_Lambda(eq, return_eq=False, lhs_handler=lambda args: []):
    lhs = eq.lhs
    args = (lhs.args[1:] if isinstance(lhs, Indexed) else 
            lhs.args if isinstance(lhs, Function) else 
            lhs_handler(lhs))
    yield Lambda(args, eq if return_eq else eq.rhs)

def save_latex_repr(term, filename):
    with open(filename, 'w') as f:
        with redirect_stdout(f):
            print('.. math::\n\n\t{}'.format(latex(term)))

