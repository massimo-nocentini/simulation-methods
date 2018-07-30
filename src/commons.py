from contextlib import contextmanager, redirect_stdout
from sympy import Eq, Lambda, Function, Indexed, latex

def define(let, be, **kwds):

    if 'evaluate' not in kwds:                                      # If `evaluate` is already given, use it as it is,
        kwds['evaluate'] = False                                    # otherwise set to `False` to preevent evaluation 
                                                                    # by `Eq`, which implicitly do simplifications;
    return Eq(let, be, **kwds)                                      # finally, return an equation object.

@contextmanager
def lift_to_Lambda(eq, return_eq=False, lhs_handler=lambda args: []):

    lhs = eq.lhs
    args = (lhs.args[1:] if isinstance(lhs, Indexed) else           # get arguments with respect to the type of `lhs`
            lhs.args     if isinstance(lhs, Function) else          # object; here we handle both function and
            lhs_handler(lhs))                                       # subscript notations. Finally, `Lambda` is the
    yield Lambda(args, eq if return_eq else eq.rhs)                 # class of callable objects in SymPy.

def save_latex_repr(term, filename, iterable=False):

    with open(filename, 'w') as f:
        with redirect_stdout(f):
            print('.. math::\n\n', end='')
            if iterable:
                for subterm in term:
                    print('\t& ' + r'{}\\'.format(latex(subterm)))
            else:
                print('\t' + latex(term))

def identity(*args):
    return args if len(args) > 1 else args[0]                       # Return the given arguments, simply.

def foldl1(f, lst, key=identity):
    
    first, *rest = lst                                              # Destructure `lst` into car and cdr.
    accumulated = key(first)                                        # Init the accumulator with the first element.
    while rest:                                                     # While there is element in the cdr,
        current, *rest = rest                                       # take the next element to combine
        accumulated = f(accumulated, key(current))                  # with the accumulated value via the function `f`.

    return accumulated                                              # Finally, return the accumulated value.
