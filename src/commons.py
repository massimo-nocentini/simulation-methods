from contextlib import contextmanager, redirect_stdout
from sympy import Eq, Lambda, Function, Indexed, latex

def define(let, be, **kwds):

    if 'evaluate' not in kwds:                                      # If `evaluate` is already given, use it as it is,
        kwds['evaluate'] = False                                    # otherwise set to `False` to preevent evaluation 
                                                                    # by `Eq`, which implicitly do simplifications;
    return Eq(let, be, **kwds)                                      # finally, return an equation object.

@contextmanager
def lift_to_Lambda(eq, return_eq=False, lhs_handler=lambda args: []):
    '''
    I lift an equational function def into a callable object.

    >>> a = IndexedBase('a')
    >>> aeq = Eq(a[n], n+a[n-1])
    >>> with lift_to_Lambda(aeq, return_eq=True) as aEQ:
    ...     arec = aEQ(n+1)
    ... arec
    Eq(a[n + 1], n + a[n] + 1)

    >>> b = Function('b')
    >>> beq = Eq(b(n), n+b(n-1))
    >>> with lift_to_Lambda(beq, return_eq=True) as bEQ:
    ...     brec = bEQ(n+1)
    ... brec
    Eq(b(n + 1), n + b(n) + 1))
                
    '''
    lhs = eq.lhs
    args = (lhs.args[1:] if isinstance(lhs, Indexed) else 
            lhs.args if isinstance(lhs, Function) else 
            lhs_handler(lhs))
    yield Lambda(args, eq if return_eq else eq.rhs)

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
