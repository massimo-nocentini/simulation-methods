
from contextlib import contextmanager, redirect_stdout
from sympy import Eq, Lambda, Function, Indexed, latex, Subs, factorial, Symbol, binomial, Integer
from functools import partial

def define(let, be, ctor=Eq, **kwds,):

    if 'evaluate' not in kwds:       # If `evaluate` is already given, use it # as it is,
        kwds['evaluate'] = False     # otherwise set to `False` to preevent evaluation
                                     # by `Eq`, which implicitly do simplifications;
    return ctor(let, be, **kwds)     # finally, return an equation object.

@contextmanager
def lift_to_Lambda(eq, return_eq=False, lhs_handler=lambda args: []):

    lhs = eq.lhs
    args = (lhs.args[1:] if isinstance(lhs, Indexed) else    # get arguments wrt the type of `lhs` object;
            lhs.args     if isinstance(lhs, Function) else   # here we handle both function and subscript
            lhs_handler(lhs))                                # notations. Finally, `Lambda` is the
    yield Lambda(args, eq if return_eq else eq.rhs)          # class of callable objects in SymPy.

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

def foldr(f, init, iterable):
    def F(l):
        try:
            car, *cdr = l
            return f(car, F(cdr))
        except ValueError:
            return init()
    return F(iterable)

def foldl1(f, lst, key=identity):

    first, *rest = lst                                              # Destructure `lst` into car and cdr.
    accumulated = key(first)                                        # Init the accumulator with the first element.
    while rest:                                                     # While there is element in the cdr,
        current, *rest = rest                                       # take the next element to combine
        accumulated = f(accumulated, key(current))                  # with the accumulated value via the function `f`.

    return accumulated                                              # Finally, return the accumulated value.


class FEq(Eq):

    def __call__(self, *args, **kwds):
        with lift_to_Lambda(self, **kwds) as feq:
            applied = feq(*args)
            if isinstance(applied, Eq):
                subs = Subs(applied, self.lhs, self.rhs)
                setattr(subs, 'substitution', partial(
                    lambda inner_self, *args, **kwds: self, subs))
                return subs
            else:
                return applied

    def swap(self):
        return self.__class__(self.rhs, self.lhs)

    def __iter__(self):
        yield self.lhs, self.rhs

    def as_substitution(self):
        return dict(self)

    def __mod__(self, arg):

        def S(term):
            return term.subs(self, simultaneous=True)

        if isinstance(arg, Eq):
            lhs = self.lhs if self.lhs == arg.lhs else S(arg.lhs)
            rhs = S(arg.rhs)
            return arg.__class__(lhs, rhs, evaluate=False)
        else:
            return S(arg)

    def series(self, *args, **kwds):
        expr = self.rhs
        n = kwds.get('n', 6)
        kwds['n'] = n # simple update to make sure that SymPy's `series` uses the same `n`.
        kernel = kwds.pop('kernel', 'ordinary')
        t, *_ = args

        if kernel == 'exponential':
            handler = factorial 
        elif kernel == 'catalan':
            handler = lambda n: (n+1) / binomial(2*n, n)
        elif kernel == 'ordinary':
            handler = lambda _: Integer(1)
        else:
            raise ValueError

        s = expr.series(*args, **kwds)
        bigO = s.getO()
        s = sum([s.coeff(t,i) * handler(i) * t**i
                 for i in range(n)]) + (bigO if bigO else Integer(0))

        return define(self.lhs, define(self.rhs, s, ctor=FEq), ctor=FEq)

def commute_symbols(eq, commutative=True):

    symbols = eq.free_symbols
    commuted_symbols = {fs: Symbol(fs.name, commutative=commutative)
                        for fs in symbols
                        if not fs.is_commutative}

    def subs(eq):
        return define(eq.lhs.subs({c:o for o,c in commuted_symbols.items()}, 
                                  simultaneous=True),
                      eq.rhs.subs({c:o for o,c in commuted_symbols.items()}, 
                                  simultaneous=True))

    return define(eq.lhs.subs(commuted_symbols, simultaneous=True),
                  eq.rhs.subs(commuted_symbols, simultaneous=True)), subs



