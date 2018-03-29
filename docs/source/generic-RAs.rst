
Generic Riordan Arrays
######################

Component polynomials
=====================

.. testcode::
 
    from sympy import *
    from sympy.abc import n, i, N, x, lamda, phi, z, j, r, k, a, t, alpha
    from sympy.functions.elementary.integers import ceiling, floor

    from commons import *
    from sequences import *
    from matrix_functions import *

    R_cal, d = IndexedBase(r'\mathcal{R}'), IndexedBase('d') # helpers bases

.. doctest::

    >>> m = 4
    >>> R = define(R_cal[m], Matrix(m, m, riordan_matrix_by_recurrence(m, lambda n, k: {(n, k):d[n, k]}, init={(0,0):d[0,0]})))
    >>> R
    Eq(\mathcal{R}[4], Matrix([[d[0, 0], 0, 0, 0], [d[1, 0], d[1, 1], 0, 0], [d[2, 0], d[2, 1], d[2, 2], 0], [d[3, 0], d[3, 1], d[3, 2], d[3, 3]]]))

    >>> eigendata = spectrum(R)
    >>> eigendata # doctest: +NORMALIZE_WHITESPACE
    Eq(\sigma(\mathcal{R}[4]), 
       ({1: (\lambda[1], m[1]), 2: (\lambda[2], m[2]), 3: (\lambda[3], m[3]), 4: (\lambda[4], m[4])}, 
        {\lambda[1]: d[0, 0], \lambda[2]: d[1, 1], \lambda[3]: d[2, 2], \lambda[4]: d[3, 3]}, 
        {m[1]: 1, m[2]: 1, m[3]: 1, m[4]: 1}))

    >>> data, eigenvals, multiplicities = eigendata.rhs

.. testcode::
    :hide:

    save_latex_repr(eigendata, './source/latex-snippets/generic-RAs-0.rst')

.. include:: latex-snippets/generic-RAs-0.rst

.. doctest::

    >>> Phi_polynomials = component_polynomials(eigendata.rhs)
    >>> Phi_polynomials
    {(1, 1): Eq(\Phi_{ 1, 1 }(z), z**3/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - z**2*(\lambda[2] + \lambda[3] + \lambda[4])/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + z*(\lambda[2]*\lambda[3] + \lambda[2]*\lambda[4] + \lambda[3]*\lambda[4])/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - \lambda[2]*\lambda[3]*\lambda[4]/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4])), (2, 1): Eq(\Phi_{ 2, 1 }(z), -z**3/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + z**2*(\lambda[1] + \lambda[3] + \lambda[4])/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - z*(\lambda[1]*\lambda[3] + \lambda[1]*\lambda[4] + \lambda[3]*\lambda[4])/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + \lambda[1]*\lambda[3]*\lambda[4]/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4])), (3, 1): Eq(\Phi_{ 3, 1 }(z), z**3/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) - z**2*(\lambda[1] + \lambda[2] + \lambda[4])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) + z*(\lambda[1]*\lambda[2] + \lambda[1]*\lambda[4] + \lambda[2]*\lambda[4])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) - \lambda[1]*\lambda[2]*\lambda[4]/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4])), (4, 1): Eq(\Phi_{ 4, 1 }(z), -z**3/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) + z**2*(\lambda[1] + \lambda[2] + \lambda[3])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) - z*(\lambda[1]*\lambda[2] + \lambda[1]*\lambda[3] + \lambda[2]*\lambda[3])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) + \lambda[1]*\lambda[2]*\lambda[3]/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3))}

.. testcode::
    :hide:

    obj = [eq.factor() for k, eq in Phi_polynomials.items()] # pretty print
    save_latex_repr(obj, './source/latex-snippets/generic-RAs-1.rst', iterable=True)

.. include:: latex-snippets/generic-RAs-1.rst

.. doctest::

    >>> m = 8
    >>> R = define(R_cal[m], Matrix(m, m, riordan_matrix_by_recurrence(m, lambda n, k: {(n, k):1 if n == k else d[n, k]})))
    >>> R # doctest: +NORMALIZE_WHITESPACE
    Eq(\mathcal{R}[8], 
       Matrix([[1, 0, 0, 0, 0, 0, 0, 0], 
               [d[1, 0], 1, 0, 0, 0, 0, 0, 0], 
               [d[2, 0], d[2, 1], 1, 0, 0, 0, 0, 0], 
               [d[3, 0], d[3, 1], d[3, 2], 1, 0, 0, 0, 0], 
               [d[4, 0], d[4, 1], d[4, 2], d[4, 3], 1, 0, 0, 0], 
               [d[5, 0], d[5, 1], d[5, 2], d[5, 3], d[5, 4], 1, 0, 0], 
               [d[6, 0], d[6, 1], d[6, 2], d[6, 3], d[6, 4], d[6, 5], 1, 0], 
               [d[7, 0], d[7, 1], d[7, 2], d[7, 3], d[7, 4], d[7, 5], d[7, 6], 1]]))

    >>> eigendata = spectrum(R)
    >>> eigendata
    Eq(\sigma(\mathcal{R}[8]), ({1: (\lambda[1], m[1])}, {\lambda[1]: 1}, {m[1]: 8}))

    >>> data, eigenvals, multiplicities = eigendata.rhs

.. testcode::
    :hide:

    save_latex_repr(eigendata, './source/latex-snippets/generic-RAs-2.rst')

.. include:: latex-snippets/generic-RAs-2.rst

.. doctest::

    >>> Phi_polynomials = component_polynomials(eigendata.rhs)
    >>> Phi_polynomials
    {(1, 1): Eq(\Phi_{ 1, 1 }(z), 1), (1, 2): Eq(\Phi_{ 1, 2 }(z), z - \lambda[1]), (1, 3): Eq(\Phi_{ 1, 3 }(z), z**2/2 - z*\lambda[1] + \lambda[1]**2/2), (1, 4): Eq(\Phi_{ 1, 4 }(z), z**3/6 - z**2*\lambda[1]/2 + z*\lambda[1]**2/2 - \lambda[1]**3/6), (1, 5): Eq(\Phi_{ 1, 5 }(z), z**4/24 - z**3*\lambda[1]/6 + z**2*\lambda[1]**2/4 - z*\lambda[1]**3/6 + \lambda[1]**4/24), (1, 6): Eq(\Phi_{ 1, 6 }(z), z**5/120 - z**4*\lambda[1]/24 + z**3*\lambda[1]**2/12 - z**2*\lambda[1]**3/12 + z*\lambda[1]**4/24 - \lambda[1]**5/120), (1, 7): Eq(\Phi_{ 1, 7 }(z), z**6/720 - z**5*\lambda[1]/120 + z**4*\lambda[1]**2/48 - z**3*\lambda[1]**3/36 + z**2*\lambda[1]**4/48 - z*\lambda[1]**5/120 + \lambda[1]**6/720), (1, 8): Eq(\Phi_{ 1, 8 }(z), z**7/5040 - z**6*\lambda[1]/720 + z**5*\lambda[1]**2/240 - z**4*\lambda[1]**3/144 + z**3*\lambda[1]**4/144 - z**2*\lambda[1]**5/240 + z*\lambda[1]**6/720 - \lambda[1]**7/5040)}

.. testcode::
    :hide:

    obj = [Phi_polynomials[1,i] for i in range(1, m+1)] # pretty print
    save_latex_repr(obj, './source/latex-snippets/generic-RAs-3.rst', iterable=True)

.. include:: latex-snippets/generic-RAs-3.rst

which can be factored as follows

.. testcode::
    :hide:

    obj = [Phi_polynomials[1,i].factor() for i in range(1, m+1)] # pretty print
    save_latex_repr(obj, './source/latex-snippets/generic-RAs-4.rst', iterable=True)

.. include:: latex-snippets/generic-RAs-4.rst

.. doctest::

    >>> generic_coeff = define(d[n,k], (-lamda_indexed[1])**(n-k)/(factorial(n-k)*factorial(k)))
    >>> generic_coeff
    Eq(d[n, k], (-\lambda[1])**(-k + n)/(factorial(k)*factorial(-k + n)))

    >>> with lift_to_Lambda(generic_coeff) as G:
    ...     res_ordinary = M_ordinary, z_ordinary, Phi_ordinary = (
    ...         Matrix(m, m, lambda n,k: G(n,k) if k <= n else 0),
    ...         Matrix([z**i for i in range(m)]),
    ...         Matrix([Function(r'\Phi_{{ {}, {} }}'.format(1, j)).__call__(z) for j in range(1, m+1)]))
    >>> res_ordinary
    (Matrix([
    [                  1,                  0,                  0,                 0,                  0,                 0,               0,      0],
    [        -\lambda[1],                  1,                  0,                 0,                  0,                 0,               0,      0],
    [    \lambda[1]**2/2,        -\lambda[1],                1/2,                 0,                  0,                 0,               0,      0],
    [   -\lambda[1]**3/6,    \lambda[1]**2/2,      -\lambda[1]/2,               1/6,                  0,                 0,               0,      0],
    [   \lambda[1]**4/24,   -\lambda[1]**3/6,    \lambda[1]**2/4,     -\lambda[1]/6,               1/24,                 0,               0,      0],
    [ -\lambda[1]**5/120,   \lambda[1]**4/24,  -\lambda[1]**3/12,  \lambda[1]**2/12,     -\lambda[1]/24,             1/120,               0,      0],
    [  \lambda[1]**6/720, -\lambda[1]**5/120,   \lambda[1]**4/48, -\lambda[1]**3/36,   \lambda[1]**2/48,   -\lambda[1]/120,           1/720,      0],
    [-\lambda[1]**7/5040,  \lambda[1]**6/720, -\lambda[1]**5/240, \lambda[1]**4/144, -\lambda[1]**3/144, \lambda[1]**2/240, -\lambda[1]/720, 1/5040]]), Matrix([
    [   1],
    [   z],
    [z**2],
    [z**3],
    [z**4],
    [z**5],
    [z**6],
    [z**7]]), Matrix([
    [\Phi_{ 1, 1 }(z)],
    [\Phi_{ 1, 2 }(z)],
    [\Phi_{ 1, 3 }(z)],
    [\Phi_{ 1, 4 }(z)],
    [\Phi_{ 1, 5 }(z)],
    [\Phi_{ 1, 6 }(z)],
    [\Phi_{ 1, 7 }(z)],
    [\Phi_{ 1, 8 }(z)]]))
    
.. testcode::
    :hide:

    save_latex_repr(res_ordinary, './source/latex-snippets/generic-RAs-5.rst')

.. include:: latex-snippets/generic-RAs-5.rst
    
Above matrix is `A098361 <http://oeis.org/A098361>`_, triangular shaped.

.. doctest::

    >>> expt_coeff = Eq(generic_coeff.lhs, generic_coeff.rhs*factorial(k))
    >>> expt_coeff
    Eq(d[n, k], (-\lambda[1])**(-k + n)/factorial(-k + n))

    >>> with lift_to_Lambda(expt_coeff) as G:
    ...     res_expt = (Matrix(m, m, lambda n,k: G(n,k) if k <= n else 0),
    ...                 Matrix([z**i / factorial(i, evaluate=i<2) for i in range(m)]),
    ...                 Matrix([Function(r'\Phi_{{ {}, {} }}'.format(1, j)).__call__(z) for j in range(1, m+1)]))
    >>> res_expt
    (Matrix([         
    [                  1,                  0,                  0,                0,                0,               0,           0, 0],
    [        -\lambda[1],                  1,                  0,                0,                0,               0,           0, 0],
    [    \lambda[1]**2/2,        -\lambda[1],                  1,                0,                0,               0,           0, 0],
    [   -\lambda[1]**3/6,    \lambda[1]**2/2,        -\lambda[1],                1,                0,               0,           0, 0],
    [   \lambda[1]**4/24,   -\lambda[1]**3/6,    \lambda[1]**2/2,      -\lambda[1],                1,               0,           0, 0],
    [ -\lambda[1]**5/120,   \lambda[1]**4/24,   -\lambda[1]**3/6,  \lambda[1]**2/2,      -\lambda[1],               1,           0, 0],
    [  \lambda[1]**6/720, -\lambda[1]**5/120,   \lambda[1]**4/24, -\lambda[1]**3/6,  \lambda[1]**2/2,     -\lambda[1],           1, 0],
    [-\lambda[1]**7/5040,  \lambda[1]**6/720, -\lambda[1]**5/120, \lambda[1]**4/24, -\lambda[1]**3/6, \lambda[1]**2/2, -\lambda[1], 1]]), Matrix([
    [                1],
    [                z],
    [z**2/factorial(2)],
    [z**3/factorial(3)],
    [z**4/factorial(4)],
    [z**5/factorial(5)],
    [z**6/factorial(6)], 
    [z**7/factorial(7)]]), Matrix([
    [\Phi_{ 1, 1 }(z)],
    [\Phi_{ 1, 2 }(z)],
    [\Phi_{ 1, 3 }(z)],
    [\Phi_{ 1, 4 }(z)],
    [\Phi_{ 1, 5 }(z)],
    [\Phi_{ 1, 6 }(z)],
    [\Phi_{ 1, 7 }(z)],
    [\Phi_{ 1, 8 }(z)]]))

    >>> M_expt, z_expt, Phi_expt = res_expt

    
.. testcode::
    :hide:

    save_latex_repr(res_expt, './source/latex-snippets/generic-RAs-6.rst')

.. include:: latex-snippets/generic-RAs-6.rst
    

.. doctest::

    >>> # we need to manipulate a bit because some terms are factored and some are left unevaluated
    >>> assert (M_ordinary*z_ordinary.applyfunc(lambda i: i.collect(z)) == 
    ...         M_expt*z_expt.applyfunc(lambda i: i.doit()))

Moreover, it is interesting the shape of :code:`M_expt**(-1)` which has
positive coefficients only:

.. doctest::

    >>> M_expt**(-1)
    Matrix([
    [                 1,                 0,                 0,                0,               0,               0,          0, 0],
    [        \lambda[1],                 1,                 0,                0,               0,               0,          0, 0],
    [   \lambda[1]**2/2,        \lambda[1],                 1,                0,               0,               0,          0, 0],
    [   \lambda[1]**3/6,   \lambda[1]**2/2,        \lambda[1],                1,               0,               0,          0, 0],
    [  \lambda[1]**4/24,   \lambda[1]**3/6,   \lambda[1]**2/2,       \lambda[1],               1,               0,          0, 0],
    [ \lambda[1]**5/120,  \lambda[1]**4/24,   \lambda[1]**3/6,  \lambda[1]**2/2,      \lambda[1],               1,          0, 0],
    [ \lambda[1]**6/720, \lambda[1]**5/120,  \lambda[1]**4/24,  \lambda[1]**3/6, \lambda[1]**2/2,      \lambda[1],          1, 0],
    [\lambda[1]**7/5040, \lambda[1]**6/720, \lambda[1]**5/120, \lambda[1]**4/24, \lambda[1]**3/6, \lambda[1]**2/2, \lambda[1], 1]])
    
.. testcode::
    :hide:

    save_latex_repr(_, './source/latex-snippets/generic-RAs-7.rst')

.. include:: latex-snippets/generic-RAs-7.rst

Functions
=========

.. doctest::

    >>> f, g = Function('f'), Function('g') # abstract function and its Hermite interpolating polynomial

:math:`f(z)=z^{r}, r\in\mathbb{Q}` 
----------------------------------

.. doctest::

    >>> f_power = define(f(z), z**r)
    >>> f_power
    Eq(f(z), z**r)

    >>> g_power = g_poly(f_power, eigendata, Phi_polynomials)
    >>> g_power
    Eq(g(z), -r**7*\lambda[1]**r/5040 + r**6*\lambda[1]**r/180 - 23*r**5*\lambda[1]**r/360 + 7*r**4*\lambda[1]**r/18 - 967*r**3*\lambda[1]**r/720 + 469*r**2*\lambda[1]**r/180 - 363*r*\lambda[1]**r/140 + z**7*(r**7*\lambda[1]**r/(5040*\lambda[1]**7) - r**6*\lambda[1]**r/(240*\lambda[1]**7) + 5*r**5*\lambda[1]**r/(144*\lambda[1]**7) - 7*r**4*\lambda[1]**r/(48*\lambda[1]**7) + 29*r**3*\lambda[1]**r/(90*\lambda[1]**7) - 7*r**2*\lambda[1]**r/(20*\lambda[1]**7) + r*\lambda[1]**r/(7*\lambda[1]**7)) + z**6*(-r**7*\lambda[1]**r/(720*\lambda[1]**6) + 11*r**6*\lambda[1]**r/(360*\lambda[1]**6) - 19*r**5*\lambda[1]**r/(72*\lambda[1]**6) + 41*r**4*\lambda[1]**r/(36*\lambda[1]**6) - 1849*r**3*\lambda[1]**r/(720*\lambda[1]**6) + 1019*r**2*\lambda[1]**r/(360*\lambda[1]**6) - 7*r*\lambda[1]**r/(6*\lambda[1]**6)) + z**5*(r**7*\lambda[1]**r/(240*\lambda[1]**5) - 23*r**6*\lambda[1]**r/(240*\lambda[1]**5) + 69*r**5*\lambda[1]**r/(80*\lambda[1]**5) - 185*r**4*\lambda[1]**r/(48*\lambda[1]**5) + 134*r**3*\lambda[1]**r/(15*\lambda[1]**5) - 201*r**2*\lambda[1]**r/(20*\lambda[1]**5) + 21*r*\lambda[1]**r/(5*\lambda[1]**5)) + z**4*(-r**7*\lambda[1]**r/(144*\lambda[1]**4) + r**6*\lambda[1]**r/(6*\lambda[1]**4) - 113*r**5*\lambda[1]**r/(72*\lambda[1]**4) + 22*r**4*\lambda[1]**r/(3*\lambda[1]**4) - 2545*r**3*\lambda[1]**r/(144*\lambda[1]**4) + 41*r**2*\lambda[1]**r/(2*\lambda[1]**4) - 35*r*\lambda[1]**r/(4*\lambda[1]**4)) + z**3*(r**7*\lambda[1]**r/(144*\lambda[1]**3) - 25*r**6*\lambda[1]**r/(144*\lambda[1]**3) + 247*r**5*\lambda[1]**r/(144*\lambda[1]**3) - 1219*r**4*\lambda[1]**r/(144*\lambda[1]**3) + 389*r**3*\lambda[1]**r/(18*\lambda[1]**3) - 949*r**2*\lambda[1]**r/(36*\lambda[1]**3) + 35*r*\lambda[1]**r/(3*\lambda[1]**3)) + z**2*(-r**7*\lambda[1]**r/(240*\lambda[1]**2) + 13*r**6*\lambda[1]**r/(120*\lambda[1]**2) - 9*r**5*\lambda[1]**r/(8*\lambda[1]**2) + 71*r**4*\lambda[1]**r/(12*\lambda[1]**2) - 3929*r**3*\lambda[1]**r/(240*\lambda[1]**2) + 879*r**2*\lambda[1]**r/(40*\lambda[1]**2) - 21*r*\lambda[1]**r/(2*\lambda[1]**2)) + z*(r**7*\lambda[1]**r/(720*\lambda[1]) - 3*r**6*\lambda[1]**r/(80*\lambda[1]) + 59*r**5*\lambda[1]**r/(144*\lambda[1]) - 37*r**4*\lambda[1]**r/(16*\lambda[1]) + 319*r**3*\lambda[1]**r/(45*\lambda[1]) - 223*r**2*\lambda[1]**r/(20*\lambda[1]) + 7*r*\lambda[1]**r/\lambda[1]) + \lambda[1]**r)

.. testcode::
    :hide:

    save_latex_repr(g_power, './source/latex-snippets/generic-RAs-8.rst')

.. include:: latex-snippets/generic-RAs-8.rst

.. doctest::

    >>> g_power = g_power.subs(eigenvals)
    >>> g_power
    Eq(g(z), -r**7/5040 + r**6/180 - 23*r**5/360 + 7*r**4/18 - 967*r**3/720 + 469*r**2/180 - 363*r/140 + z**7*(r**7/5040 - r**6/240 + 5*r**5/144 - 7*r**4/48 + 29*r**3/90 - 7*r**2/20 + r/7) + z**6*(-r**7/720 + 11*r**6/360 - 19*r**5/72 + 41*r**4/36 - 1849*r**3/720 + 1019*r**2/360 - 7*r/6) + z**5*(r**7/240 - 23*r**6/240 + 69*r**5/80 - 185*r**4/48 + 134*r**3/15 - 201*r**2/20 + 21*r/5) + z**4*(-r**7/144 + r**6/6 - 113*r**5/72 + 22*r**4/3 - 2545*r**3/144 + 41*r**2/2 - 35*r/4) + z**3*(r**7/144 - 25*r**6/144 + 247*r**5/144 - 1219*r**4/144 + 389*r**3/18 - 949*r**2/36 + 35*r/3) + z**2*(-r**7/240 + 13*r**6/120 - 9*r**5/8 + 71*r**4/12 - 3929*r**3/240 + 879*r**2/40 - 21*r/2) + z*(r**7/720 - 3*r**6/80 + 59*r**5/144 - 37*r**4/16 + 319*r**3/45 - 223*r**2/20 + 7*r) + 1)

.. testcode::
    :hide:

    save_latex_repr(g_power, './source/latex-snippets/generic-RAs-9.rst')

.. include:: latex-snippets/generic-RAs-9.rst

.. doctest::

    >>> g_power_theo_jk = define(g(z), sum((-1)**(j-1)*binomial(r, j-1)*binomial(j-1, k)*(-z)**k
    ...                                    for j in range(1, m+1)
    ...                                    for k in range(j)).collect(z))
    >>> g_power_theo_kj = define(g(z), sum((-1)**(j)*binomial(r, j)*binomial(j, k)*(-z)**k 
    ...                                    for k in range(m) 
    ...                                    for j in range(k, m)).collect(z))
    >>> assert (g_power.rhs.collect(z) == 
    ...         g_power_theo_jk.rhs.combsimp().collect(z) == 
    ...         g_power_theo_kj.rhs.combsimp().collect(z))

.. testcode::
    :hide:

    save_latex_repr(g_power_theo_kj, './source/latex-snippets/generic-RAs-10.rst')

.. include:: latex-snippets/generic-RAs-10.rst



