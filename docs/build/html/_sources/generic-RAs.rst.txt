
Generic Riordan Arrays
======================

.. testcode::
 
    from sympy import *
    from sympy.abc import n, i, N, x, lamda, phi, z, j, r, k, a, t, alpha
    from sympy.functions.elementary.integers import ceiling, floor

    from sequences import *
    from matrix_functions import *

.. doctest::

    >>> m = 4
    >>> d = IndexedBase('d')
    >>> R = Matrix(m, m, riordan_matrix_by_recurrence(m, lambda n, k: {(n, k):d[n, k]}, init={(0,0):d[0,0]}))
    >>> R
    Matrix([
    [d[0, 0],       0,       0,       0],
    [d[1, 0], d[1, 1],       0,       0],
    [d[2, 0], d[2, 1], d[2, 2],       0],
    [d[3, 0], d[3, 1], d[3, 2], d[3, 3]]])

    >>> data, eigenvals, multiplicities = eigendata = eigen_data(R)
    >>> eigendata # doctest: +NORMALIZE_WHITESPACE
    ({1: (\lambda[1], m[1]), 2: (\lambda[2], m[2]), 3: (\lambda[3], m[3]), 4: (\lambda[4], m[4])}, 
     {\lambda[1]: d[0, 0], \lambda[2]: d[1, 1], \lambda[3]: d[2, 2], \lambda[4]: d[3, 3]}, 
     {m[1]: 1, m[2]: 1, m[3]: 1, m[4]: 1})

.. testcode::
    :hide:

    save_latex_repr(eigendata, './source/latex-snippets/generic-RAs-0.tex')

.. include:: latex-snippets/generic-RAs-0.tex

.. doctest::

    >>> Phi_polynomials = component_polynomials(eigendata)
    >>> Phi_polynomials
    {(1, 1): Eq(\Phi_{ 1, 1 }(z), z**3/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - z**2*(\lambda[2] + \lambda[3] + \lambda[4])/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + z*(\lambda[2]*\lambda[3] + \lambda[2]*\lambda[4] + \lambda[3]*\lambda[4])/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - \lambda[2]*\lambda[3]*\lambda[4]/(\lambda[1]**3 - \lambda[1]**2*\lambda[2] - \lambda[1]**2*\lambda[3] - \lambda[1]**2*\lambda[4] + \lambda[1]*\lambda[2]*\lambda[3] + \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4])), (2, 1): Eq(\Phi_{ 2, 1 }(z), -z**3/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + z**2*(\lambda[1] + \lambda[3] + \lambda[4])/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) - z*(\lambda[1]*\lambda[3] + \lambda[1]*\lambda[4] + \lambda[3]*\lambda[4])/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4]) + \lambda[1]*\lambda[3]*\lambda[4]/(\lambda[1]*\lambda[2]**2 - \lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]**3 + \lambda[2]**2*\lambda[3] + \lambda[2]**2*\lambda[4] - \lambda[2]*\lambda[3]*\lambda[4])), (3, 1): Eq(\Phi_{ 3, 1 }(z), z**3/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) - z**2*(\lambda[1] + \lambda[2] + \lambda[4])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) + z*(\lambda[1]*\lambda[2] + \lambda[1]*\lambda[4] + \lambda[2]*\lambda[4])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4]) - \lambda[1]*\lambda[2]*\lambda[4]/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]**2 + \lambda[1]*\lambda[3]*\lambda[4] - \lambda[2]*\lambda[3]**2 + \lambda[2]*\lambda[3]*\lambda[4] + \lambda[3]**3 - \lambda[3]**2*\lambda[4])), (4, 1): Eq(\Phi_{ 4, 1 }(z), -z**3/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) + z**2*(\lambda[1] + \lambda[2] + \lambda[3])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) - z*(\lambda[1]*\lambda[2] + \lambda[1]*\lambda[3] + \lambda[2]*\lambda[3])/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3) + \lambda[1]*\lambda[2]*\lambda[3]/(\lambda[1]*\lambda[2]*\lambda[3] - \lambda[1]*\lambda[2]*\lambda[4] - \lambda[1]*\lambda[3]*\lambda[4] + \lambda[1]*\lambda[4]**2 - \lambda[2]*\lambda[3]*\lambda[4] + \lambda[2]*\lambda[4]**2 + \lambda[3]*\lambda[4]**2 - \lambda[4]**3))}

.. testcode::
    :hide:

    obj = [eq.factor() for k, eq in Phi_polynomials.items()] # pretty print
    save_latex_repr(obj, './source/latex-snippets/generic-RAs-1.tex')

.. include:: latex-snippets/generic-RAs-1.tex

