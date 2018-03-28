
Title
=====


.. testcode::
    :hide:

    from commons import save_latex_repr

Some text

.. doctest::

    >>> from sympy import *

    >>> t = symbols('t')
    >>> term = Sum(t, (t, 0, oo))
    >>> term
    Sum(t, (t, 0, oo))
    
.. testcode::
    :hide:

    save_latex_repr(term, './source/my-sum.rst')

Produces

.. include:: my-sum.rst


.. doctest::

    >>> 3
    3

another test added three words.  
Lorem ipsum [#f1]_ dolor sit amet ... [#f2]_

.. rubric:: Footnotes

.. [#f1] Text of the first footnote.
.. [#f2] Text of the second footnote.
