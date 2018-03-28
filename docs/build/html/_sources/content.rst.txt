
Title
=====

Some text

.. doctest::

    >>> from contextlib import redirect_stdout
    >>> from sympy import *

    >>> t = symbols('t')
    >>> term = Sum(t, (t, 0, oo))
    >>> term
    Sum(t, (t, 0, oo))
    
.. testcode::
    :hide:

    with open('./source/my-sum.rst', 'w') as f:
        with redirect_stdout(f):
            print('.. math::\n\n\t{}'.format(latex(term)))

Produces

.. include:: my-sum.rst

.. sidebar:: Sidebar Title
    :subtitle: Optional Sidebar Subtitle

    Subsequent indented lines comprise
    the body of the sidebar, and are
    interpreted as body elements.

.. doctest::

    >>> 3
    3

another test added three words.  
Lorem ipsum [#f1]_ dolor sit amet ... [#f2]_

.. rubric:: Footnotes

.. [#f1] Text of the first footnote.
.. [#f2] Text of the second footnote.
