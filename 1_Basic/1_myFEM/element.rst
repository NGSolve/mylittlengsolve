My Element
============

We want to start with a simple linear triangular finite element for the space :math:`H^1`.
Therefore we can use the ``ScalarFiniteElement`` as a base class:

.. literalinclude:: myElement.hpp
   :language: cpp
   :start-at: class MyLinearTrig
   :end-at: };

So we need to define 3 functions: The constructor, ``CalcShape`` and ``CalcDShape`` to compute
the shape functions and their derivatives.

.. literalinclude:: myElement.cpp
   :language: cpp
   :start-at: MyLinearTrig :: MyLinearTrig
   :end-before: MyQuadraticTrig :: MyQuad

Finally we can export it to Python, as we have done it with the CoefficientFunction. Note that this is
not necessary if we just want to build a finite element space with it.

.. literalinclude:: myElement.cpp
   :language: cpp
   :start-at: void Export

.. todo::

   Explain AutoDiff with second order elements

Next we want to build this finite element into a space, so we can use it.
