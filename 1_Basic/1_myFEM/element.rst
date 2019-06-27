My Element
============

We want to start with a simple linear triangular finite element for the space :math:`H^1`.
First we create a base class, so that we can later create higher order elements as well. This base class requires elements to be able to compute values and gradients:

.. literalinclude:: myElement.hpp
   :language: cpp
   :start-at: class MyBaseElement
   :end-at: };

The linear trig implements this functionality:

.. literalinclude:: myElement.hpp
   :language: cpp
   :start-at: class MyLinearTrig
   :end-at: };

So we need to define 3 functions: The constructor, ``CalcShape`` and ``CalcDShape`` to compute
the shape functions and their derivatives.

.. literalinclude:: myElement.cpp
   :language: cpp
   :start-at: void MyLinearTrig :: CalcShape
   :end-before: void MyQuadraticTrig :: CalcShape

Automatic Differentiation
-------------------------

Next we will create second order elements for our space. When creating higher order spaces the derivatives of the polynomials will get quite ugly, so we don't want to implement them by hand. NGSolve provides automatic exact differentiation by the `AutoDiff` class. A `AutoDiff` object knows its value and the value of its derivative. Combining `AutoDiff` objects automatically uses differentiation rules.

We implement nodal quadratic elements, note `CalcShape` and `CalcDShape` look quite similar, in `CalcDShape` we only use the diff value of the `Autodiff` variable:

.. literalinclude:: myElement.cpp
   :language: cpp
   :start-at: void MyQuadraticTrig :: CalcShape

The creation of boundary elements will be left as an exercise.
Next we will create (differential) operators working with these finite elements.
