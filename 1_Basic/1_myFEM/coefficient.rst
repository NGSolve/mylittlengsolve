My CoefficientFunction
========================

This tutorial shows how to create arbitrary CoeffientFunctions. A can be done in Python but sometimes the functionality there is not enough. 

A new CF can be created pretty easily: Derive from the base class provide a constructor and override the pure virtual method `Evaluate`. The CF in our example should be able to read a value from a file and return that value. If the value in the file changes then the CF should change as well.

The simplest implementation of the CF would look like this:

.. literalinclude:: myCoefficient.hpp
   :start-after: #define __FILE_MYCOEFFICIENT_HPP
   :end-before: #ifdef NGS_PYTHON

We can make our CF accessible as a Python class in the usual way:

.. literalinclude:: myCoefficient.hpp
   :start-at: py::class
   :end-at: .def
