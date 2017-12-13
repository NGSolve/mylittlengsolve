My CoefficientFunction
========================

Createing your own CoefficientFunctions is probably the most common thing someone needs to do when
you have to dive into NGSolve C++. There's a lot you can do from Python, but sometimes you want
something special that's not possible. This tutorial shows how to create arbitrary CoeffientFunctions.

For it's basic functionality a CF does not need much: Derive from the base class, provide a
constructor and override the pure virtual method `Evaluate`. The CF in our example will read it's
value (live) from a given file.

A simple implementation of the CF would look like this:

.. literalinclude:: myCoefficient.hpp
   :language: cpp
   :start-after: #define __FILE_MYCOEFFICIENT_HPP
   :end-before: void ExportMyCoefficient

Note that there are a lot of other virtual functions a CF can override, mostly for performance
(i.e. evaluate in all IntegrationPoints of an IntegrationRule at once or use AVX vectorization
with our SIMD class), but we want to keep it simple for now.

Next we want to make our CF accessible from Python. Therefore we use pybind11 to export the class
by createing a ``py::class_`` instance. This template class takes the exported class, then the
holder type (a shared pointer to the exported class) and the base class(es) as template
arguments. The constructor takes the module where we export it, it's Python name and
the docstring. Further we define an ``__init__`` function, where we use the default pybind
``py::init`` method, which just calls the class constructor with the same arguments we give as
template arguments. Note that the semicolon is in the next line: ``.def`` always returns the
object again, so we can use it to define multiple functions for the class and then after all
definitions close the statment with the semicolon.

.. literalinclude:: myCoefficient.hpp
   :language: cpp
   :start-at: py::class
   :end-at: .def

Our next step will be to create our first finite element.
