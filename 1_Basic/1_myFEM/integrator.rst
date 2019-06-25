Custom Integrators
==================

Sometimes new `DifferentialOperators` and `CoefficientFunctions` are not enough and one needs to implement special element matrices or vectors. This can be done by creating `BilinearFormIntegrators` and `LinearFormIntegrators`.

The integrators need to implement some informations about what dimension it is defined on, if it is a boundary or volume integrator and so on. 

.. literalinclude:: myIntegrator.hpp
   :language: cpp
   :start-at: class MyLaplaceIntegrator
   :end-at: };

The main function to be implemented is `CalcElementMatrix`, in the implementation the element matrix is constructed by creating an integration rule and iterating over it using the `CalcDShape` function of our elements.

.. literalinclude:: myIntegrator.cpp
   :language: cpp
   :start-at: void MyLaplaceIntegrator ::
   :end-before: /*\n    Calculates

For a `BilinearFormIntegrator` the main function to be implemented is
