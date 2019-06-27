Custom Integrators
==================

Sometimes new `DifferentialOperators` and `CoefficientFunctions` are not enough and one needs to implement special element matrices or vectors. This can be done by creating `BilinearFormIntegrators` and `LinearFormIntegrators`.

The integrators need to implement some informations about what dimension it is defined on, if it is a boundary or volume integrator and so on. 

.. literalinclude:: myIntegrator.hpp
   :language: cpp
   :start-at: class MyLaplaceIntegrator
   :end-at: };

For a `BilinearFormIntegrator` the main function to be implemented is the construction of element matrices. This is done by creating an integration rule and iterating over it using the `CalcDShape` function of our elements.

.. literalinclude:: myIntegrator.cpp
   :language: cpp
   :start-at: void MyLaplaceIntegrator ::
   :end-before: Calculates the element vector

Similarily a `LinearFormIntegrator` has to implement the construction of element vectors:

.. literalinclude:: myIntegrator.cpp
   :language: cpp
   :start-at: void MySourceIntegrator ::
   :end-before: void ExportMyIntegrator

Again we export the integrators to Python:

.. literalinclude:: myIntegrator.cpp
   :language: cpp
   :start-at: void ExportMyIntegrator
