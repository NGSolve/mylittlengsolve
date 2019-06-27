My DiffOp
=========

DiffOps are the static polymorphism classes implementing Operators for GridFunctions and Test- and TrialFunctions. They only have static member functions and are never created but only used as a template argument for DifferentialOperators. They use `CRTP <https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern>`_ for "static virtual" methods. DifferentialOperator is the runtime polymorphism class wrapping them. Use T_DifferentialOperator<DiffOp> as a coupling mechanism.

We will need the identity and the gradient operator. These operators have to define the following constexpr values:

==============  ===========================
 DIM              Dimension of input
 DIM_SPACE        Dimension of FESpace
 DIM_ELEMENT      Dimension of Element
 DIM_DMAT         Dimension of range
 DIFFORDER        Order of differentiation
==============  ===========================

The order of differentiation is used to choose the correct integration rules for integrating the polynomials. Further we need to implement a function that generates the result of the Diffop for one mapped integration point. This function gets a matrix because the result could be matrix valued. A minimal implementation of the identity would be this:

.. literalinclude:: myDiffOp.hpp
   :language: cpp
   :start-at: class MyDiffOpId
   :end-before:     // ******************** Performance improvements *******************
   :append:
        };
      }

Again some more functions can be implemented to improve performance, for example evaluation in the whole integration rule and SIMD implementations.

The gradient operator looks similar:

.. literalinclude:: myDiffOp.hpp
   :language: cpp
   :start-at: class MyDiffOpGradient
   :end-before:     // ******************** Performance improvements *******************
   :append:
        };
      }

Again the implementation of a boundary evaluation operator will be left as an exercise.

To make the DiffOp available at runtime we have to put it into a T_DifferentialOperator class as template argument. This happens next in the FESpace implementation.
