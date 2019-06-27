
High Order Spaces
===================

Creating the Element
--------------------

For higher order spaces most of the times some type of continuity of the bases functions over edges or faces must be fullfilled. Because of this the order of the basis functions, for example on an edge is important. So the edges and faces must be oriented consistently. For this NGSolve provides the interface class `VertexOrientedFE<ELEMENT_TYPE>` defined in ``fem/fe_interfaces.hpp``.

This class stores the vertex numbers of the element and when asked for the vertices of an edge it returns them sorted according to the global orientation. With this information we can build basis functions which are consistent over element boundaries. The definition of the segment looks similar to the definition of low order elements:

.. literalinclude:: myHOElement.hpp
   :language: cpp
   :start-at: class MyHighOrderSegm
   :end-at: };

We will use the `AutoDiff` class from the implementation of the second order elements again to not have to compute the gradients of the shape function manually. For this we define the template function `T_CalcShape` which can either take an arbitrary type (we will use `double` and `AutoDiff`) and compute the value of the basis function.

For the segment the first two basis functions are the vertex ("hat") functions, then the higher order edge basis functions (integrated legendre polynomials):

.. literalinclude:: myHOElement.cpp
   :language: cpp
   :start-at: void MyHighOrderSegm :: T_CalcShape
   :end-before: MyHighOrderTrig
   :prepend: template <class T>

`CalcShape` and `CalcDShape` just use this helper function:

.. literalinclude:: myHOElement.cpp
   :language: cpp
   :start-at: void MyHighOrderSegm :: CalcShape
   :end-before: template <class T>

The implementation of trigs is similar, only inner dofs (bubbles) are added:

.. literalinclude:: myHOElement.cpp
   :language: cpp
   :start-at: void MyHighOrderTrig :: T_CalcShape
   :prepend: template <class T>

Creating the FESpace
---------------------

The definition of the finite element space looks similar to the low order case, we only store an array of the dof numbers of edge and cell dofs:

.. literalinclude:: myHOFESpace.hpp
   :language: cpp
   :start-at: class MyHighOrderFESpace
   :end-at: };

In the `Update` function we fill this array and compute some global variables:

.. literalinclude:: myHOFESpace.cpp
   :language: cpp
   :start-at: void MyHighOrderFESpace :: Update
   :end-before: void MyHighOrderFESpace :: GetDofNrs

The `GetDofNrs` function collects the dof nrs into the output array:

.. literalinclude:: myHOFESpace.cpp
   :language: cpp
   :start-at: void MyHighOrderFESpace :: GetDofNrs
   :end-before: FiniteElement &

When returning the finite element we need to set the vertex numbers, so that `VertexOrientedFE` has the correct orientation:

.. literalinclude:: myHOFESpace.cpp
   :language: cpp
   :start-at: FiniteElement & MyHighOrderFESpace
   :end-before: static RegisterFESpace

Finally we have to register the new space to the NGSolve register:

.. literalinclude:: myHOFESpace.cpp
   :language: cpp
   :start-at: static RegisterFESpace
   :end-at: static RegisterFESpace
