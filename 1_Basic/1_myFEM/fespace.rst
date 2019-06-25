My FESpace
===========

Next we want to use our created finite elements and put them into a space. While the finite element
handles all the local element wise computations, the space knows the global mapping of the dofs,
has access to the mesh and connects local to global computations.

Again there are only a few virtual function overrides necessary to get our space working:

.. literalinclude:: myFESpace.hpp
   :language: cpp
   :start-at: class MyFESpace
   :end-at: };

In the constructor we query the flags (if we want a second order space) and define evaluators for
functions on the space, their derivatives and
their boundary values (which will be left as an exercise):

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: MyFESpace :: MyFESpace
   :end-before: DocInfo MyFESpace

In the ``Update`` function, the number of degrees of freedom are set:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: void MyFESpace :: Update
   :end-before: void MyFESpace :: GetDofNrs

Further we must define the mapping from local dofs for each ``ElementId`` to the global
dof numbers:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: void MyFESpace :: GetDofNrs
   :end-before: FiniteElement & MyFESpace

as well as a function that returns our new element:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: FiniteElement & MyFESpace
   :end-before: /*

NGSolve has an internal register of finite element spaces, where we need to register our new space:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: static Register
   :end-at: static Register

To register FESpace for Python we provide the helper function `ExportFESpace` in `python_comp.hpp`. This automatically creates the constructor and pickling support for the space. The function returns the Python class and additional functionality can be appended.

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: ExportFESpace
   :end-at: ;

If we want to use the `secondorder` flag now, we would get a warning, that we use an undocumented flag. This is a safety guide to prevent typos for keyword arguments, what we need to do is to provide documentation for our flag:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: DocInfo MyFESpace
   :end-before: void MyFESpace :: 

..
   TODO: Remove the derived flags that are not available.

..
   If we want to have a nice docstring including our flag we have to call another helper function that searches our module for all classes and appends the docstring:

..
   .. literalinclude:: ../myngspy.cpp
      :language: cpp
      :start-at: ngs.attr("_add_flags_doc")(m);
      :end-at: ngs.attr

Next we will show how to implement custom integrators.
