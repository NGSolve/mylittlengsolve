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
their boundary values:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: void MyFESpace :: Update
   :end-at: }

In the ``Update`` function, the degrees of freedom are set:

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

Registering you FESpace for Python is quite similar as before, only the init (optional) function
is different:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: auto myfes
   :end-at: }), py::arg

The init function takes the mesh as argument and all the flags as keyword arguments,
then we call the ``CreateFlagsFromKwArgs`` function with the mesh as additional information
(you can just copy paste this from the base class constructor). This is necessary, because the
base FESpace may have some flags we don't know anything about (and do not have to care about)
in the derived class, i.e. dirichlet boundaries or it is defined on only a part of the domain.
Then we create our new space and call the ``__initialize__`` function from the base class,
which calls amongst other things the ``Update`` function of our space.
Instead of this init function you can create the registered FESpace in Python as well just with

.. code:: python

   FESpace("myfespace", mesh, ...)

with the name we registered it.

If we want to use the `secondorder` flag now, we would get a warning, that we use an undocumented
flag. This is a safety guide to prevent typos, so we want to register our flag for the space:

.. literalinclude:: myFESpace.cpp
   :language: cpp
   :start-at: .def_static
   :end-before: .def("GetNVert"

Here we query the flags_doc from the base class and append our secondorder flag to the dictionary.

