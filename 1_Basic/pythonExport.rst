Cmake & Python Export
=======================

Cmake
-------

We use `cmake <https://cmake.org/>`_ to build NGSolve, therefore we recommend, that you set up your
extension using cmake as well.
First we add NGSolve as a dependency and give the default install dirs as hints:

.. literalinclude:: CMakeLists.txt
   :language: cmake
   :end-before: add_ngsolve

then we create our extension library `myngspy` which contains all the cpp files we will create in this
tutorial:

.. literalinclude:: CMakeLists.txt
   :language: cmake
   :start-at: add_ngsolve
   :end-before: # check if

We recommend, that you set the default install (if `CMAKE_INSTALL_PREFIX` is not set by the user)
dir to the NGSolve install dir, then the include paths and the `PYTHONPATH` will be set correctly:

.. literalinclude:: CMakeLists.txt
   :language: cmake
   :start-at: if(CMAKE_INSTALL

Thats all you need to build your extension! :)


Python Export
---------------

We use `Pybind11 <https://github.com/pybind/pybind11>`_ to export C++ classes to Python. With pybind we create a Python module with wrapper classes and functions to access the C++ code. Note that this is only a very short introduction and focused on NGSolve classes, for a more detailed introduction into C++ binding have a look at the `Pybind11 documentation <http://pybind11.readthedocs.io/en/latest/>`_ .
In your project there must be one .cpp file defining the macro ``PYBIND11_MODULE`` with the module name and the variable the module will be stored in.
We recommend importing the NGSolve module here, then the dependencies on the NGSolve classes will be available.

.. code-block:: C++

   PYBIND11_MODULE(myngspy,m) {
     // import ngsolve such that python base classes are defined
     py::module::import("ngsolve");

     /*
     all the functionality you want to export...
     */
   }

Inside this macro we export all the functions and classes we need.
Functions can be defined in two ways, first, if all the arguments and the return value can be directly converted to Python objects (They are either references to exported types or ``std::shared_ptr`` to an exported type with holder type ``std::shared_ptr``, or raw pointer to exported types) then the function can simply exported by handing pybind the function pointer.

.. code-block:: C++

   m.def("SomeFunction", &SomeCppFunction, "some documentation")

the documentation string will be displayed if you call

.. code-block:: Python

   help(SomeFunction)

in Python.

If the arguments can't be directly parsed, or some pre-/post-processing needs to be done before they can be parsed to the C++ function it is possible to do this in a lambda function:

.. code-block:: C++

   m.def("SomeFunction", [] (First & first,
                             shared_ptr<Second> second,
                             py::object obj)
     {
       shared_ptr<MyExportedObject> myObject;
       if(py::is_none(obj))
         myObject = do_some_stuff(first, second);
       else
       {
         auto third = py::cast<shared_ptr<Third>>(obj);
         myObject = do_some_other_stuff(first, second, third);
       }
       return myObject;
     }, arg("first"), arg("second"), arg("third") = py::none(),
     "some documentation");

Note, we can either pass references, shared_ptr or pybind11 objects to our function. Further we can use utility functions from pybind, like py::cast before calling our C++ functions. The shared_ptr return type will be automatically converted to a Python object sharing lifetime with the C++ object.
We can additionally name our arguments, so they can be called from Python as keyword arguments and we can prescribe default arguments.

But let's start with some NGSolve stuff, as a first example we will create our own CoefficientFunction ->
