Python Export
===============

We use `Pybind11 <https://github.com/pybind/pybind11>`_ to export C++ classes to Python. With pybind we create a Python module with wrapper classes and functions to access the C++ code. Note that this is only a very short introduction and focused on NGSolve classes, for a more detailed introduction into C++ binding have a look at the `Pybind11 documentation <http://pybind11.readthedocs.io/en/latest/>`_ .
In your project there must be one .cpp file defining the macro ``PYBIND11_PLUGIN`` with the module name and define the module.
We recommend importing the NGSolve module here, if you don't do this here, the NGSolve module must be imported before your module all the time to have the Python base classes available.

.. code-block:: C++

   PYBIND11_PLUGIN(myngspy) {
     // import ngsolve such that python base classes are defined
     py::module::import("ngsolve");
     py::module m("myngspy", "myngspy doc-string");

     return m.ptr();
   }

In between the module definition and the return statement we export all the functions and classes.
Functions can be defined in two ways, first if all the arguments and the return value can be directly converted to Python objects (They are either references to exported types or ``std::shared_ptr`` to an exported type with holder type ``std::shared_ptr`` then the function can simply exported by giving pybind the function pointer.

.. code-block:: C++

   m.def("SomeFunction", &SomeCppFunction, "some documentation")

the documentation string will be displayed if you call

.. code-block:: Python

   help(SomeFunction)

in Python.

If the arguments can't be directly parsed, or some pre-/post-processing needs to be done before they can be parsed to the C++ function we have to do this in a lambda function

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

As you can see, we can either pass references, shared_ptr or pybind11 objects to our function, then we can use utility functions from pybind, like py::cast and call our C++ functions. The shared_ptr return type will be automatically converted to a Python object sharing lifetime with the C++ object.
We can additionally name our arguments, so they can be called from Python as keyword arguments and we can prescribe default arguments.

