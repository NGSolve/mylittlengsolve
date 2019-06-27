NGSolve Basics
=================

In section one of `mylittlengsolve <https://github.com/NGSolve/mylittlengsolve>`_ we want to show how
you can modify NGSolve classes to your needs and export them to Python for proper use.
This is not intended as an introduction to C++ programming, but to guide you
to the C++ parts of NGSolve. Some basic concepts of object orientated programming, generic programming
(templates and that stuff) will make life a lot easier here ;) It is as well recommended to have
a basic knowledge about C++11 features like smart pointers and lambda functions.

Later parts of this tutorial may depend on earlier parts so it's recommended to read it in order,
however with some C++ experience it should be ok to start anywhere.

Build instructions
------------------

This project is built using CMake, just create a build directory, execute cmake and make:

.. code-block:: bash

   mkdir build
   cd build
   cmake ..
   make -j
   make install

On Windows you need to have Visual Studio/build tools (at least 2017) installed:

.. code-block:: bash

   mkdir build
   cd build
   cmake .. -G "Visual Studio 15 2017 Win64"
   cmake --build . --config release --target install --config RelWithDebInfo

As default mylittlengsolve is installed where your NGSolve installation is. If you want to change
that behaviour, set the `CMAKE_INSTALL_PREFIX` variable:

.. code-block:: bash

   cmake .. -DCMAKE_INSTALL_PREFIX=your_install_path

In this chapter we will look into how to implement new CoefficientFunctions and FESpaces in NGSolve and export them to Python. For the spaces we will need to create FiniteElements and DifferentialOperators.

.. toctree::
   :maxdepth: 1
   :caption: Basic

   pythonExport.rst

   1_myFEM/coefficient.rst

   1_myFEM/element.rst

   1_myFEM/diffop.rst

   1_myFEM/fespace.rst

   1_myFEM/integrators.rst

   2_myHOFEM/hofem.rst

