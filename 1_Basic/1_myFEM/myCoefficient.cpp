
#include <comp.hpp>
#include <python_ngstd.hpp>
#include "myCoefficient.hpp"

namespace ngcomp
{
  // Export cf to Python
  void ExportMyCoefficient(py::module m)
  {
    py::class_<MyCoefficientFunction, shared_ptr<MyCoefficientFunction>, CoefficientFunction>
      (m, "MyCoefficient", "CoefficientFunction that returns x*y.")
      .def(py::init<>())
      ;

    m.def("MyIntegrate", &MyIntegrate, 
	py::arg("cf"),
        py::arg("mesh"),
        py::arg("order")=5
        );
  }

  // Register cf for pickling/archiving
  // Create static object with template parameter function and base class.
  static RegisterClassForArchive<MyCoefficientFunction, CoefficientFunction> regmycf;
}
