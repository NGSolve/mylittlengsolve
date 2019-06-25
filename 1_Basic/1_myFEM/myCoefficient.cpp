
#include <fem.hpp>
#include <python_ngstd.hpp>
#include "myCoefficient.hpp"

namespace ngfem
{
  // Export cf to Python
  void ExportMyCoefficient(py::module m)
  {
    py::class_<MyCoefficientFunction, shared_ptr<MyCoefficientFunction>, CoefficientFunction>
      (m, "MyCoefficient", "CoefficientFunction that returns x*y.")
      .def(py::init<>())
      ;
  }

  // Register cf for pickling/archiving
  // Create static object with template parameter function and base class.
  static RegisterClassForArchive<MyCoefficientFunction, CoefficientFunction> regmycf;
}
