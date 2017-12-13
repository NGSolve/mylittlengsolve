#ifndef __FILE_MYCOEFFICIENT_HPP
#define __FILE_MYCOEFFICIENT_HPP

namespace ngfem
{
  class MyCoefficientFunction : public CoefficientFunction
  {
    string filename;
  public:
    MyCoefficientFunction(string afilename)
      : CoefficientFunction(/*dimension = */ 1),filename(afilename) { ; }

    virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
    {
      ifstream is(filename);
      double val;
      is >> val;
      return val;
    }
  };
}

void ExportMyCoefficient(py::module m)
{
  py::class_<MyCoefficientFunction, shared_ptr<MyCoefficientFunction>, CoefficientFunction>
    (m, "MyCoefficient", "CoefficientFunction that reads value from file")
    .def(py::init<string>())
    ;
}

#endif //  __FILE_MYCOEFFICIENT_HPP
