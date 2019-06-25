#ifndef __FILE_MYCOEFFICIENT_HPP
#define __FILE_MYCOEFFICIENT_HPP

namespace ngfem
{
  class MyCoefficientFunction : public CoefficientFunction
  {
  public:
    MyCoefficientFunction()
      : CoefficientFunction(/*dimension = */ 1) { ; }

    // ********************* Necessary **********************************
    virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
    {
      // maybe do something with mip and trafo?
      auto pnt = mip.GetPoint();
      return pnt[0] * pnt[1];
    }

    // ******************* Performance improvements **********************
  };

  void ExportMyCoefficient(py::module m);
}


#endif //  __FILE_MYCOEFFICIENT_HPP
