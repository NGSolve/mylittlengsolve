#ifndef __FILE_MYCOEFFICIENT_HPP
#define __FILE_MYCOEFFICIENT_HPP

#include <comp.hpp>

namespace ngcomp
{
  class MyCoefficientFunction : public CoefficientFunction
  {
  public:
    MyCoefficientFunction()
      : CoefficientFunction(/*dimension = */ 1) { ; }

    // ********************* Necessary **********************************
    virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
    {
      static Timer timer("MyCoefficientFunction - Evaluate(BaseMappedIntegrationPoint)");
      RegionTimer reg_timer(timer);
      // maybe do something with mip and trafo?
      auto pnt = mip.GetPoint();
      return pnt[0] * pnt[1];
    }

    // ******************* Performance improvements **********************
  };

  void ExportMyCoefficient(py::module m);

  inline double MyIntegrate( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    LocalHeap lh(1000000, "localheap");

    double hsum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      auto element = ma->GetElement(ElementId(VOL, elnum));
      auto & trafo = ma->GetTrafo (element, lh);
      IntegrationRule ir(trafo.GetElementType(), order);
      BaseMappedIntegrationRule & mir = trafo(ir, lh);
      Matrix<> values(ir.Size(), 1);
      cf -> Evaluate (mir, values);
      for (int i = 0; i < values.Height(); i++)
        hsum += mir[i].GetWeight() * values(i,0);
    }
    return hsum;
  }
}


#endif //  __FILE_MYCOEFFICIENT_HPP
