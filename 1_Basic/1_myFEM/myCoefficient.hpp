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
      RegionTracer r_timer(TaskManager::GetThreadId(), timer);
      // maybe do something with mip and trafo?
      auto pnt = mip.GetPoint();
      return pnt[0] * pnt[1];
    }

    // ******************* Performance improvements **********************

    virtual void Evaluate (const BaseMappedIntegrationRule & mir, BareSliceMatrix<double> hvalues) const override
    {
      static Timer timer("MyCoefficientFunction - Evaluate(BaseMappedIntegrationRule)");
      RegionTracer r_timer(TaskManager::GetThreadId(), timer);
      auto values = hvalues.AddSize(mir.Size(), Dimension());
      auto points = mir.GetPoints();
      for (int i = 0; i < mir.Size(); i++)
        values(i,0) = points(i,0) * points(i,1);
    }

  };

  void ExportMyCoefficient(py::module m);

  inline double MyIntegrate( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    LocalHeap lh(100000000, "localheap");

    double hsum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      auto element = ma->GetElement(ElementId(VOL, elnum));
      auto & trafo = ma->GetTrafo (element, lh);
      IntegrationRule ir(trafo.GetElementType(), order);
      Matrix<> values(ir.Size(), 1);

      BaseMappedIntegrationRule & mir = trafo(ir, lh);
      cf -> Evaluate (mir, values);

      for (int i = 0; i < values.Height(); i++)
        hsum += mir[i].GetWeight() * values(i,0);
    }
    return hsum;
  }
}

  /*
  // Version including Timers
  inline double MyIntegrate( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    static Timer t_all("MyIntegrate");
    static Timer t_evaluate("Evaluate");
    static Timer t_sum("Sum");
    static Timer t_after_init("After init");

    RegionTimer r_all(t_all);

    LocalHeap lh(100000000, "localheap");

    double hsum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      RegionTimer r_int_el(t_all);
      auto element = ma->GetElement(ElementId(VOL, elnum));
      auto & trafo = ma->GetTrafo (element, lh);
      IntegrationRule ir(trafo.GetElementType(), order);
      Matrix<> values(ir.Size(), 1);

      RegionTimer r_after_init( t_after_init);
      BaseMappedIntegrationRule & mir = trafo(ir, lh);
      {
        RegionTimer r_evaluate( t_evaluate);
        cf -> Evaluate (mir, values);
      }

      {
        RegionTimer r_sum( t_sum);
        for (int i = 0; i < values.Height(); i++)
          hsum += mir[i].GetWeight() * values(i,0);
      }
    }
    return hsum;
  }
  */
}


#endif //  __FILE_MYCOEFFICIENT_HPP
