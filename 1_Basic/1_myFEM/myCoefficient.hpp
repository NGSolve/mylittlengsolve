#ifndef __FILE_MYCOEFFICIENT_HPP
#define __FILE_MYCOEFFICIENT_HPP

#include <comp.hpp>

namespace ngcomp
{
  template <typename T>
  auto myFunction( T point )
  {
    auto x = point[0];
    auto y = point[1];
    auto z = point[2];

    return x*y;

//     auto res = x*x+y*y+z*z;
// 
//     for (auto i : Range(100))
//       res = (res+1.0)/res;
// 
//     return res*x*y*z;
  }

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
      return myFunction(pnt);
    }

    // ******************* Performance improvements **********************

    /*
    // Evaluate whole integration rule at once (saves function calls)
    virtual void Evaluate (const BaseMappedIntegrationRule & mir, BareSliceMatrix<double> hvalues) const override
    {
      static Timer timer("MyCoefficientFunction - Evaluate(BaseMappedIntegrationRule)");
      RegionTracer r_timer(TaskManager::GetThreadId(), timer);
      auto values = hvalues.AddSize(mir.Size(), Dimension());
      auto points = mir.GetPoints();
      for (int i = 0; i < mir.Size(); i++)
        values(i,0) = myFunction(points.Row(i));
    }

    // Use SIMD evaluation (multiple multiplications in one CPU instruction)
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> hvalues) const override
    {
      static Timer timer("MyCoefficientFunction - Evaluate(SIMD_BaseMappedIntegrationRule)");
      RegionTracer r_timer(TaskManager::GetThreadId(), timer);
      auto values = hvalues.AddSize(Dimension(), mir.Size());
      auto points = mir.GetPoints();
      for (int i = 0; i < mir.Size(); i++)
        values(0,i) = myFunction(points.Row(i));
    }
    */
  };

  void ExportMyCoefficient(py::module m);

}


#endif //  __FILE_MYCOEFFICIENT_HPP
