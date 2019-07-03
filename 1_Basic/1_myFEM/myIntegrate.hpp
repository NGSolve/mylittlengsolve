#ifndef __FILE_MYINTEGRATE_HPP
#define __FILE_MYINTEGRATE_HPP

#include <comp.hpp>

namespace ngcomp
{

  /*

  inline double MyIntegrate( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    LocalHeap lh(1000000, "localheap");

    double sum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      HeapReset hr(lh);
      auto element = ma->GetElement(ElementId(VOL, elnum));
      auto & trafo = ma->GetTrafo (element, lh);
      IntegrationRule ir(trafo.GetElementType(), order);
      Matrix<> values(ir.Size(), 1);

      BaseMappedIntegrationRule & mir = trafo(ir, lh);
      cf -> Evaluate (mir, values);

      for (int i = 0; i < values.Height(); i++)
        sum += mir[i].GetWeight() * values(i,0);
    }
    return sum;
  }
  */

  // Version including Timers
  inline double MyIntegrate( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    static Timer t_all("MyIntegrate");
    static Timer t_evaluate("Evaluate");
    static Timer t_sum("Sum");
    static Timer t_after_init("After init");

    RegionTimer r_all(t_all);

    LocalHeap lh(1000000, "localheap");

    double sum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      HeapReset hr(lh);
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
          sum += mir[i].GetWeight() * values(i,0);
      }
    }
    return sum;
  }

  // SIMD Version
  inline double MyIntegrateSIMD( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    static Timer t_all("MyIntegrateSIMD");
    static Timer t_evaluate("Evaluate");
    static Timer t_sum("Sum");
    static Timer t_after_init("After init");

    RegionTimer r_all(t_all);

    LocalHeap lh(100000000, "localheap");

    SIMD<double> hsum = 0.0;
    for (auto elnum : Range(ma->GetNE(VOL)))
    {
      HeapReset hr(lh);
      RegionTimer r_int_el(t_all);
      auto element = ma->GetElement(ElementId(VOL, elnum));
      auto & trafo = ma->GetTrafo (element, lh);
      SIMD_IntegrationRule ir(trafo.GetElementType(), order);
      Matrix<SIMD<double>> values(1, ir.Size());

      RegionTimer r_after_init( t_after_init);
      SIMD_BaseMappedIntegrationRule & mir = trafo(ir, lh);
      {
        RegionTimer r_evaluate( t_evaluate);
        cf -> Evaluate (mir, values);
      }

      {
        RegionTimer r_sum( t_sum);
        for (int i = 0; i < values.Width(); i++)
          hsum += mir[i].GetWeight() * values(0,i);
      }
    }
    return HSum(hsum);
  }

  // Parallel Version
  inline double MyIntegrateParallel( shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, int order )
  {
    static Timer t_all("MyIntegrateParallel");
    static Timer t_evaluate("Evaluate");
    static Timer t_sum("Sum");
    static Timer t_after_init("After init");

    RegionTimer r_all(t_all);

    LocalHeap glh(1000000, "localheap");

    atomic<double> sum = 0.0;

    ParallelForRange ( Range(ma->GetNE(VOL)), [&] (auto myrange)
    {
      auto lh = glh.Split();
      SIMD<double> local_sum = 0.0;
      for (auto elnum : myrange)
      {
        HeapReset hr(lh);
        RegionTracer r_int_el(TaskManager::GetThreadId(), t_all);
        auto element = ma->GetElement(ElementId(VOL, elnum));
        auto & trafo = ma->GetTrafo (element, lh);
        SIMD_IntegrationRule ir(trafo.GetElementType(), order);
        FlatMatrix<SIMD<double>> values(1, ir.Size(), lh);

        RegionTracer r_after_init(TaskManager::GetThreadId(),  t_after_init);
        SIMD_BaseMappedIntegrationRule & mir = trafo(ir, lh);
        {
          RegionTracer r_evaluate(TaskManager::GetThreadId(),  t_evaluate);
          cf -> Evaluate (mir, values);
        }

        {
          RegionTracer r_sum(TaskManager::GetThreadId(),  t_sum);
          for (int i = 0; i < values.Width(); i++)
            local_sum += mir[i].GetWeight() * values(0,i);
        }
      }
      sum += HSum(local_sum);
    }, 100);

    return sum;
  }
}

#endif // __FILE_MYINTEGRATE_HPP
