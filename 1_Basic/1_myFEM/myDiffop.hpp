#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

namespace ngfem
{
  class MyPow3 : public DifferentialOperator
  {
  public:
    MyPow3() : DifferentialOperator(1,1,VOL,0) { ; }

    virtual void Apply(const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   FlatVector<double> x,
	   FlatVector<double> pow3,
	   LocalHeap & lh) const
    {
      double u = static_cast<const ScalarFiniteElement<2>&> (fel).Evaluate(mip.IP(),x);
      pow3(0) = u*u*u;
    }
  };
}


#ifdef NGS_PYTHON
#include <python_ngstd.hpp>
void ExportMyDiffOp(py::module m)
{
  py::class_<MyPow3, shared_ptr<MyPow3>, DifferentialOperator>
    (m, "MyPow3")
    .def(py::init<>())
    ;
}
#endif // NGS_PYTHON


#endif // FILE_MYDIFFOP_HPP
