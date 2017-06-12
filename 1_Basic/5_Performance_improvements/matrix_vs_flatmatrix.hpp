#ifndef FILE_MATRIX_VS_FLATMATRIX_HPP
#define FILE_MATRIX_VS_FLATMATRIX_HPP

/*
  This file contains the same integrator as in ../1_myFEM/myIntegrator but using
  FlatMatrices on the LocalHeap instead of Matrices with memory management
 */


namespace ngfem
{
  class MyBetterLaplaceIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_lambda;
  public:
    MyBetterLaplaceIntegrator (shared_ptr<CoefficientFunction> coef)
      : coef_lambda(coef) { ; }

    virtual string Name () const { return "MyLaplace"; }

    virtual int DimElement () const { return 2; }
    virtual int DimSpace () const { return 2; }
    virtual bool IsSymmetric () const { return true; }

    // a volume integral (ngsolve 6.2)
    virtual VorB VB() const { return VOL; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & base_fel,
                       const ElementTransformation & eltrans,
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const
    {
      const ScalarFiniteElement<2> & fel =
        dynamic_cast<const ScalarFiniteElement<2> &> (base_fel);
      int ndof = fel.GetNDof();
      elmat = 0;

      // Allocate dshape matrices on the LocalHeap
      FlatMatrix<> dshape_ref(ndof, 2, lh);
      FlatMatrix<> dshape(ndof, 2, lh);
      IntegrationRule ir(fel.ElementType(), 2*fel.Order());

      for (int i = 0 ; i < ir.GetNIP(); i++)
        {
          MappedIntegrationPoint<2,2> mip(ir[i], eltrans);
          double lam = coef_lambda -> Evaluate (mip);
          fel.CalcDShape (ir[i], dshape);
          dshape = dshape_ref * mip.GetJacobianInverse();
          double fac = mip.IP().Weight() * mip.GetMeasure();
          elmat += (fac*lam) * dshape * Trans(dshape);
        }
    }
  };
}

#endif // FILE_MATRIX_VS_FLATMATRIX_HPP
