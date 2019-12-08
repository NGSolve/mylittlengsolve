#ifndef FILE_MYINTEGRATOR_HPP
#define FILE_MYINTEGRATOR_HPP

/*********************************************************************/
/* File:   myIntegrator.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple integrators for the Poisson Equation

*/


namespace ngfem
{
  
  // integrator for \int \lambda(x) \nabla u \nabla v dx
  class MyLaplaceIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_lambda;
  public:
    MyLaplaceIntegrator (shared_ptr<CoefficientFunction> coef)
      : coef_lambda(coef) { ; }

    string Name () const override { return "MyLaplace"; }

    int DimElement () const override { return 2; }
    int DimSpace () const override { return 2; }
    xbool IsSymmetric () const override { return true; }
    VorB VB() const override { return VOL; }

    // Calculates the element matrix
    void CalcElementMatrix (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatMatrix<double> elmat,
                            LocalHeap & lh) const override;
  };

  // integrator for \int f v dx
  class MySourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_f;
  public:
    MySourceIntegrator (shared_ptr<CoefficientFunction>coef)
      : coef_f(coef)
    { ; }

    string Name () const override { return "MySource"; }

    bool BoundaryForm () const override { return false; }
    VorB VB() const override { return VOL; }
    
    // Calculates the right hand side element vector
    void CalcElementVector (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatVector<double> elvec,
                            LocalHeap & lh) const override;
  };
}

#include <python_ngstd.hpp>
void ExportMyIntegrator(py::module m);
#endif
