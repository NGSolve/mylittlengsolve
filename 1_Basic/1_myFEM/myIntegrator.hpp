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

    virtual string Name () const { return "MyLaplace"; }

    virtual int DimElement () const { return 2; }
    virtual int DimSpace () const { return 2; }
    virtual xbool IsSymmetric () const { return true; }

    // a volume integral (ngsolve 6.2)
    virtual VorB VB() const { return VOL; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans, 
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;
  };

  // integrator for \int f v dx
  class MySourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_f;
  public:
    MySourceIntegrator (shared_ptr<CoefficientFunction>coef)
      : coef_f(coef)
    { ; }

    virtual string Name () const { return "MySource"; }

    virtual bool BoundaryForm () const { return false; }
    virtual VorB VB() const { return VOL; }
    
    // Calculates the right hand side element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;
  };
}

#include <python_ngstd.hpp>
void ExportMyIntegrator(py::module m);
#endif
