#ifndef FILE_MECHSOLVE
#define FILE_MECHSOLVE

/*********************************************************************/
/* File:   mechsolve.hpp                                             */
/* Author: START                                                     */
/* Date:   25. Mar. 2003                                             */
/*********************************************************************/


using namespace ngsolve;


template <int D>
class NonlinMechIntegrator : public BilinearFormIntegrator
{
  shared_ptr<CoefficientFunction> coeff_e;
  shared_ptr<CoefficientFunction> coeff_nu;
public:
  NonlinMechIntegrator (const Array<shared_ptr<CoefficientFunction>> & coefs)
    : coeff_e(coefs[0]), coeff_nu(coefs[1])
  { ; }

  virtual ~NonlinMechIntegrator ()
  { ; }
 
  virtual bool BoundaryForm () const { return 0; }
  virtual int DimFlux () const { return D*(D+1)/2; }
  virtual int DimElement () const { return D; }
  virtual int DimSpace () const { return D; }



  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans,
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    Vector<double> vec0 (D*fel.GetNDof());
    vec0 = 0;
    CalcLinearizedElementMatrix (fel, eltrans, vec0, elmat, lh);
  }

  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> & elveclin,
			       FlatMatrix<double> & elmat,
			       LocalHeap & lh) const;
  
  virtual void
  ApplyElementMatrix (const FiniteElement & bfel,
		      const ElementTransformation & eltrans,
		      const FlatVector<double> & elx,
		      FlatVector<double> & ely,
		      void * precomputed,
		      LocalHeap & lh) const;

  virtual void
  ApplyLinearizedElementMatrix (const FiniteElement & fel,
				const ElementTransformation & eltrans,
				const FlatVector<double> & ellin,
				const FlatVector<double> & elx,
				FlatVector<double> & ely,
				LocalHeap & lh) const
  {
    ApplyElementMatrix (fel, eltrans, elx, ely, 0, lh);
  }

  virtual double Energy (const FiniteElement & bfel,
			 const ElementTransformation & eltrans,
			 const FlatVector<double> & elx,
			 LocalHeap & lh) const;


  virtual void
  CalcFlux (const FiniteElement & fel,
	    const BaseMappedIntegrationPoint & bmip,
	    const FlatVector<double> & elx, 
	    FlatVector<double> & flux,
	    bool applyd,
	    LocalHeap & lh) const;
};



#endif
