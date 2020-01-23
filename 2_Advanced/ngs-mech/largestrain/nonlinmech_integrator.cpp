// integrator for non-linear elasticity

#include <solve.hpp>
using namespace ngsolve;

template <int D>
class LinearMaterialLaw
{
  shared_ptr<CoefficientFunction> coef_E;
  shared_ptr<CoefficientFunction> coef_nu;

public:
  LinearMaterialLaw (shared_ptr<CoefficientFunction> acoef_E,
		     shared_ptr<CoefficientFunction> acoef_nu)
    : coef_E (acoef_E), coef_nu(acoef_nu) { ; }

  
  double Energy (const MappedIntegrationPoint<D,D> & mip,
		 const Mat<D,D> & strain) const
  {
    double e = coef_E -> Evaluate (mip);
    double nu = coef_nu -> Evaluate (mip);
    double mu = e / 2 / (1+nu);
    double lam = e * nu / ((1+nu)*(1-2*nu));

    return mu*Trace (strain * strain) + lam/2*sqr (Trace(strain));
  }


  // 2nd Piola-Kirchhoff tensor
  Mat<D> CalcStress (const MappedIntegrationPoint<D,D> & mip,
		     const Mat<D,D> & strain) const
  {
    double e = coef_E -> Evaluate (mip);
    double nu = coef_nu -> Evaluate (mip);
    double mu = e / 2 / (1+nu);
    double lam = e * nu / ((1+nu)*(1-2*nu));

    Mat<D> piola2 = (2*mu) * strain + (lam * Trace(strain)) * Id<D>();
    return piola2;
  }


  // directional derivative 
  Mat<D> CalcDStress (const MappedIntegrationPoint<D,D> & mip,
		      const Mat<D,D> & strain,
		      const Mat<D,D> & dstrain) const
  {
    double e = coef_E -> Evaluate (mip);
    double nu = coef_nu -> Evaluate (mip);
    double mu = e / 2 / (1+nu);
    double lam = e * nu / ((1+nu)*(1-2*nu));

    Mat<D> piola2 = (2*mu) * strain + (lam * Trace(strain)) * Id<D>();
    Mat<D> dpiola2 = (2*mu) * dstrain + (lam * Trace(dstrain)) * Id<D>();

    return dpiola2;
  }
  
};



template <int D>
class NonlinMechIntegrator : public BilinearFormIntegrator
{
  LinearMaterialLaw<D> matlaw;
public:
  NonlinMechIntegrator (const Array<shared_ptr<CoefficientFunction>> & coefs)
    : matlaw(coefs[0], coefs[1]) { ; }

  virtual ~NonlinMechIntegrator () { ; }
 
  virtual bool BoundaryForm () const { return 0; }
  virtual int DimFlux () const { return D*(D+1)/2; }
  virtual int DimElement () const { return D; }
  virtual int DimSpace () const { return D; }

  xbool IsSymmetric () const override
  {
    return true;
  }

  VorB VB () const override
  {
    return VOL;
  }

  void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans,
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const override
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










template <int D>
void NonlinMechIntegrator<D> ::
CalcLinearizedElementMatrix (const FiniteElement & bfel,
			     const ElementTransformation & eltrans,
			     FlatVector<double> & elveclin,
			     FlatMatrix<double> & elmat,
			     LocalHeap & lh) const
{
  static Timer timer ("NonLinMech::Calclinearized");
  static Timer timer2 ("NonLinMech::Calclinearized");
  RegionTimer reg (timer);


  const ScalarFiniteElement<D> & fel =
    dynamic_cast<const ScalarFiniteElement<D>&> (bfel);
  
  int ndof = fel.GetNDof();
  
  elmat = 0;
  
  FlatMatrixFixHeight<D*D> bmat (ndof * D, lh);
  FlatMatrixFixHeight<D*D> dbmat (ndof * D, lh);
  Mat<D*D,D*D> dmat;
  
  FlatMatrixFixWidth<D> grad (ndof, lh), gradr(ndof, lh);
  FlatMatrixFixWidth<D> mlin(ndof, &elveclin[0]);
  
  int order = 2 * fel.Order();
  const IntegrationRule & ir = SelectIntegrationRule (fel.ElementType(), order);
  
  for (int i = 0; i < ir.GetNIP(); i++)
    {
      MappedIntegrationPoint<D,D> mip(ir[i], eltrans);
      
      fel.CalcDShape (mip.IP(), gradr);
      grad = gradr * mip.GetJacobianInverse();

      bmat = 0;
      for (int k = 0; k < ndof; k++)
	for (int j=0; j < D; j++)
	  for (int m = 0; m < D; m++)
	    bmat(m+D*j, D*k+j) = grad(k, m);
      
      // gradient in reference coordinates:
      Mat<D,D> gradref = Trans (mlin) * gradr; 
      
      // gradient in real coordinates:
      Mat<D> grad1 = gradref * mip.GetJacobianInverse();
      
      // deformation gradient:
      Mat<D> matf = grad1 + Id<D>();
      Mat<D> strain = 0.5 * (Trans ( matf) * matf - Id<D>());
      Mat<D> piola2 = matlaw.CalcStress (mip, strain);
      
      
      for (int i1=0; i1 < D; i1++)
	for (int i2=0; i2 < D; i2++)
	  {
	    Mat<D> df = 0.0;
	    df(i1,i2) = 1;

	    Mat<D> strain = 0.5 * (Trans ( matf) * matf - Id<D>());
	    Mat<D> dstrain = 0.5 * ( Trans (df) * matf + Trans (matf) * df);

	    Mat<D> dpiola2 = matlaw.CalcDStress (mip, strain, dstrain);
	    Mat<D> dpiola1 = matf * dpiola2 + df * piola2;

	    for (int j1=0; j1 < D; j1++)
	      for (int j2=0; j2 < D; j2++)
		dmat(i1*D+i2,j1*D+j2) = dpiola1(j1,j2);
	  }


      dmat *= mip.GetWeight();

      dbmat = dmat * bmat;
      elmat += Trans (dbmat) * bmat;
    }
}



template <int D>
void NonlinMechIntegrator<D> ::
ApplyElementMatrix (const FiniteElement & bfel,
		    const ElementTransformation & eltrans,
		    const FlatVector<double> & elx,
		    FlatVector<double> & ely,
		    void * precomputed,
		    LocalHeap & lh) const
{
  static Timer timer ("NonLinMech::Apply");
  RegionTimer reg (timer);

  const ScalarFiniteElement<D> & fel =
    dynamic_cast<const ScalarFiniteElement<D>&> (bfel);
  int ndof = fel.GetNDof ();

  ely = 0;
  FlatMatrixFixWidth<D> melx(ndof, &elx[0]);
  FlatMatrixFixWidth<D> mely(ndof, &ely[0]);
  FlatMatrixFixWidth<D> dshape(ndof, lh);

  int order = 2 * fel.Order();
  const IntegrationRule & ir = SelectIntegrationRule (fel.ElementType(), order);

  for (int i = 0; i < ir.GetNIP(); i++)
    {
      MappedIntegrationPoint<D,D> mip (ir[i], eltrans);

      fel.CalcDShape (ir[i], dshape);
      Mat<D> gradref = Trans (melx) * dshape;
      Mat<D> grad = gradref * mip.GetJacobianInverse();
      Mat<D> matf = grad + Id<D>();
      Mat<D> strain = 0.5 * (Trans(matf) * matf - Id<D>());

      Mat<D> piola2 = matlaw.CalcStress (mip, strain);
      Mat<D> piola1 = matf * piola2;
      double fac = mip.GetWeight();

      Mat<D> trans = mip.GetJacobianInverse() * Trans (piola1);
      mely += fac * dshape * trans;
    }
}


template <int D>
double NonlinMechIntegrator<D> ::
Energy (const FiniteElement & bfel,
	const ElementTransformation & eltrans,
	const FlatVector<double> & elx,
	LocalHeap & lh) const
{
  static Timer timer ("NonLinMech::Energy");
  RegionTimer reg (timer);
  
  const ScalarFiniteElement<D> & fel =
    dynamic_cast<const ScalarFiniteElement<D>&> (bfel);
  int ndof = fel.GetNDof ();

  FlatMatrixFixWidth<D> melx(ndof, &elx(0));
  FlatMatrixFixWidth<D> dshape(ndof, lh);

  int order = 2 * fel.Order();
  const IntegrationRule & ir = SelectIntegrationRule (fel.ElementType(), order);

  double energy = 0;

  for (int i = 0; i < ir.GetNIP(); i++)
    {
      MappedIntegrationPoint<D,D> mip (ir[i], eltrans);
      fel.CalcDShape (ir[i], dshape);
      Mat<D> gradref = Trans (melx) * dshape;
      Mat<D> grad = gradref * mip.GetJacobianInverse();
      Mat<D> matf = grad + Id<D>();
      Mat<D> strain = 0.5 * (Trans ( matf) * matf - Id<D>());

      energy += mip.GetWeight() * matlaw.Energy (mip, strain);
    }

  return energy;
}


template <int D>
void NonlinMechIntegrator<D> ::
CalcFlux (const FiniteElement & bfel,
	  const BaseMappedIntegrationPoint & bmip,
	  const FlatVector<double> & elx,
	  FlatVector<double> & flux,
	  bool applyd,
	  LocalHeap & lh) const
{
  const ScalarFiniteElement<D> & fel =
    dynamic_cast<const ScalarFiniteElement<D>&> (bfel);

  int ndof = fel.GetNDof ();

  FlatMatrix<double> melx(ndof, D, &elx(0));

  const MappedIntegrationPoint<D,D> & mip = 
    static_cast<const MappedIntegrationPoint<D,D> & > (bmip);

  FlatMatrixFixWidth<D> dshape(ndof, lh);
  fel.CalcDShape (mip.IP(), dshape);

  Mat<D> gradref = Trans (melx) * dshape;
  Mat<D> grad = gradref * mip.GetJacobianInverse();
  Mat<D> matf = grad + Id<D>();
  Mat<D> strain = 0.5 * (Trans ( matf) * matf - Id<D>());

  Mat<D> piola2 = matlaw.CalcStress (mip, strain);
  Mat<D> cauchy = matf * piola2 * Trans (matf);

  Mat<D> mflux = applyd ? cauchy : strain;

  if (D == 2)
    {
      flux(0) = mflux(0,0);
      flux(1) = mflux(1,1);
      flux(2) = mflux(0,1);
    }
  else
    {
      flux(0) = mflux(0,0);
      flux(1) = mflux(1,1);
      flux(2) = mflux(2,2);
      flux(3) = mflux(1,2);
      flux(4) = mflux(0,2);
      flux(5) = mflux(0,1);
    }
}



static RegisterBilinearFormIntegrator<NonlinMechIntegrator<2> > initnlm2 ("nonlinmech", 2, 2);
static RegisterBilinearFormIntegrator<NonlinMechIntegrator<3> > initnlm3 ("nonlinmech", 3, 2);


