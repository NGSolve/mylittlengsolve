/*********************************************************************/
/* File:   plasticity_simple.cpp                                     */
/* Author: AFEM09, RWTH Aachen                                       */
/* Date:   9. May. 2009                                              */
/*********************************************************************/


/*

A simple solver for a model problem for elasto-plasticity

yield function is

f(sigma) = | sigma | - sigma_y

by dualization we obtain the minimization problem in the
displacement u and the plastic strain p

min \int  1/2 | eps(u) - p |_D^2  + | p |  - b u  dx

*/

#include <solve.hpp>

using namespace ngsolve;

namespace plasticity_simple
{
  
  class PlasticityFESpace : public CompoundFESpace
  {
  public:
    PlasticityFESpace (const MeshAccess & ama, 		   
                       const Array<const FESpace*> & aspaces, 
                       const Flags & flags)
      : CompoundFESpace (ama, aspaces, flags)
    { ; }
    
    virtual ~PlasticityFESpace () { ; }
    
    static FESpace * Create (const MeshAccess & ma, const Flags & flags)
    {
      // space has 2 displacement and 3 strain components
      Array<const FESpace*> spaces(2+3);

      int order = int (flags.GetNumFlag ("order", 2));
      
      Flags uflags;
      uflags.SetFlag ("order", order);
      spaces[0] = new H1HighOrderFESpace (ma, uflags);
      spaces[1] = new H1HighOrderFESpace (ma, uflags);

      Flags pflags;
      pflags.SetFlag ("order", order-1);
      spaces[2] = new L2HighOrderFESpace (ma, pflags);
      spaces[3] = new L2HighOrderFESpace (ma, pflags);
      spaces[4] = new L2HighOrderFESpace (ma, pflags);

      PlasticityFESpace * fes = new PlasticityFESpace (ma, spaces, flags);
      return fes;
    }
  
    virtual string GetClassName () const { return "Displacement-Strain"; }
  };


  




  class PlasticityIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction *coef_E, *coef_nu, *coef_sigma_yield;
  public:
    
    PlasticityIntegrator (Array<CoefficientFunction*> & coefs)
    {
      coef_E = coefs[0];
      coef_nu = coefs[1];
      coef_sigma_yield = coefs[2];
    }
    
    static Integrator * Create (Array<CoefficientFunction*> & coefs)
    {
      return new PlasticityIntegrator (coefs);
    }
  
    virtual string Name () const { return "Plasticity"; }

    virtual bool BoundaryForm () const { return 0; }
    
    virtual void
    AssembleElementMatrix (const FiniteElement & fel,
                           const ElementTransformation & eltrans, 
                           FlatMatrix<double> & elmat,
                           LocalHeap & lh) const
    {
      FlatVector<> elveclin(fel.GetNDof(), lh);
      elveclin = 0;
      AssembleLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
    }
    

    double PhiStar (const Vec<3> & p) const
    {
      double reg = 1e-3;
      double reg2 = 0.1;

      double norm = sqrt (sqr(p(0))+sqr(p(1))+2*sqr(p(2)));  // norm of symmetric tensor
      if (norm > reg)
        return norm + (norm-1)*(norm-1)*reg2;
      else
        return 0.5 + norm*norm/reg + 0.5 * reg;
    }


    double PhiStarPrime (const Vec<3> & p, Vec<3> & prime) const
    {
      double reg = 1e-3;
      double reg2 = 0.1;
      double norm = sqrt (sqr(p(0))+sqr(p(1))+2*sqr(p(2)));

      double fac = 0;
      if (norm > reg)
        fac = 1.0/norm + 2*reg2*(norm-1) / norm;
      else
        fac = 2/reg;

      prime(0) = p(0) * fac;
      prime(1) = p(1) * fac;
      prime(2) = 2*p(2) * fac;
    }



    
    // compute the functional
    virtual double Energy (const FiniteElement & bfel, 
                           const ElementTransformation & eltrans, 
                           const FlatVector<double> & elx, 
                           LocalHeap & lh) const
    {
      const CompoundFiniteElement & cfel = dynamic_cast<const CompoundFiniteElement&> (bfel);

      const ScalarFiniteElement<2> & felu = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
      const ScalarFiniteElement<2> & felp = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[2]);

      int ndofu = felu.GetNDof();
      int ndofp = felp.GetNDof();
      
      FlatMatrix<> dshapeu(ndofu, 2, lh);
      FlatVector<> shapep(ndofp, lh);

      const IntegrationRule & ir = 
        SelectIntegrationRule (felu.ElementType(), 2*felu.Order());
      
      double energy = 0;
      
      for (int i = 0 ; i < ir.GetNIP(); i++)
        {
          HeapReset hr(lh);
          
          SpecificIntegrationPoint<2,2> sip(ir[i], eltrans, lh); 
          
          felu.CalcMappedDShape (sip, dshapeu);
          felp.CalcShape (ir[i], shapep);          

          Mat<2> gradu = 0.0;
          for (int j = 0; j < ndofu; j++)
            for (int k = 0; k < 2; k++)
              for (int l = 0; l < 2; l++)
                gradu(k,l) += dshapeu(j, l) * elx(k*ndofu+j);
          
          Mat<2> strain;
          strain = 0.5 * (gradu + Trans(gradu));
          
          Vec<3> p = 0.0;
          for (int j = 0; j < ndofp; j++)
            for (int k = 0; k < 3; k++)
              p(k) += shapep(j) * elx(2*ndofu+k*ndofp+j);
          
          strain(0,0) -= p(0);
          strain(1,1) -= p(1);
          strain(0,1) -= p(2);
          strain(1,0) -= p(2);

          double valE = coef_E -> Evaluate (sip);
          double fac = 0.5 * fabs (sip.GetJacobiDet()) * sip.IP().Weight();
          
          double sum = 0;
          for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
              sum += sqr (strain(k,l));

          energy += fac * valE * sum;
          energy += fac * PhiStar (p);
        }
      return energy;
    }
    
    
    // compute the gradient at a given point
    virtual void 
    ApplyElementMatrix (const FiniteElement & bfel, 
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> & elx,    // element vector
                        FlatVector<double> & ely,          // element gradient
                        void * precomputed,
                        LocalHeap & lh) const
    {
      /*
      double eps = 1e-6;

      int ndof = bfel.GetNDof();
      FlatVector<> xr(ndof, lh), xl(ndof, lh);

      for (int i = 0; i < ndof; i++)
        {
          xr = elx;
          xl = elx;
          xr(i) += eps;
          xl(i) -= eps;
          double energyr = Energy(bfel, eltrans, xr, lh);
          double energyl = Energy(bfel, eltrans, xl, lh);

          ely(i) = (energyr-energyl) / (2*eps);
        }
      */


      const CompoundFiniteElement & cfel = dynamic_cast<const CompoundFiniteElement&> (bfel);

      const ScalarFiniteElement<2> & felu = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
      const ScalarFiniteElement<2> & felp = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[2]);

      int ndofu = felu.GetNDof();
      int ndofp = felp.GetNDof();
      
      FlatMatrix<> dshapeu(ndofu, 2, lh);
      FlatVector<> shapep(ndofp, lh);

      const IntegrationRule & ir = 
        SelectIntegrationRule (felu.ElementType(), 2*felu.Order());
      
      ely = 0;

      for (int i = 0 ; i < ir.GetNIP(); i++)
        {
          HeapReset hr(lh);
          
          SpecificIntegrationPoint<2,2> sip(ir[i], eltrans, lh); 
          
          felu.CalcMappedDShape (sip, dshapeu);
          felp.CalcShape (ir[i], shapep);          

          Mat<2> gradu = 0.0;
          for (int j = 0; j < ndofu; j++)
            for (int k = 0; k < 2; k++)
              for (int l = 0; l < 2; l++)
                gradu(k,l) += dshapeu(j, l) * elx(k*ndofu+j);
          
          Mat<2> strain;
          strain = 0.5 * (gradu + Trans(gradu));
          
          Vec<3> p = 0.0;
          for (int j = 0; j < ndofp; j++)
            for (int k = 0; k < 3; k++)
              p(k) += shapep(j) * elx(2*ndofu+k*ndofp+j);
          
          strain(0,0) -= p(0);
          strain(1,1) -= p(1);
          strain(0,1) -= p(2);
          strain(1,0) -= p(2);

          double valE = coef_E -> Evaluate (sip);
          double fac = fabs (sip.GetJacobiDet()) * sip.IP().Weight();

          Mat<2> stress = valE * strain;


          Vec<3> dp;
          PhiStarPrime (p, dp);

          dp(0) -= stress(0,0);
          dp(1) -= stress(1,1);
          dp(2) -= 2 * stress(0,1);

          dp *= fac;
          stress *= fac;

          for (int j = 0; j < ndofu; j++)
            for (int k = 0; k < 2; k++)
              for (int l = 0; l < 2; l++)
                ely(k*ndofu+j) += stress(k,l) * dshapeu(j, l);

          for (int j = 0; j < ndofp; j++)
            for (int k = 0; k < 3; k++)
              ely(2*ndofu+k*ndofp+j) += dp(k) * shapep(j);
        }
    }

    
    
    
    // compute the Hesse Matrix at point elveclin
    virtual void
    AssembleLinearizedElementMatrix (const FiniteElement & bfel,
                                     const ElementTransformation & eltrans,
                                     FlatVector<double> & elveclin,
                                     FlatMatrix<double> & elmat,
                                     LocalHeap & lh) const
    {
      // numerical differentiation

      double eps = 1e-6;

      int ndof = bfel.GetNDof();

      // elmat.AssignMemory (ndof, ndof, lh);
      elmat = 0;

      FlatVector<> xr(ndof, lh), xl(ndof, lh);
      FlatVector<> yr(ndof, lh), yl(ndof, lh);

      for (int i = 0; i < ndof; i++)
        {
          xr = elveclin;
          xl = elveclin;
          xr(i) += eps;
          xl(i) -= eps;
          ApplyElementMatrix (bfel, eltrans, xr, yr, 0, lh);
          ApplyElementMatrix (bfel, eltrans, xl, yl, 0, lh);

          elmat.Col(i) = (0.5/eps) * (yr-yl);
        }


      /*
      const NodalFiniteElement<2> & fel = static_cast<const NodalFiniteElement<2>&> (bfel);
      int ndof = fel.GetNDof();

      elmat.AssignMemory (ndof, ndof, lh);
      elmat = 0;
    

      FlatVector<> shape(ndof, lh);
      const IntegrationRule & ir = 
        SelectIntegrationRule (fel.ElementType(), 2*fel.Order());


      for (int i = 0 ; i < ir.GetNIP(); i++)
        {
          HeapReset hr(lh);
          SpecificIntegrationPoint<2,2> sip(ir[i], eltrans, lh); 

          fel.CalcShape (ir[i], shape);

          double uilin = InnerProduct(shape, elveclin);
          double phiprimeprime = 12 * pow(uilin,2);

          double fac = fabs (sip.GetJacobiDet()) * sip.IP().Weight();

          elmat += (fac * phiprimeprime) * shape * Trans(shape);
        }
      */
    }


    virtual int DimElement () const { return 2; }
    virtual int DimSpace () const { return 2; }
    
    virtual int DimFlux () const { return 3; }

  virtual void
  CalcFlux (const FiniteElement & bfel,
            const BaseSpecificIntegrationPoint & bsip,
            const FlatVector<double> & elx, 
            FlatVector<double> & flux,
            bool applyd,
            LocalHeap & lh) const
  {
    flux.AssignMemory (3, lh);

    HeapReset hr(lh);

    
    const CompoundFiniteElement & cfel = dynamic_cast<const CompoundFiniteElement&> (bfel);
    
    const ScalarFiniteElement<2> & felu = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
    const ScalarFiniteElement<2> & felp = dynamic_cast<const ScalarFiniteElement<2>&> (cfel[2]);
    
    int ndofu = felu.GetNDof();
    int ndofp = felp.GetNDof();
    
    FlatVector<> shapeu(ndofu, lh);
    FlatMatrix<> dshapeu(ndofu, 2, lh);
    FlatVector<> shapep(ndofp, lh);
    
    
    const SpecificIntegrationPoint<2,2> & sip = static_cast<const SpecificIntegrationPoint<2,2> &> (bsip);
          
    felu.CalcShape (sip.IP(), shapeu);
    felu.CalcMappedDShape (sip, dshapeu);
    felp.CalcShape (sip.IP(), shapep);          
    
    Mat<2> gradu = 0.0;
    for (int j = 0; j < ndofu; j++)
      for (int k = 0; k < 2; k++)
        for (int l = 0; l < 2; l++)
          gradu(k,l) += dshapeu(j, l) * elx(k*ndofu+j);
    
    Mat<2> strain;
    strain = 0.5 * (gradu + Trans(gradu));
    
    Vec<3> p = 0.0;
    for (int j = 0; j < ndofp; j++)
      for (int k = 0; k < 3; k++)
        p(k) += shapep(j) * elx(2*ndofu+k*ndofp+j);
    
    strain(0,0) -= p(0);
    strain(1,1) -= p(1);
    strain(0,1) -= p(2);
    strain(1,0) -= p(2);
    
    double valE = coef_E -> Evaluate (sip);
    Mat<2> stress = valE * strain;

    flux(0) = stress(0,0);
    flux(1) = stress(1,1);
    flux(2) = stress(0,1);
  }
  
       
  };





  class NumProcNonlinearSolve : public NumProc
  {
  protected:
    BilinearForm * bfa;
    LinearForm * lff;
    GridFunction * gfu;

    int maxit;

  public:
    NumProcNonlinearSolve (PDE & apde, const Flags & flags)
      : NumProc (apde)
    { 
      bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
      lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));
      gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

      maxit = int ( flags.GetNumFlag ("maxit", 30));
    }

    virtual ~NumProcNonlinearSolve()
    { ; }

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcNonlinearSolve (pde, flags);
    }

    virtual void Do(LocalHeap & lh)
    {
      cout << "nonlinmag solver called" << endl;

      BaseVector & vecu = gfu->GetVector();
      const BaseVector & vecf = lff->GetVector();

      BaseVector & uold = *vecu.CreateVector();
      BaseVector & d = *vecu.CreateVector();
      BaseVector & w = *vecu.CreateVector();

      BilinearFormApplication applya(bfa);

      double err, errold, err0;
      double energy, energyold;

      d = vecf - applya * vecu;
      err0 = L2Norm(d);

      for (int i = 1; i <= maxit; i++)
        {
          cout << "newton it " << i << endl;
	
          bfa -> AssembleLinearization (vecu, lh);
          BaseMatrix & inva = 
            *dynamic_cast<BaseSparseMatrix&> (bfa -> GetMatrix()).InverseMatrix();


          d = vecf - applya * vecu;
          err = L2Norm(d);
          energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

          cout << " err = " << err/err0;
          cout << " energy = " << energy << endl;

          errold = err;
          energyold = energy;

          w = inva * d;
          uold = vecu;
          int lin_its = 0;
          double tau = 1;

          do
            {
              vecu = uold + tau * w;
              energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

              cout << "tau = " << tau
                   << " energy = " << energy << endl;

              tau *= 0.5;
            }
          while (energy > energyold && lin_its++ < 30 && err > 1e-7*err0);


          delete &inva;

          // if (err < 1e-7*err0) break;
        }
    }
  };











  /*
    Extract the (ux, uy) components
  */
  class DiffOpIdU : public DiffOp<DiffOpIdU>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 2 };
    enum { DIFFORDER = 0 };
    
    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      const CompoundFiniteElement & cfel = 
        dynamic_cast<const CompoundFiniteElement&> (bfel);
      const ScalarFiniteElement<2> & fel_u = 
        dynamic_cast<const ScalarFiniteElement<2>&> (cfel[0]);
      
      int nd_u = fel_u.GetNDof();
      
      FlatVector<> vecu(nd_u, lh);
      vecu = fel_u.GetShape (sip.IP(), lh);
      
      mat = 0;
      for (int j = 0; j < nd_u; j++)
        {
          mat(0,j) = vecu(j);
          mat(1,nd_u+j) = vecu(j);
        }
    }
  };
  
  
  class ElastoPlasticUIntegrator 
    : public T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
  {
  public:
    ///
    ElastoPlasticUIntegrator ()
      :  T_BDBIntegrator<DiffOpIdU, DiagDMat<2>, FiniteElement>
    (DiagDMat<2> (new ConstantCoefficientFunction(1)))
    { ; }
    
    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new ElastoPlasticUIntegrator ();
    }
    
    ///
    virtual string Name () const { return "ElastoPlastic IdU"; }
  };







  class Init
  { 
  public: 
    Init ();
  };
  
  Init::Init()
  {
    GetNumProcs().AddNumProc ("plasticitysimple", NumProcNonlinearSolve::Create);
    GetFESpaceClasses().AddFESpace ("plasticity", PlasticityFESpace::Create);
    GetIntegrators().AddBFIntegrator ("plasticity", 2, 3,
                                      PlasticityIntegrator::Create);

    GetIntegrators().AddBFIntegrator ("plasticity_idu", 2, 0,
                                      ElastoPlasticUIntegrator::Create);
  }
  
  Init init;
}
