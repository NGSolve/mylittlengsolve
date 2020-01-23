/*
  A tutorial for geometrically non-linear elasticity.

  The non-linear iteration
*/


#include <solve.hpp>
using namespace ngsolve;



class NumProcNonlinElast : public NumProc
{
protected:
  BilinearForm * bfa;
  LinearForm * lff;
  GridFunction * gfu;
  Preconditioner * pre;

  int maxsteps;
  double prec;

public:
  NumProcNonlinElast (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = pde.GetPreconditioner (flags.GetStringFlag ("preconditioner", ""), 1);
    maxsteps = int(flags.GetNumFlag ("maxsteps", 200));
    prec = flags.GetNumFlag ("prec", 1e-12);
  }
    
  virtual ~NumProcNonlinElast()  { ; }

  virtual void Do(LocalHeap & lh);

  virtual string GetClassName () const
  {
    return "NonlinElast";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Bilinear-form = " << bfa->GetName() << endl
	<< "Linear-form   = " << lff->GetName() << endl
	<< "Gridfunction  = " << gfu->GetName() << endl
	<< "Preconditioner = " << ((pre) ? pre->ClassName() : "None") << endl
	<< "precision     = " << prec << endl
	<< "maxsteps      = " << maxsteps << endl;
  }

  static void PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc NonlinElast:\n" \
      "------------\n" \
      "Solves the nonlinear system resulting from a boundary value problem\n\n" \
      "Required flags:\n" 
      "-bilinearform=<bfname>\n" 
      "    bilinear-form providing the matrix\n" \
      "-linearform=<lfname>\n" \
      "    linear-form providing the right hand side\n" \
      "-gridfunction=<gfname>\n" \
      "    grid-function to store the solution vector\n" 
      "\nOptional flags:\n"
      "-predoncitioner=<prename>\n"
      "-maxsteps=n\n"
      "-prec=eps\n"
	<< endl;
  }
};





void NumProcNonlinElast :: Do(LocalHeap & lh)
{
  cout << "Solve nonlinear elasticity" << endl;

  const FESpace & fes = bfa->GetFESpace();

  const BaseVector & vecf = lff->GetVector();
  BaseVector & vecu = gfu->GetVector();

  BaseVector & uold = *vecu.CreateVector();
  BaseVector & d = *vecu.CreateVector();
  BaseVector & w = *vecu.CreateVector();

  BilinearFormApplication applya(bfa);

  double err0;
    
  for (int i = 1; i <= maxsteps; i++)
    {
      bfa->AssembleLinearization (vecu, lh);
      bfa->GetMatrix().SetInverseType(SPARSECHOLESKY);
      BaseMatrix & inva = *bfa->GetMatrix().InverseMatrix(fes.GetFreeDofs());

      d = vecf - applya * vecu;
      w = inva * d;

      // err = L2Norm(d);
      double err = L2Norm(w);
      if (i == 1) err0 = err;

      double energyold = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

      cout << "newton it " << i;
      cout << " err = " << err/err0;
      cout << " energyold = " << energyold << endl;

      uold = vecu;

      int lin_its = 0;
      double tau = 1;

      do
	{
	  vecu = uold + tau * w;
	  double energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

	  cout << "tau = " << tau;
	  cout << " energy = " << energy << endl;

	  if (energy < energyold) break;

	  tau *= 0.5;
	}
      while (lin_its++ < 20);

      delete &inva;

      if (err < 1e-7*err0) break;
    }
}




static RegisterNumProc<NumProcNonlinElast> npinitbvp("nonlinelast");


