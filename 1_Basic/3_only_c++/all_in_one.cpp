/*********************************************************************/
/* File:   all_in_one.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   4. May. 2009                                              */
/*********************************************************************/

// ngscxx -c all_in_one.cpp; ngsld all_in_one.o -lngcomp -o all_in_one

/*

In this example we show how to use NGSolve C++ code only, without any interface.
We will solve a poisson problem here.
For this the mesh must be created by Netgen and stored in the file square.vol

*/

#include <comp.hpp>
using namespace ngcomp;

int main(int argc, char** argv)
{
  cout << "Start AllInOne" << endl;

  auto ma = make_shared<MeshAccess>("square.vol"); 
  Flags flags_fes;
  flags_fes.SetFlag ("order", 4);
  auto fes = make_shared<H1HighOrderFESpace> (ma, flags_fes);

  Flags flags_gfu;
  auto gfu = make_shared<T_GridFunction<double>> (fes, "u", flags_gfu);

  Flags flags_bfa;
  auto bfa = make_shared<T_BilinearFormSymmetric<double>> (fes, "a", flags_bfa);

  auto u = fes->GetTrialFunction();
  auto v = fes->GetTestFunction();

  bfa -> AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(u->Deriv() * v->Deriv(), VOL, VOL));

  // shared_ptr<BilinearFormIntegrator> bfi;
  // bfi = make_shared<LaplaceIntegrator<2>> (make_shared<ConstantCoefficientFunction> (1));
  // bfa -> AddIntegrator (bfi);

  Array<double> penalty(ma->GetNBoundaries());
  penalty = 0.0;
  penalty[0] = 1e10;

  // bfi = make_shared<RobinIntegrator<2>> (make_shared<DomainConstantCoefficientFunction> (penalty));
  // bfa -> AddIntegrator (bfi);
  auto pen = make_shared<DomainConstantCoefficientFunction> (penalty);
  bfa -> AddIntegrator(make_shared<SymbolicBilinearFormIntegrator>(pen*u*v, BND, VOL));

  Flags flags_lff;
  auto lff = make_shared<T_LinearForm<double>> (fes, "f", flags_lff);

  auto lfi = make_shared<SourceIntegrator<2>> (make_shared<ConstantCoefficientFunction> (5));
  lff -> AddIntegrator (lfi);

  LocalHeap lh(100000);
  fes -> Update();
  fes -> FinalizeUpdate();

  gfu -> Update();
  bfa -> Assemble(lh);
  lff -> Assemble(lh);

  const BaseMatrix & mata = bfa -> GetMatrix();
  const BaseVector & vecf = lff -> GetVector();
  BaseVector & vecu = gfu -> GetVector();

  auto inverse = mata.InverseMatrix(fes->GetFreeDofs());

  vecu = *inverse * vecf;

  cout << "<u,f> = " << InnerProduct(vecu, vecf) << endl;
  // cout << "Solution vector = " << endl << vecu << endl;
  return 0;
}
