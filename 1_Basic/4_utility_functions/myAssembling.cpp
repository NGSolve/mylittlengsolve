/*********************************************************************/
/* File:   myAssembling.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   3. May. 2010                                              */
/*********************************************************************/


/*

  In this example we pass a fespace, a bilinearformintegrator and a linearformintegrator
  to our function. We will build the system matrix and vector in parallel together and
  solve the system. The result will be stored in a new created GridFunction which we will
  return. This function is exported to Python.

*/

#include <solve.hpp>
using namespace ngsolve;
#include "utility_functions.hpp"


namespace myassemble
{
  shared_ptr<GridFunction> MyAssemble(shared_ptr<FESpace> fes,
                                      shared_ptr<BilinearFormIntegrator> bfi,
                                      shared_ptr<LinearFormIntegrator> lfi)
  {
    cout << "We assemble matrix and rhs vector" << endl;

    auto ma = fes->GetMeshAccess();

    int ndof = fes->GetNDof();
    int ne = ma->GetNE();
    

    // setup element->dof table:

    // first we get the number of dofs per element ...
    Array<int> dnums;
    Array<int> cnt(ne);

    for (auto ei : ma->Elements(VOL))
      {
        fes->GetDofNrs (ei, dnums);
        cnt[ei.Nr()] = dnums.Size();
      }
      
    // allocate the table in compressed form ...
    Table<int> el2dof(cnt);

    // and fill it
    for (auto ei : ma->Elements(VOL))
      {
        fes->GetDofNrs (ei, dnums);
        el2dof[ei.Nr()] = dnums;
      }
    cout << "el2dof - table: " << el2dof << endl;

    // generate sparse matrix from element-to-dof table
    auto mat = make_shared<SparseMatrixSymmetric<double>> (ndof, el2dof);

    VVector<double> vecf (fes->GetNDof());

    *mat = 0.0;
    vecf = 0.0;

    /*
      Parallel iteration over elements using element coloring
      This automatically splits the lh into partial localheaps and
      calls the given lambda function with the element and the splitted lh
      We can create flatmatrices and -vectors on the local heap without
      memory management efficiently - the local heap will be reseted by
      IterateElements automatically.
    */
    LocalHeap lh(100000);
    IterateElements(*fes, VOL, lh, [&] (FESpace::Element el, LocalHeap &lh)
                    {
                      const ElementTransformation& eltrans = ma->GetTrafo(el,lh);
                      const FiniteElement& fel = fes->GetFE(el,lh);
                      auto dofs = el.GetDofs();
                      FlatMatrix<> elmat(dofs.Size(),lh);
                      bfi->CalcElementMatrix(fel,eltrans,elmat,lh);
                      mat->AddElementMatrixSymmetric(dofs,elmat);
                      FlatVector<> elvec(dofs.Size(),lh);
                      lfi->CalcElementVector(fel, eltrans, elvec, lh);
                      vecf.AddIndirect(dofs,elvec);
                    });

    *testout << "mat = " << *mat << endl;
    *testout << "vecf = " << vecf << endl;

    shared_ptr<BaseMatrix> inv = mat->InverseMatrix (fes->GetFreeDofs());

    Flags gfuFlags;
    auto gfu = make_shared<T_GridFunction<double>> (fes, "u", gfuFlags);
    gfu->Update();
    gfu -> GetVector() = (*inv) * vecf;
    return gfu;
  }
}
