/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

My own FESpace for linear and quadratic triangular elements.  

A fe-space provides the connection between the local reference
element, and the global mesh.

*/


#include <comp.hpp>    // provides FESpace, ...
#include <h1lofe.hpp>
#include <regex>

#include "myElement.hpp"
#include "myFESpace.hpp"
#include "myDiffOp.hpp"


namespace ngcomp
{

  MyFESpace :: MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyFESpace" << endl;
    cout << "Flags = " << flags << endl;

    secondorder = flags.GetDefineFlag ("secondorder");

    if (!secondorder)
      cout << "You have chosen first order elements" << endl;
    else
      cout << "You have chosen second order elements" << endl;

    
    // needed for symbolic integrators and to draw solution
    evaluator[VOL] = make_shared<T_DifferentialOperator<MyIdentity>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    // (still) needed to draw solution
    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(), 
                                                 make_shared<ConstantCoefficientFunction>(1));
  }
    
  
  
  void MyFESpace :: Update(LocalHeap & lh)
  {
    // some global update:
    cout << "Update MyFESpace, #vert = " << ma->GetNV() 
         << ", #edge = " << ma->GetNEdges() << endl;

    // number of vertices
    nvert = ma->GetNV();
    
    // number of dofs:
    ndof = nvert;
    if (secondorder)
      ndof += ma->GetNEdges();  // num vertics + num edges
  }

  void MyFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment
    
    dnums.SetSize(0);

    // first dofs are vertex numbers:
    for (auto v : ma->GetElVertices(ei))
      dnums.Append (v);

    if (secondorder)
      {
        // more dofs on edges:
        for (auto e : ma->GetElEdges(ei))
          dnums.Append (nvert+e);
      }
  }
  
  FiniteElement & MyFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if (ei.IsVolume())
      {
        if (!secondorder)
          return * new (alloc) MyLinearTrig;
        else
          return * new (alloc) MyQuadraticTrig;          
      }
    else
      {
        if (!secondorder)
          return * new (alloc) FE_Segm1;
        else
          return * new (alloc) FE_Segm2;          
      }
  }


  SymbolTable<shared_ptr<DifferentialOperator>>
  MyFESpace :: GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> additional;
    additional.Set ("defaultId", make_shared<T_DifferentialOperator<DiffOpId<2>>>());
    return additional;
  }


  /*
    register fe-spaces
    Object of type MyFESpace can be defined in the pde-file via
    "define fespace v -type=myfespace"
  */

  static RegisterFESpace<MyFESpace> initifes ("myfespace");
}

#ifdef NGS_PYTHON

void ExportMyFESpace(py::module m)
{
  using namespace ngcomp;
  /*
    We just export the class here and use the FESpace constructor to create our space.
    This has the advantage, that we do not need to specify all the flags to parse (like
    dirichlet, definedon,...), but we can still append new functions only for that space.
   */
  py::class_<MyFESpace, shared_ptr<MyFESpace>, FESpace>
    (m, "MyFESpace", "FESpace with first order and second order trigs on 2d mesh")
    .def("GetNVert", &MyFESpace::GetNVert)
    ;
}

#endif // NGS_PYTHON
