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


#include "myElement.hpp"
#include "myFESpace.hpp"

/*
#include <diffop_impl.hpp>
#ifdef WIN32
      template ngcomp::T_DifferentialOperator<ngcomp::DiffOpId<2> >;
#endif
*/

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

    
    // needed to draw solution function
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(), 
                                                 make_shared<ConstantCoefficientFunction>(1));
  }
    
  
  MyFESpace :: ~MyFESpace ()
  {
    // nothing to do
  }

  
  void MyFESpace :: Update(LocalHeap & lh)
  {
    // some global update:

    cout << "Update MyFESpace, #vert = " << ma->GetNV() 
         << ", #edge = " << ma->GetNEdges() << endl;

    nvert = ma->GetNV();
    // number of dofs:
    if (!secondorder)
      ndof = nvert;  // number of vertices
    else
      ndof = nvert + ma->GetNEdges();  // num vertics + num edges
  }

  void MyFESpace :: GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment
    
    dnums.SetSize(0);

    Array<int> vert_nums;
    ma->GetElVertices (ei, vert_nums);  // global vertex numbers

    // first dofs are vertex numbers:
    for (auto v : vert_nums)
      dnums.Append (v);

    if (secondorder)
      {
        // more dofs on edges:

        Array<int> edge_nums;
        ma->GetElEdges (ei, edge_nums);    // global edge numbers

        for (auto e : edge_nums)
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



  /*
    register fe-spaces
    Object of type MyFESpace can be defined in the pde-file via
    "define fespace v -type=myfespace"
  */

  static RegisterFESpace<MyFESpace> initifes ("myfespace");
}
