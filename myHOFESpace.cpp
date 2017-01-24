/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

My own FESpace for high order finite elements

*/


#include <comp.hpp>    // provides FESpace, ...

#include "myHOElement.hpp"
#include "myHOFESpace.hpp"


namespace ngcomp
{

  MyHighOrderFESpace :: MyHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyHighOrderFESpace" << endl;

    order = int(flags.GetNumFlag ("order", 2));

    // needed to draw solution function
    /*
    evaluator = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();

    integrator = GetIntegrators().CreateBFI("mass", ma->GetDimension(), 
                                            make_shared<ConstantCoefficientFunction>(1));
    */
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
    
    integrator[VOL] = GetIntegrators().CreateBFI("mass", ma->GetDimension(), 
                                                 make_shared<ConstantCoefficientFunction>(1));
  }
    
  
  MyHighOrderFESpace :: ~MyHighOrderFESpace ()
  {
    // nothing to do
  }

  
  void MyHighOrderFESpace :: Update(LocalHeap & lh)
  {
    // some global update:

    int n_vert = ma->GetNV();  
    int n_edge = ma->GetNEdges(); 
    int n_cell = ma->GetNE();  

    first_edge_dof.SetSize (n_edge+1);
    int ii = n_vert;
    for (int i = 0; i < n_edge; i++, ii+=order-1)
      first_edge_dof[i] = ii;
    first_edge_dof[n_edge] = ii;
      
    first_cell_dof.SetSize (n_cell+1);
    for (int i = 0; i < n_cell; i++, ii+=(order-1)*(order-2)/2)
      first_cell_dof[i] = ii;
    first_cell_dof[n_cell] = ii;

    // cout << "first_edge_dof = " << endl << first_edge_dof << endl;
    // cout << "first_cell_dof = " << endl << first_cell_dof << endl;

    ndof = ii;
  }


  void MyHighOrderFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    // returns dofs of element number elnr

    dnums.SetSize(0);

    Ngs_Element ngel = ma->GetElement (ei);

    // vertex dofs
    for (auto v : ngel.Vertices())
      dnums.Append(v);

    // edge dofs
    for (auto e : ngel.Edges())
      {
        int first = first_edge_dof[e];
        int next  = first_edge_dof[e+1];
        for (int j = first; j < next; j++)
          dnums.Append (j);
      }

    if (ei.IsVolume())
      {
        int first = first_cell_dof[ei.Nr()];
        int next  = first_cell_dof[ei.Nr()+1];
        for (int j = first; j < next; j++)
          dnums.Append (j);
      }
  }

  
  FiniteElement & MyHighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TRIG:
        {
          MyHighOrderTrig * trig = new (alloc) MyHighOrderTrig(order);
          for (int i = 0; i < 3; i++)
            trig->SetVertexNumber (i, ngel.vertices[i]);
          return *trig;
        }
      case ET_SEGM:
        {
          MyHighOrderSegm * segm = new (alloc) MyHighOrderSegm(order);
          for (int i = 0; i < 2; i++)
            segm->SetVertexNumber (i, ngel.vertices[i]);
          return *segm;
        }
      default:
        throw Exception (string("Element type ")+ToString(ngel.GetType())+" not supported");
      }
  }
  



  static RegisterFESpace<MyHighOrderFESpace> initifes ("myhofespace");
}
