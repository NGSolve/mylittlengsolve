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
  py::class_<MyFESpace, shared_ptr<MyFESpace>, FESpace>
    (m, "MyFESpace", "FESpace with first order and second order trigs on 2d mesh")
  /*
    for more complicated constructors, or if we have to polish some arguments coming from
    python before we can put them in the c++ constructor we can define the __init__ function.
    In the init function we need to construct the new object into the existing pointer.
   */
    .def("__init__", [] (MyFESpace * instance, shared_ptr<MeshAccess> mesh, py::object dirichlet,
                         py::object definedon, py::kwargs kwargs)
         {
          cout << "in init" << endl;
          auto flags = py::extract<Flags>(kwargs)();
          if (py::isinstance<py::list>(dirichlet)) {
            flags.SetFlag("dirichlet", makeCArray<double>(py::list(dirichlet)));
          }
          if (py::isinstance<py::str>(dirichlet))
            {
              cout << "dir is string" << endl;
              std::regex pattern(dirichlet.cast<string>());
              Array<double> dirlist;
              for (int i = 0; i < mesh->GetNBoundaries(); i++)
                if (std::regex_match (mesh->GetMaterial(BND, i), pattern))
                  dirlist.Append (i+1);
              cout << "dirlist = " << dirlist << endl;
              flags.SetFlag("dirichlet", dirlist);
            }

          if (py::isinstance<py::str>(definedon))
            {
              std::regex pattern(definedon.cast<string>());
              Array<double> defonlist;
              for (int i = 0; i < mesh->GetNDomains(); i++)
                if (regex_match(mesh->GetMaterial(VOL,i), pattern))
                  defonlist.Append(i+1);
              flags.SetFlag ("definedon", defonlist);
            }

          if (py::isinstance<py::list> (definedon))
            flags.SetFlag ("definedon", makeCArray<double> (definedon));
          py::extract<Region> definedon_reg(definedon);
          if (definedon_reg.check() && definedon_reg().IsVolume())
            {
              Array<double> defonlist;
              for (int i = 0; i < definedon_reg().Mask().Size(); i++)
                if (definedon_reg().Mask().Test(i))
                  defonlist.Append(i+1);
              flags.SetFlag ("definedon", defonlist);
            }

           auto fes = new (instance) MyFESpace(mesh,flags);
           LocalHeap lh(100000, "myfes - lh");
           fes->Update(lh);
           fes->FinalizeUpdate(lh);
         },
         /*
           pybind allows us to give the arguments keywords and default values. Secondorder is
           default false and the flags are defaulted to an empty dictionary. Python dicts
           are automatically converted to flags objects if passed. So you can always pass
           a flags object.
          */
         py::arg("mesh"), py::arg("dirichlet")=DummyArgument(), py::arg("definedon")=DummyArgument())
    ;
}

#endif // NGS_PYTHON
