/*********************************************************************/
/* File:   myPreconditioner.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

  We define a preconditioner ...

*/


#include <solve.hpp>    // provides FESpace, ...
#include <python_ngstd.hpp>
#include "myPreconditioner.hpp"


namespace ngcomp
{
  MyPreconditioner :: MyPreconditioner (shared_ptr<BilinearForm> abfa, Flags& flags)
    : Preconditioner(abfa, flags, "MyPreconditioner"), bfa(abfa)
  { ; }

  
  void MyPreconditioner :: Update()
  {
    const BaseSparseMatrix & mat = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
    shared_ptr<BitArray> freedofs = bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal());

    jacobi = mat.CreateJacobiPrecond (freedofs);

    if (test) Test();
  }
}

void ExportMyPreconditioner(py::module m)
{
  using namespace ngcomp;
  py::class_<MyPreconditioner, shared_ptr<MyPreconditioner>, Preconditioner>
    (m, "MyPreconditioner")
    .def(py::init<shared_ptr<BilinearForm>, Flags&>());
}

