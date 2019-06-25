#ifndef FILE_MYHOFESPACE_HPP
#define FILE_MYHOFESPACE_HPP

/*********************************************************************/
/* File:   myHOFESpace.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

  My own FESpace for high order finite elements ...

*/


namespace ngcomp
{
  class MyHighOrderFESpace : public FESpace
  {
    int order;
    int ndof;
    
    Array<int> first_edge_dof;
    Array<int> first_cell_dof;
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    MyHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "MyHighOrderFESpace";
    }

    virtual void Update(LocalHeap & lh);
    virtual size_t GetNDof () const { return ndof; }

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const;
  };
}    

#endif
