/*
  We construct a high order H1 finite element space
 */

#include <comp.hpp>
using namespace ngcomp;
using ngfem::INT;


class PeriodicH1HighOrderFESpace : public H1HighOrderFESpace
{
  Array<int> dofmap; // mapping of dofs
  int idnr;   // number of periodic identification
public:
  PeriodicH1HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : H1HighOrderFESpace (ama, flags)
  {
    idnr = int(flags.GetNumFlag ("idnr", 1));
    cout << " initialize periodic fespace" << endl;
  }

  virtual ~PeriodicH1HighOrderFESpace () { ; }

  virtual void Update (LocalHeap & lh)
  {
    H1HighOrderFESpace::Update (lh);
    
    dofmap.SetSize (GetNDof());
    for (int i = 0; i < dofmap.Size(); i++)
      dofmap[i] = i;

    
    Array<INT<2> > per_verts; 
    ma->GetPeriodicVertices(idnr, per_verts);

    // first dofs are vertex dofs
    for (auto pair : per_verts)
      dofmap[pair[1]] = pair[0];


    // find periodic edges
    // build vertex-pair to edge hashtable:
    HashTable<INT<2>, int> vp2e(ma->GetNEdges());

    for (int enr = 0; enr < ma->GetNEdges(); enr++)
      {
        int v1, v2;
        ma->GetEdgePNums (enr, v1, v2);
        if (v1 > v2) Swap (v1, v2);
        vp2e[INT<2>(v1,v2)] = enr;
      }

    for (int enr = 0; enr < ma->GetNEdges(); enr++)
      {
        int v1, v2;
        ma->GetEdgePNums (enr, v1, v2);
        // number of master-vertices
        // use that dofmap[0:nv] is exactly the vertex-map
        int mv1 = dofmap[v1];   // 
        int mv2 = dofmap[v2];
        if (v1 != mv1 && v2 != mv2) // edge shall be mapped
          {
            if (mv1 > mv2) Swap (mv1, mv2);
            int menr = vp2e[INT<2>(mv1,mv2)];  // the master edge-nr

            IntRange edofs = GetEdgeDofs (enr);   // dofs on slave edge
            IntRange medofs = GetEdgeDofs (menr); // dofs on master edge
            for (int i = 0; i < edofs.Size(); i++)
              dofmap[edofs[i]] = medofs[i];
          }
      }


    HashTable<INT<3>, int> v2f(ma->GetNFaces());

    Array<int> pnums;
    for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
      {
        ma->GetFacePNums (fnr, pnums);
        INT<3> i3(pnums[0], pnums[1], pnums[2]);
        i3.Sort();
        v2f[i3] = fnr;
      }

    for (int fnr = 0; fnr < ma->GetNFaces(); fnr++)
      {
        ma->GetFacePNums (fnr, pnums);
        INT<3> i3(dofmap[pnums[0]], dofmap[pnums[1]], dofmap[pnums[2]]);
        if (i3[0] != pnums[0] && i3[1] != pnums[1] && i3[2] != pnums[2])
          {
            i3.Sort();
            int mfnr = v2f[i3];

            IntRange fdofs = GetFaceDofs (fnr);
            IntRange mfdofs = GetFaceDofs (mfnr);
            for (int i = 0; i < fdofs.Size(); i++)
              dofmap[fdofs[i]] = mfdofs[i];
          }
      }

    
    
    for (int i = 0; i < dofmap.Size(); i++)
      if (dofmap[i] != i)
        ctofdof[i] = UNUSED_DOF;
  }
  
  virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const
  {
    Ngs_Element ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TRIG:
        {
          auto fe = new(alloc) H1HighOrderFE<ET_TRIG> (order);
          fe->SetVertexNumbers (dofmap[ngel.Vertices()]);
          return *fe;
        }
      case ET_TET:
        {
          auto fe = new(alloc) H1HighOrderFE<ET_TET> (order);
          fe->SetVertexNumbers (dofmap[ngel.Vertices()]);
          return *fe;
        }
      default:
        throw Exception (string ("PeriodicH1, element not implemented: ")+
                         ElementTopology::GetElementName (ngel.GetType()));
      }
    throw Exception("GetFE: undefined element");
  }
  

  virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    H1HighOrderFESpace::GetDofNrs (ei, dnums);
    for (size_t i = 0; i < dnums.Size(); i++)
      dnums[i] = dofmap[dnums[i]];
  }
};

static RegisterFESpace<PeriodicH1HighOrderFESpace> myinitifes ("periodic_h1ho");
