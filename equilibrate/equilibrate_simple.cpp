// ngscxx -shared equilibrate_simple.cpp -lngcomp -o libequilibrate.so
#include <comp.hpp>
#include <python_ngstd.hpp>
using namespace ngcomp;

void EquilibratePatches (CoefficientFunction & flux,
                         CoefficientFunction & source,
                         BilinearForm & bf,
                         GridFunction & equflux,
                         FESpace & fespace_u)
{
  equflux.GetVector() = 0.0;
  MeshAccess & ma = *equflux.GetMeshAccess();
  shared_ptr<CompoundFESpace> fespace_x = dynamic_pointer_cast<CompoundFESpace> (bf.GetFESpace());

  LocalHeap lh(100000, "equilibration-lh");

  /*
  cout << "hdiv-dofs   = " << fespace_x->GetRange(0) << endl;
  cout << "L2-dofs     = " << fespace_x->GetRange(1) << endl;
  cout << "facet-dofs  = " << fespace_x->GetRange(2) << endl;
  cout << "number-dofs = " << fespace_x->GetRange(3) << endl;
  */
  
  IntRange numberdofs = fespace_x->GetRange(3);

  int nv = ma.GetNV();
  int ne = ma.GetNE();
  int nedges = ma.GetNEdges();

  // find vertices and edges on Dirichlet boundary
  BitArray dirichlet_vertex(nv);
  BitArray dirichlet_edge(nedges);
  dirichlet_vertex.Clear();
  dirichlet_edge.Clear();

  for (Ngs_Element el : ma.Elements(BND))
    if (fespace_u.IsDirichletBoundary (el.GetIndex()))
      {
        for (int v : el.Vertices())
          dirichlet_vertex.Set(v);
        for (int e : el.Edges())
          dirichlet_edge.Set(e);
      }

  cout << "diredge: " << endl << dirichlet_edge << endl;
  cout << "dirvert: " << endl << dirichlet_vertex << endl;

  // build vertex -> element table:
  TableCreator<int> creator(nv);
  for (; !creator.Done(); creator++)
    for (Ngs_Element el : ma.Elements(VOL))
      for (int v : el.Vertices())
        creator.Add (v, el.Nr());

  Table<int> patch_elements = creator.MoveTable();

  // build vertex -> edge table:  
  TableCreator<int> creatoredges(nv);
  for (; !creatoredges.Done(); creatoredges++)
    {
      for (int i = 0; i < ne; i++)
        {
          Array<int> enums, vnums, fvnums;
          vnums = ma.GetElVertices(i);
          ma.GetElEdges(i, enums);
          for (int e : enums)
            {
              ma.GetEdgePNums(e, fvnums);
              for (int v : vnums)
                if (!fvnums.Contains(v))
                  creatoredges.Add(v, e); // edge on patch boundary
            }
        }
      
      for (int i = 0; i < nedges; i++)
        if (!dirichlet_edge[i])
          {
            Array<int> vnums;
            ma.GetEdgePNums(i, vnums);
            for (int v : vnums)
              creatoredges.Add(v, i);
          }
    }
  
  Table<int> patch_edges = creatoredges.MoveTable();

  // cout << "patch_elements:" << endl << patch_elements << endl;
  // cout << "patch_edges:" << endl << patch_edges << endl;

  const SparseMatrix<double> & matrix = dynamic_cast<SparseMatrix<double> &>(bf.GetMatrix());

  for (int i = 0; i < nv; ++i)
    {
      Array<int> globaldofs, dofs;

      for (int j = 0; j < patch_elements[i].Size(); j++)
        {
          fespace_x->GetInnerDofNrs(patch_elements[i][j], dofs);
          globaldofs += dofs;
        }

      for (int j = 0; j < patch_edges[i].Size(); j++)
        {
          fespace_x->GetEdgeDofNrs(patch_edges[i][j], dofs);
          globaldofs += dofs;
        }

      if (!dirichlet_vertex[i])
        globaldofs += numberdofs;

      // cout << "globaldofs = " << endl << globaldofs << endl;

      int ldofs = globaldofs.Size();
      QuickSort(globaldofs);

      Matrix<> localmata(ldofs);
      for (int j = 0; j < ldofs; ++j)
        for (int k = 0; k < ldofs; ++k)
          localmata(j, k) = matrix(globaldofs[j], globaldofs[k]);

      Vector<> localrhs(ldofs);
      localrhs = 0.0;

      for (int j = 0; j < patch_elements[i].Size(); j++)
        {
          HeapReset hr(lh);
          int el = patch_elements[i][j];
          ElementId ei(VOL, el);

          const ElementTransformation & eltrans = ma.GetTrafo(ei, lh);

          const FiniteElement & felx = fespace_x->GetFE(ei, lh);
          const CompoundFiniteElement & cfelx = dynamic_cast<const CompoundFiniteElement&>(felx);
          const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelx[1]);

          int nd_vr = l2fel.GetNDof();
          Array<int> dnumss(felx.GetNDof(), lh);
          fespace_x->GetDofNrs(el, dnumss);

          FlatVector<> shape(nd_vr, lh), elf(nd_vr, lh);
          FlatMatrix<> dshape(nd_vr, 2, lh);

          IntegrationRule ir(eltrans.GetElementType(), 2*l2fel.Order()+2);

          ScalarFE<ET_TRIG, 1> p1fe;
          FlatVector<> p1vec(p1fe.GetNDof(), lh);
          Array<int> vnums;
          ma.GetElVertices(el, vnums);
          for (int k = 0; k < p1fe.GetNDof(); k++)
            p1vec(k) = (vnums[k] == i) ? 1 : 0;

          elf = 0.0;
          for (int k = 0; k < ir.GetNIP(); k++)
            {
              MappedIntegrationPoint<2, 2> mip(ir[k], eltrans);
              double hatfunction = p1fe.Evaluate(ir[k], p1vec);
              Vec<2> gradhatfunction = Trans(mip.GetJacobianInverse()) * p1fe.EvaluateGrad(ir[k], p1vec);
              Vec<1> fi;
              Vec<2> fluxi;
              source.Evaluate(mip, fi);
              flux.Evaluate(mip, fluxi);
              l2fel.CalcShape (ir[k], shape);
              l2fel.CalcMappedDShape (mip, dshape);
              elf += mip.GetWeight() * (hatfunction*fi(0) - InnerProduct(gradhatfunction, fluxi)) * shape;
              elf -= (hatfunction*mip.GetWeight()) * dshape * fluxi;
            }

          const FacetVolumeFiniteElement<2> & facetfel =
            dynamic_cast<const FacetVolumeFiniteElement<2> &> (cfelx[2]);
          
          ngfem::ELEMENT_TYPE eltype = cfelx.ElementType();
          Facet2ElementTrafo transform(eltype);

          // outward normals on ref-element
          FlatVector<Vec<2>> normals = ElementTopology::GetNormals<2>(eltype);
          int nfacet = ElementTopology::GetNFacets(eltype);
          
          int nd_facet = facetfel.GetNDof();
          FlatVector<> shape_facet(nd_facet, lh), f_facet(nd_facet, lh);
          f_facet = 0.0;

          for (int k = 0; k < nfacet; k++)
            {
              // outward normal of reference element facet
              Vec<2> normal_ref = normals[k];
              
              // integration rule on facet
              ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType(eltype, k);
              IntegrationRule ir_facet(etfacet, facetfel.Order() * 2+2);
              
              for (int l = 0; l < ir_facet.GetNIP(); l++)
                {
                  IntegrationPoint ip = transform(k, ir_facet[l]);   // point in ref element
                  MappedIntegrationPoint<2, 2> mip(ip, eltrans);     // point in physical element
                  
                  // normal of physical element
                  Vec<2> normal = fabs(mip.GetJacobiDet()) * Trans(mip.GetJacobianInverse()) * normal_ref;
                  double len = L2Norm(normal);
                  double weight = ir_facet[l].Weight() * len;
                  normal /= len;

                  shape_facet = 0.0;
                  facetfel.Facet(k).CalcShape(ip, shape_facet.Range(facetfel.GetFacetDofs(k)));
                  l2fel.CalcShape (ip, shape);
                  
                  Vec<2> fluxi;
                  flux.Evaluate(mip, fluxi);
                  double fac = InnerProduct (fluxi, normal) * weight * p1fe.Evaluate(ip, p1vec);
                  
                  f_facet -= fac * shape_facet;
                  elf += fac * shape;
                }
            }

          FlatArray<int> dnums_l2 = dnumss.Range(cfelx.GetRange(1));
          for (int k = 0; k < dnums_l2.Size(); k++)
            localrhs[globaldofs.Pos(dnums_l2[k])] = elf[k];
          FlatArray<int> dnums_facet = dnumss.Range(cfelx.GetRange(2));
          for (int k = 0; k < dnums_facet.Size(); k++)
            if (globaldofs.Pos(dnums_facet[k]) != -1)
              localrhs[globaldofs.Pos(dnums_facet[k])] += f_facet[k];
        }
      
      // *testout << "localmata = " << endl << localmata << endl;
      // *testout << "localrhs = " << endl << localrhs << endl;
      
      Vector<> sol(ldofs);
      Matrix<> inv = localmata;
      CalcInverse(inv);
      sol = inv * localrhs;

      // *testout << "sol = " << endl << sol << endl;
      
      IntRange sigmadofs = fespace_x->GetRange(0);
      int end_sigma_dofs = 0;
      for (int i = 0; i < globaldofs.Size(); i++)
        if (globaldofs[i] < sigmadofs.end())
          end_sigma_dofs = i+1;
      equflux.GetVector().AddIndirect(globaldofs.Range(0, end_sigma_dofs), sol.Range(0, end_sigma_dofs));
    }  
}


BOOST_PYTHON_MODULE(libequilibrate) {
  
  bp::def("EquilibratePatches", &EquilibratePatches);

}
