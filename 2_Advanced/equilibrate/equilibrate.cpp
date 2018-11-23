// ngscxx  -c equilibrate.cpp ; ngsld -shared equilibrate.o -lngcomp -lngfem -lngstd -o libequilibrate.so

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
    for (int i = 0; i < nedges; i++)
      if (!dirichlet_edge[i])
        for (int v : ma.GetEdgePNums(i))
          creatoredges.Add(v, i);

  Table<int> patch_edges = creatoredges.MoveTable();

  // cout << "patch_elements:" << endl << patch_elements << endl;
  // cout << "patch_edges:" << endl << patch_edges << endl;

  const SparseMatrix<double> & matrix = dynamic_cast<SparseMatrix<double> &>(bf.GetMatrix());

  for (int i = 0; i < nv; ++i)
    {
      HeapReset hr(lh);      
      Array<int> globaldofs, dofs;

      for (int el : patch_elements[i])
        {
          fespace_x->GetInnerDofNrs(el, dofs);
          globaldofs += dofs;
        }

      for (int edge :  patch_edges[i])
        {
          fespace_x->GetEdgeDofNrs(edge, dofs);
          globaldofs += dofs;
        }

      if (!dirichlet_vertex[i])
        globaldofs += numberdofs;


      // cout << "globaldofs1 = " << endl << globaldofs << endl;
      
      // remove external sigma dofs:
      Array<int> facetdofs;
      for (int el : patch_elements[i])
        {
          ElementId ei(VOL, el);
          const FiniteElement & felx = fespace_x->GetFE(ei, lh);
          auto & cfelx = dynamic_cast<const CompoundFiniteElement&>(felx);
          auto & hdivfel = dynamic_cast<const HDivFiniteElement<2>&> (cfelx[0]);
          
          Array<int> dnums(felx.GetNDof(), lh);
          fespace_x->GetDofNrs(el, dnums);
          FlatArray<int> dnums_hdiv = dnums.Range(cfelx.GetRange(0));
          
          auto edges = ma.GetElement(ElementId(VOL,el)).Edges();
          for (int k = 0; k < 3; k++)
            if (ma.GetEdgePNums(edges[k])[0]!=i && ma.GetEdgePNums(edges[k])[1]!=i)
              {
                hdivfel.GetFacetDofs(k, facetdofs);
                for (int d : facetdofs)
                  globaldofs.DeleteElement(globaldofs.Pos(dnums_hdiv[d]));
              }

        }
      
      // cout << "globaldofs2 = " << endl << globaldofs << endl;

      int ldofs = globaldofs.Size();
      QuickSort(globaldofs);

      Matrix<> localmata(ldofs);
      for (int j = 0; j < ldofs; ++j)
        for (int k = 0; k < ldofs; ++k)
          localmata(j, k) = matrix(globaldofs[j], globaldofs[k]);

      FlatVector<> localrhs(ldofs, lh);
      localrhs = 0.0;

      for (int el : patch_elements[i])
        {
          HeapReset hr(lh);
          ElementId ei(VOL, el);

          const ElementTransformation & eltrans = ma.GetTrafo(ei, lh);

          const FiniteElement & felx = fespace_x->GetFE(ei, lh);
          const CompoundFiniteElement & cfelx = dynamic_cast<const CompoundFiniteElement&>(felx);
          const ScalarFiniteElement<2> & l2fel = dynamic_cast<const ScalarFiniteElement<2> &> (cfelx[1]);

          Array<int> dnums(felx.GetNDof(), lh);
          fespace_x->GetDofNrs(el, dnums);

          // compute right hand side (hat is vertex basis-function):
          //    int_T hat (source + div flux) v
          // 
          // use integration by parts:
          //    int_T hat source v - flux grad(hat v) + int_boundary(T) ...
          
          FlatVector<> shape(l2fel.GetNDof(), lh), elf(l2fel.GetNDof(), lh);
          IntegrationRule ir(eltrans.GetElementType(), 2*l2fel.Order()+2);


          MappedIntegrationRule<2,2> mir(ir, eltrans, lh);
          FlatVector<> hat(ir.Size(), lh), ri(ir.Size(), lh);
          FlatMatrix<> dhat(ir.Size(), 2, lh);
          FlatMatrix<> dri(ir.Size(), 2, lh), fi(ir.Size(), 1, lh), fluxi(ir.Size(), 2, lh);

          // vertex hat-basis function and gradient
          ScalarFE<ET_TRIG, 1> p1fe;
          FlatVector<> p1vec(p1fe.GetNDof(), lh);
          p1vec = 0;
          p1vec(ma.GetElement(ei).Vertices().Pos(i)) = 1;
          
          p1fe.Evaluate(ir, p1vec, hat);
          p1fe.EvaluateGrad(ir, p1vec, dhat);

          source.Evaluate(mir, fi);
          flux.Evaluate(mir, fluxi);
          
          for (int k = 0; k < ir.GetNIP(); k++)
            {
              auto & mip = mir[k];
              Vec<2> gradhat = Trans(mip.GetJacobianInverse()) * dhat.Row(k);
              ri(k) = mip.GetWeight() * (hat(k)*fi(k,0) - InnerProduct(gradhat, fluxi.Row(k)));
              dri.Row(k) = -hat(k)*mip.GetWeight() * mip.GetJacobianInverse()*fluxi.Row(k);
            }

          l2fel.EvaluateTrans(ir, ri, elf);
          l2fel.EvaluateGradTrans(ir, dri, shape);
          elf += shape;


          // compute element boundary terms:
          // int_boundary(T)  hat * sigma_n * (v - vfacet)
          
          const FacetVolumeFiniteElement<2> & facetfel =
            dynamic_cast<const FacetVolumeFiniteElement<2> &> (cfelx[2]);
          
          ngfem::ELEMENT_TYPE eltype = cfelx.ElementType();
          Facet2ElementTrafo transform(eltype);

          // outward normals on ref-element
          FlatVector<Vec<2>> normals = ElementTopology::GetNormals<2>(eltype);
          int nfacet = ElementTopology::GetNFacets(eltype);
          
          FlatVector<> f_facet(facetfel.GetNDof(), lh);
          for (int k = 0; k < nfacet; k++)
            {
              FlatVector<> shape_facet(facetfel.Facet(k).GetNDof(), lh);
              
              ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType(eltype, k);
              IntegrationRule ir_facet(etfacet, facetfel.Order() * 2+2);

              FlatVector<> rfi(ir_facet.Size(), lh);

              IntegrationRule & ir_vol = transform(k, ir_facet, lh);
              MappedIntegrationRule<2,2> mir(ir_vol, eltrans, lh);
              FlatMatrix<> fluxi(mir.Size(), 2, lh);
              FlatVector<> hat(mir.Size(), lh);

              p1fe.Evaluate (ir_vol, p1vec, hat);
              flux.Evaluate (mir, fluxi);
              
              for (int l = 0; l < ir_facet.GetNIP(); l++)
                {
                  auto & mip = mir[l];
                  Vec<2> normal = mip.GetMeasure() * Trans(mip.GetJacobianInverse()) * normals[k];
                  rfi(l) = -InnerProduct (fluxi.Row(l), normal) * ir_facet[l].Weight() * hat(l); 
                }
              
              l2fel.EvaluateTrans (ir_vol, rfi, shape);
              elf -= shape;
              facetfel.Facet(k).EvaluateTrans(ir_vol, rfi, f_facet.Range(facetfel.GetFacetDofs(k)));
            }

          FlatArray<int> dnums_l2 = dnums.Range(cfelx.GetRange(1));
          for (int k = 0; k < dnums_l2.Size(); k++)
            localrhs[globaldofs.Pos(dnums_l2[k])] = elf[k];
          FlatArray<int> dnums_facet = dnums.Range(cfelx.GetRange(2));
          for (int k = 0; k < dnums_facet.Size(); k++)
            if (globaldofs.Pos(dnums_facet[k]) != -1)
              localrhs[globaldofs.Pos(dnums_facet[k])] += f_facet[k];
        }
      
      FlatVector<> sol(ldofs, lh);
      CalcInverse(localmata);
      sol = localmata * localrhs;

      IntRange sigmadofs = fespace_x->GetRange(0);
      int end_sigma_dofs = 0;
      for (int i = 0; i < globaldofs.Size(); i++)
        if (globaldofs[i] < sigmadofs.end())
          end_sigma_dofs = i+1;
      equflux.GetVector().AddIndirect(globaldofs.Range(0, end_sigma_dofs), sol.Range(0, end_sigma_dofs));
    }  
}

PYBIND11_MODULE(libequilibrate, m) {
  m.def("EquilibratePatches", &EquilibratePatches)
    ;
}
