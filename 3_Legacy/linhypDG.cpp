
/*
  
Solver for the linear hyperbolic equation

du/dt  +  div (b u) = 0

by an explicit time-stepping method

*/



#include <solve.hpp>

using namespace ngsolve;
using ngfem::ELEMENT_TYPE;



template <int D>
class NumProcLinearHyperbolic : public NumProc
{
protected:
  shared_ptr<CoefficientFunction> cfflow;
  shared_ptr<GridFunction> gfu;

  double dt;
  double tend;

  Timer timer_element, timer_facet, timer_mass;

  class FacetData
  {
  public:
    int elnr[2];
    int facetnr[2];
    Vector<> flown;

    FacetData() = default;
    // FacetData(const FacetData & ed2) = default;
    // FacetData(FacetData && ed2) = default;
    FacetData (int nip) : flown(nip) { ; }
    // ~FacetData() = default;
    // FacetData & operator= (const FacetData & el2) = default;
    // FacetData & operator= (FacetData && el2) = default;
  };


  class ElementData
  {
  public:
    MatrixFixWidth<D> flowip;
    Matrix<> invmass;
    
    ElementData() = default;
    // ElementData(const ElementData & ed2) = default;
    // ElementData(ElementData && ed2) = default;
    ElementData (int ndof, int nip) : flowip(nip), invmass(ndof) { ; }
    // ~ElementData() = default;
    // ElementData & operator= (const ElementData & el2) = default;
    // ElementData & operator= (ElementData && el2) = default;
  };

  Array<FacetData> facetdata;         // new move semantics
  Array<ElementData> elementdata;     // new move semantics

public:
    
  NumProcLinearHyperbolic (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde), 
      timer_element("convection - time element"), 
      timer_facet("convection - time facet"),
      timer_mass("convection - time mass")
  {
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    cfflow = apde->GetCoefficientFunction (flags.GetStringFlag ("flow", "flow"));

    dt = flags.GetNumFlag ("dt", 0.001);
    tend = flags.GetNumFlag ("tend", 1);
  }



  virtual void Do(LocalHeap & lh)
  {
    cout << "solve conservation equation" << endl;



    // prepare ...

    auto & fes = dynamic_cast<const L2HighOrderFESpace&> (*gfu->GetFESpace());
    
    elementdata.SetAllocSize (ma->GetNE());

    MassIntegrator<D> bfi(make_shared<ConstantCoefficientFunction> (1));

    for (auto ei : ma->Elements())
      {
	HeapReset hr(lh);
	
	auto & fel = dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (ei, lh));
	const IntegrationRule ir(fel.ElementType(), 2*fel.Order());

        const_cast<DGFiniteElement<D>&> (fel).PrecomputeShapes (ir);
        const_cast<DGFiniteElement<D>&> (fel).PrecomputeTrace ();

        auto & trafo = ma->GetTrafo(ei, lh);
	MappedIntegrationRule<D,D> mir(ir, trafo, lh);
	
	ElementData edi(fel.GetNDof(), ir.Size());
        cfflow -> Evaluate (mir, edi.flowip);
			    
	for (int j = 0; j < ir.Size(); j++)
	  {
	    Vec<D> flow = mir[j].GetJacobianInverse() * edi.flowip.Row(j);
            flow *= mir[j].GetWeight();	 // weight times Jacobian	
	    edi.flowip.Row(j) = flow;
	  }

	bfi.CalcElementMatrix (fel, trafo, edi.invmass, lh);
	CalcInverse (edi.invmass);

        elementdata.Append (move(edi));
      }




    Array<int> elnums, fnums, vnums;
    
    facetdata.SetAllocSize (ma->GetNFacets());

    for (auto i : Range(ma->GetNFacets()))
      {
	HeapReset hr(lh);
	
	const DGFiniteElement<D-1> & felfacet = 
	  dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));
	IntegrationRule ir (felfacet.ElementType(), 2*felfacet.Order());
        const_cast<DGFiniteElement<D-1>&> (felfacet).PrecomputeShapes (ir);


        FacetData fai(ir.Size());

	ma->GetFacetElements (i, elnums);

	fai.elnr[1] = -1;
	for (int j : Range(elnums))
	  {
	    fai.elnr[j] = elnums[j];
	    auto fnums = ma->GetElFacets (ElementId(VOL,elnums[j]));
	    for (int k : Range(fnums))
	      if (fnums[k] == i) fai.facetnr[j] = k;
	  }

	
	ELEMENT_TYPE eltype = ma->GetElType(ElementId(VOL,elnums[0]));

	vnums = ma->GetElVertices (ElementId(VOL,elnums[0]));
	Facet2ElementTrafo transform(eltype, vnums); 
	FlatVec<D> normal_ref = ElementTopology::GetNormals(eltype) [fai.facetnr[0]];
	
	int nip = ir.Size();
	
	// transform facet coordinates to element coordinates
	IntegrationRule & irt = transform(fai.facetnr[0], ir, lh);  
	MappedIntegrationRule<D,D> mir(irt, ma->GetTrafo(ElementId(VOL, elnums[0]), lh), lh);
	
	FlatMatrixFixWidth<D> flowir(nip, lh);
	cfflow -> Evaluate (mir, flowir);
	
	for (int j = 0; j < nip; j++)
	  {
	    Vec<D> normal = Trans (mir[j].GetJacobianInverse()) * normal_ref;       
	    
	    fai.flown(j) = InnerProduct (normal, flowir.Row(j));
	    fai.flown(j) *= ir[j].Weight() * mir[j].GetJacobiDet();
	  }

	facetdata.Append (move(fai));
      }
    







    BaseVector & vecu = gfu->GetVector();
    auto conv = vecu.CreateVector();
    auto w = vecu.CreateVector();
    auto hu = vecu.CreateVector();

    for (double t = 0; t < tend; t += dt)
      {
        cout << "\rt = " << setw(6) << t << flush;
        
        CalcConvection (vecu.FV<double>(), conv.FV<double>(), lh);
        SolveM (conv.FV<double>(), w.FV<double>(), lh);
          
        hu = vecu + (0.5*dt) * w;

        CalcConvection (hu.FV<double>(), conv.FV<double>(), lh);
        SolveM (conv.FV<double>(), w.FV<double>(), lh);
          
        vecu += dt * w;
        
        /*
          cout << " time T/F/M [us] = "
          << 1e6 * timer_element.GetTime()/timer_element.GetCounts()/vecu.Size() << " / "
          << 1e6 * timer_facet.GetTime()/timer_facet.GetCounts()/vecu.Size() << " / "
          << 1e6 * timer_mass.GetTime()/timer_mass.GetCounts()/vecu.Size() 
          << "\r";
        */
        Ng_Redraw();
      }
  }




  void SolveM (FlatVector<double> res, FlatVector<double> vecu,
	       LocalHeap & lh)
  {
    auto fes = dynamic_pointer_cast<L2HighOrderFESpace> (gfu->GetFESpace());

    timer_mass.Start();

    ParallelFor (Range(ma->GetNE()), 
                 [&] (int i)
                 {
                   IntRange dn = fes->GetElementDofs (i);
                   vecu.Range(dn) = elementdata[i].invmass * res.Range(dn);
                 });

    timer_mass.Stop();
  }






  void CalcConvection (FlatVector<double> vecu, FlatVector<double> conv,
		       LocalHeap & lh)
  {
    
    const L2HighOrderFESpace & fes = 
      dynamic_cast<const L2HighOrderFESpace&> (*gfu->GetFESpace());
    
    
    timer_element.Start();
    
    ParallelFor
      (Range(ma->GetNE()), [&] (int i)
       {
         LocalHeap slh = lh.Split(), &lh = slh;
	 
         auto & fel = static_cast<const ScalarFiniteElement<D>&> (fes.GetFE (ElementId(VOL,i), lh));
         const IntegrationRule ir(fel.ElementType(), 2*fel.Order());
         
         FlatMatrixFixWidth<D> flowip = elementdata[i].flowip;
         
	  /*
	  // use this for time-dependent flow
	  MappedIntegrationRule<D,D> mir(ir, ma->GetTrafo (i, 0, lh), lh);
	  FlatMatrixFixWidth<D> flowip(mir.Size(), lh);
	  cfflow -> Evaluate (mir, flowip);
	  for (int j = 0; j < ir.Size(); j++)
	    {
	      Vec<D> flow = mir[j].GetJacobianInverse() * flowip.Row(j);
	      flow *= mir[j].GetWeight();		
	      flowip.Row(j) = flow;
	    }
	  */

	  IntRange dn = fes.GetElementDofs (i);
	  
	  int nipt = ir.Size();
	  FlatVector<> elui(nipt, lh);
	  FlatMatrixFixWidth<D> flowui (nipt, lh);
	  
	  fel.Evaluate (ir, vecu.Range (dn), elui);
	  
	  for (auto k : Range(nipt))
	    flowui.Row(k) = elui(k) * flowip.Row(k);
	  
	  fel.EvaluateGradTrans (ir, flowui, conv.Range(dn));
       });


      timer_element.Stop();
      timer_facet.Start();

      ParallelFor 
        (Range(ma->GetNFacets()), [&] (int i)
         {
           LocalHeap slh = lh.Split(), &lh = slh;
           
           const FacetData & fai = facetdata[i];
           if (fai.elnr[1] != -1)
             {
               const DGFiniteElement<D> & fel1 = 
                 static_cast<const DGFiniteElement<D>&> (fes.GetFE (ElementId(VOL,fai.elnr[0]), lh));
	      const DGFiniteElement<D> & fel2 = 
		static_cast<const DGFiniteElement<D>&> (fes.GetFE (ElementId(VOL,fai.elnr[1]), lh));
	      const DGFiniteElement<D-1> & felfacet = 
		static_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));

	      IntRange dn1 = fes.GetElementDofs (fai.elnr[0]);
	      IntRange dn2 = fes.GetElementDofs (fai.elnr[1]);

	      int ndoffacet = felfacet.GetNDof();
	      int ndof1 = fel1.GetNDof();
	      int ndof2 = fel2.GetNDof();

	      FlatVector<> aelu1(ndof1, lh), aelu2(ndof2, lh);
	      FlatVector<> trace1(ndoffacet, lh), trace2(ndoffacet, lh);

	      fel1.GetTrace (fai.facetnr[0], vecu.Range (dn1), trace1);
	      fel2.GetTrace (fai.facetnr[1], vecu.Range (dn2), trace2);

	      IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	      int nip = ir.Size();

	      FlatVector<> flown = fai.flown;
	    
	      FlatVector<> tracei1(nip, lh), tracei2(nip, lh);
	      FlatVector<> tracei(nip, lh);

	      felfacet.Evaluate (ir, trace1, tracei1);
	      felfacet.Evaluate (ir, trace2, tracei2);
		    
	      for (int j = 0; j < nip; j++)
		tracei(j) = flown(j) * ( (flown(j) > 0) ? tracei1(j) : tracei2(j) );

	      felfacet.EvaluateTrans (ir, tracei, trace1);
	      fel1.GetTraceTrans (fai.facetnr[0], trace1, aelu1);
	      fel2.GetTraceTrans (fai.facetnr[1], trace1, aelu2);

#pragma omp critical (addres)
	      {
		conv.Range (dn1) -= aelu1;
		conv.Range (dn2) += aelu2;
	      }
	    }
	  else
	    {
	      const DGFiniteElement<D> & fel1 = 
		dynamic_cast<const DGFiniteElement<D>&> (fes.GetFE (ElementId(VOL,fai.elnr[0]), lh));
	      const DGFiniteElement<D-1> & felfacet = 
		dynamic_cast<const DGFiniteElement<D-1>&> (fes.GetFacetFE (i, lh));

	      IntRange dn1 = fes.GetElementDofs (fai.elnr[0]);

	      int ndoffacet = felfacet.GetNDof();
	      int ndof1 = fel1.GetNDof();

	      FlatVector<> elu1(ndof1, lh);
	      FlatVector<> trace1(ndoffacet, lh);

	      fel1.GetTrace (fai.facetnr[0], vecu.Range (dn1), trace1);

	      IntegrationRule ir(felfacet.ElementType(), 2*felfacet.Order());
	      int nip = ir.Size();

	      FlatVector<> flown = fai.flown; 
	      FlatVector<> tracei1(nip, lh), tracei(nip, lh);

	      felfacet.Evaluate (ir, trace1, tracei1);
		    
	      for (int j = 0; j < nip; j++)
		tracei(j) = flown(j) * ( (flown(j) > 0) ? tracei1(j) : 0 );

	      felfacet.EvaluateTrans (ir, tracei, trace1);
	      fel1.GetTraceTrans (fai.facetnr[0], trace1, elu1);
	    
#pragma omp critical (addres)
	      {
		conv.Range (dn1) -= elu1;
	      }
	    }
         });

      timer_facet.Stop(); 
  }
};




static RegisterNumProc<NumProcLinearHyperbolic<2> > npinit1("linhyp", 2);
static RegisterNumProc<NumProcLinearHyperbolic<3> > npinit2("linhyp", 3);

  










