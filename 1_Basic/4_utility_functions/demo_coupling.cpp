/*
  

In this demo, we take the one solution (GridFunction) as right
hand side of a second one.
Note that this can be done easier in Python by just using the GridFunction as a
CoefficientFunction in the integrator.

I.e. we  want to solve the first PDE,

-\Delta u = f

the second one is

-\Delta w = u


Similar problems occur when coupling different physical fields.

The Function 'MyCoupling' gets a gridfunction from the solution space and
a LinearForm.
It will use the GridFunction to build the right hand side of the equation.


*/

#include <solve.hpp>
using namespace ngsolve;
#include "utility_functions.hpp"


namespace mycoupling
{
  void MyCoupling(shared_ptr<GridFunction> gfu, shared_ptr<LinearForm> lfu)
  {
    cout << "Compute coupling terms" << endl;
    Flags lfuFlags;
    auto fes_gfu = gfu->GetFESpace();
    auto fes_lfu = lfu->GetFESpace();
    auto ma = fes_gfu->GetMeshAccess();
    lfu->GetVector() = 0.0;

    // create local heap for efficient memory management
    LocalHeap lh(100000);

    for(auto el : ma->Elements(VOL))
      {
        HeapReset hr(lh);    // reset the local heap memory at the end of the loop
		
        const ElementTransformation & eltrans = ma->GetTrafo (el, lh);

        const ScalarFiniteElement<2> & fel_gfu =
	  dynamic_cast<const ScalarFiniteElement<2>&> (fes_gfu->GetFE (el, lh));
        const ScalarFiniteElement<2> & fel_lfu =
	  dynamic_cast<const ScalarFiniteElement<2>&> (fes_lfu->GetFE (el, lh));

	int nd_gfu = fel_gfu.GetNDof();
	int nd_lfu = fel_lfu.GetNDof();

	Array<int> dnums_gfu(nd_gfu, lh), dnums_lfu(nd_lfu, lh);  // the dof-numbes for gfu and lfu
		      
        fes_gfu->GetDofNrs (el, dnums_gfu);
        fes_lfu->GetDofNrs (el, dnums_lfu);

        FlatVector<> el_gfu (nd_gfu, lh), shape (nd_lfu, lh), el_lfu (nd_lfu, lh);

        gfu -> GetElementVector (dnums_gfu, el_gfu);

        IntegrationRule ir(eltrans.GetElementType(), fel_gfu.Order()+fel_lfu.Order() );

        
        el_lfu = 0.0;
        for (int j = 0; j < ir.GetNIP(); j++)   // integration points
          {
            MappedIntegrationPoint<2, 2> mip(ir[j], eltrans);  // computes Jacobi matrix etc

            Vec<1> ui;   // value of u in point

            // compute value in point (= el_gfu * shape)
            DiffOpId<2>::Apply (fel_gfu, mip, el_gfu, ui, lh);
            
            // could use also other differential operators such as
            // DiffOpGradient<2>, DiffOpCurl<2>, ....
 
            double fac = ir[j].Weight() * mip.GetMeasure();   // integration weights
	    
	    fel_lfu.CalcShape (ir[j], shape);
            el_lfu += (fac*ui(0)) * shape;
          }

        lfu -> AddElementVector (dnums_lfu, el_lfu);
      }
  }
}
