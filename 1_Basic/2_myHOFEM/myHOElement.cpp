/*********************************************************************/
/* File:   myElement.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
It is also simple to implement high order elements ...

*/


#include <fem.hpp>
#include "myHOElement.hpp"


namespace ngfem
{

  MyHighOrderSegm :: MyHighOrderSegm (int order)
    : MyBaseElement (order+1, order)
  { ; }


  void MyHighOrderSegm :: CalcShape (const IntegrationPoint & ip, 
                                     BareSliceVector<> shape) const
  {
    double x = ip(0);
    T_CalcShape (x, shape);
  }


  void MyHighOrderSegm :: CalcDShape (const IntegrationPoint & ip, 
                                      BareSliceMatrix<> dshape) const
  {
    AutoDiff<1> adx (ip(0), 0);
    Vector<AutoDiff<1> > shapearray(ndof);
    T_CalcShape<AutoDiff<1>> (adx, shapearray);
    for (int i = 0; i < ndof; i++)
      dshape(i, 0) = shapearray[i].DValue(0);
  }

  template <class T>
  void MyHighOrderSegm :: T_CalcShape (const T& x, BareSliceVector<T> shape) const
  {
    T lami[2] = { x, 1-x };
    
    for (int i = 0; i < 2; i++)
      shape[i] = lami[i];

    int ii = 2;
    
    ArrayMem<T, 20> polx(order+1);

    if (order >= 2)  
      {
        INT<2> edge = GetVertexOrientedEdge(0);
        IntegratedLegendrePolynomial (order, 
                                      lami[edge[1]]-lami[edge[0]],
                                      polx);
        for (int j = 2; j <= order; j++)
          shape[ii++] = polx[j];
      }
  }

  MyHighOrderTrig :: MyHighOrderTrig (int order)
    : MyBaseElement ((order+1)*(order+2)/2, order)
  { ; }

  void MyHighOrderTrig :: CalcShape (const IntegrationPoint & ip, 
                                     BareSliceVector<> shape) const
  {
    double x = ip(0);
    double y = ip(1);
    T_CalcShape (x, y, shape);
  }


  void MyHighOrderTrig :: CalcDShape (const IntegrationPoint & ip, 
                                      BareSliceMatrix<> dshape) const
  {
    AutoDiff<2> adx (ip(0), 0);
    AutoDiff<2> ady (ip(1), 1);
    Vector<AutoDiff<2> > shapearray(ndof);
    T_CalcShape<AutoDiff<2>> (adx, ady, shapearray);
    for (int i = 0; i < ndof; i++)
      {
        dshape(i, 0) = shapearray[i].DValue(0);
        dshape(i, 1) = shapearray[i].DValue(1);
      }
  }




  template <class T>
  void MyHighOrderTrig :: T_CalcShape (const T & x, const T & y, BareSliceVector<T> shape) const
  {
    T lami[3] = { x, y, 1-x-y };
    
    for (int i = 0; i < 3; i++)
      shape[i] = lami[i];

    int ii = 3;
    
    ArrayMem<T, 20> polx(order+1), poly(order+1);

    for (int i = 0; i < 3; i++)
      if (order >= 2)   // more general: order_edge[i] 
	{
          auto edge = GetVertexOrientedEdge(i);
          ScaledIntegratedLegendrePolynomial (order, 
                                              lami[edge[1]]-lami[edge[0]],
                                              lami[edge[0]]+lami[edge[1]], polx);
          for (int j = 2; j <= order; j++)
            shape[ii++] = polx[j];
	}
    
    // inner dofs
    if (order >= 3)      // more general: cell order
      {
        T bub = x * y * (1-x-y);
        ScaledLegendrePolynomial (order-2, lami[1]-lami[0], lami[1]+lami[0], polx);
        LegendrePolynomial (order-1, 2*lami[2]-1, poly);

        for (int i = 0; i <= order-3; i++)
          for (int j = 0; j <= order-3-i; j++)
            shape[ii++] = bub * polx[i] * poly[j];
      }
  }
}
