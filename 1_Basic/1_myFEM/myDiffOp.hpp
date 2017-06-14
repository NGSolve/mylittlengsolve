#ifndef FILE_MYDIFFOP_HPP
#define FILE_MYDIFFOP_HPP

#include <solve.hpp>

namespace ngfem
{
  class MyIdentity : public DiffOp<MyIdentity>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      HeapReset hr(lh);
      mat.Row(0) = static_cast<const ScalarFiniteElement<2>&>(fel).GetShape(mip.IP(), lh);
    }

  };

}

#endif // FILE_MYDIFFOP_HPP
