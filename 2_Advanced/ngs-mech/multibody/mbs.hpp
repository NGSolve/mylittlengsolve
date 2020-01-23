#ifndef FILE_MBS
#define FILE_MBS

class MBS_Body
{
private:
  ///
  int domnum;
  ///
  Array<int> nodenums; ///node numbers of domain of body
  ///
  int refpts[3];
  ///
  Mat<3,3> rot;
  ///
  Vec<3> trans;

public:
  ///
  MBS_Body ()
    : nodenums(0)
  {
    ;
  }

  ///
  void SetDomainNr(int d)
  { domnum = d; }

  ///
  int DomainNr() const
  { return domnum; }
  
  ///
  void AddNode(int n)
  {
    for (int i = 0; i < nodenums.Size(); i++)
      if (nodenums[i] == n) { return; }
    nodenums.Append(n);
  }
    
  ///
  int NONodes() const {return nodenums.Size();}

  ///
  int GetNode(int i) const {return nodenums[i]; }

  ///
  void SetRefPoint(int i, int pn) {refpts[i] = pn;}

  ///
  int GetRefPoint(int i) const {return refpts[i];}

  ///
  void SetRot(const Mat<3,3> & r) {rot = r;}
  ///
  void GetRot(Mat<3,3> & r) const {r = rot;}
  ///
  void SetTrans(const Vec<3> & t) {trans = t;}
  ///
  void SetTrans(const Vec<2> & t) 
  {
    trans(0) = t(0);
    trans(1) = t(1);
    trans(2) = 0;
  }
  ///
  const Vec<3> & GetTrans() const { return trans; }
};

  
  
///
typedef enum {TI_MP, TI_R2, TI_R3, TI_R4, TI_IE, TI_L2, TI_L3} TimeIntType;
///
class NumProcDynamic : public NumProc
{
  enum { DIM = 3 };
  ///
  typedef Array<MBS_Body*> MBS_System;
  
  ///
  BilinearForm * bfa;
  ///
  BilinearForm * bfm;
  ///
  GridFunction * gfu;
  ///
  GridFunction * gfusmall;
  ///
  GridFunction * gfv;
  ///
  double te, dt;
  ///
  TimeIntType tit;
  ///
  int nbnds;
  ///
  Array<int> point2domain;
  ///
  Array<int> bnd2domain;
  ///
  Array<int[2]> connections;
  ///
  Array<int[2]> connectionsx;
  ///
  Array<int[2]> connectionsy;
  ///
  Array<int[2]> connectionsz;

  /// linear functions on individual bodies:  
  Array<VVector<Vec<DIM> > *> pol1;
  ///
  Array<VVector<Vec<DIM> > *> mpol1;
  ///
  Matrix<> masslin, invmasslin;

  ///
  Array<VVector<Vec<DIM> > *> cnts;
  ///
  Matrix<> cnts2bnds;
  
  ///
  VVector<Vec<DIM> > points; /// faster access to mesh points

public:
///
  NumProcDynamic (PDE & apde, const Flags & flags);
  ///
  virtual ~NumProcDynamic();

  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcDynamic (pde, flags);
  }
  ///
  virtual void Do(LocalHeap & lh);
  ///
  virtual void PrintReport (ostream & ost);

protected:

  ///
  void PrepareSystem (MBS_System & mbssystem);

  /// computes V'(u)
  double ComputePotentialPrime (MBS_System & mbssystem,
                                const VVector<Vec<DIM> > & u,
                                const VVector<Vec<DIM> > & ulin,
                                VVector<Vec<DIM> > & prime,
                                int simplified = 0);
    
  /// (A+\tau^2 M) - orthogonal projections into sub-space of constraints
  void Project2Constraints (MBS_System & mbssystem,
                            const BaseMatrix & invhmat,
                            Matrix<double> & bnd_int_schur,
                            VVector<Vec<DIM> > & v,
                            Vector<double> & lami);
  ///
  void CalcLinearApproximation (const VVector<Vec<3> > & u, 
                                VVector<Vec<3> > & ulin);

  ///
  void CalcRMatrix (const VVector<Vec<DIM> > & u, 
                    MBS_System & mbssystem);

  ///
  void SetLargeDeformation (VVector<Vec<DIM> > & u,
                            const MBS_System & mbssystem); 
  ///
  void ApplyRotation (const VVector<Vec<DIM> > & x,
                      VVector<Vec<DIM> > & y, 
                      const MBS_System & mbssystem,
                      int inverse = 0);

  ///
  void ComputeFCorrection (MBS_System & mbssystem,
                           const VVector<Vec<DIM> > & ularge,
                           const VVector<Vec<DIM> > & kusmall,
                           const VVector<Vec<DIM> > & rkusmall,
                           const VVector<Vec<DIM> > & usmallrot,
                           VVector<Vec<DIM> > & f2);

  ///
  void WriteData (double t,
                  const VVector<Vec<DIM> > & u,
                  MBS_System & mbssystem);



  void ComputeU(const FlatVector<double> & point,
                FlatVector<double> & up,
                const S_GridFunction<double>& bu,
                LocalHeap & lh);
    
};



#endif
