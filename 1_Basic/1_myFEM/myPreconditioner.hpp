#ifndef MY_PRECONDITIONER_HPP
#define MY_PRECONDITIONER_HPP


namespace ngcomp
{
  class MyPreconditioner : public Preconditioner
  {
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BaseJacobiPrecond> jacobi;

  public:
    MyPreconditioner (shared_ptr<BilinearForm> abfa, Flags& flags);
    void Update() override;
    
    void Mult (const BaseVector & f, BaseVector & u) const override
    {
      jacobi -> Mult (f, u);

      /*
      u = 0.0;
      jacobi -> GSSmooth (u, f);
      jacobi -> GSSmoothBack (u, f);
      */
    }

    const BaseMatrix& GetMatrix() const override
    {
      return *this;
    }

    const BaseMatrix & GetAMatrix() const override
    {
      return bfa -> GetMatrix();
    }
  };
}

void ExportMyPreconditioner(py::module m);

#endif // MY_PRECONDITIONER_HPP
