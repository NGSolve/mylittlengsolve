#include <solve.hpp>
using namespace ngsolve;

#include <python_ngstd.hpp>


// a global function ...
void Hello()
{
  cout << "Hello everybody" << endl;
}




class NumProcPyDemo : public NumProc
{

public:
    
  NumProcPyDemo (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde) { ; }

  virtual void Do(LocalHeap & lh)
  {
    cout << "solving NumProcPyDemo" << endl;
  }
  
  virtual void PrintReport (ostream & ost)
  {
    ost << "I am NumProcPyDemo" << endl;
  }

  virtual void Hello (string name)
  {
    cout << "Hi " << name << endl;
  }

  virtual double Sum (double a, double b)
  {
    return a+b;
  }
};


static RegisterNumProc<NumProcPyDemo> npinit1("demopy");




PYBIND11_PLUGIN(myngspy) {
  py::module m("myngspy", "myngspy doc-string");
  
  m.def("Hello", &Hello);

  py::class_<NumProcPyDemo,shared_ptr<NumProcPyDemo>>
    (m, "NumProcPyDemo")

    .def ("Hello", &NumProcPyDemo::Hello)
    .def ("Sum", &NumProcPyDemo::Sum)
    ;

  return m.ptr();
}



struct Init {
  Init() 
  { 
    cout << endl << endl 
         << "************ please execute in py-shell: *************" 
         << endl << endl << endl

         << "import myngspy\n"
         << "np = pde.numprocs['np1']\n"
         << "print (np)\n"
         << "np.Hello ('Matthias')\n"
         << "np.Sum (3,8)"

         << endl << endl << endl;
  }
};
static Init init;





