#include <solve.hpp>
using namespace ngsolve;
using namespace ngfem;
#include <python_comp.hpp>
#include "1_myFEM/myElement.hpp"
#include "1_myFEM/myFESpace.hpp"
#include "1_myFEM/myPreconditioner.hpp"
#include "1_myFEM/myCoefficient.hpp"
#include "1_myFEM/myIntegrator.hpp"
#include "2_myHOFEM/myHOFESpace.hpp"
#include "4_utility_functions/utility_functions.hpp"


class MyBase
{
protected:
  int count;
public:
  MyBase(int acount) : count(acount) { ; }
  int HowMany() { return count; }
};

class MyClass : public MyBase
{
private:
  string name;
public:
  MyClass(string aname) : MyBase(1), name(aname) { ; }
  string me() { return name; }
};

// a global function ...
void Hello(shared_ptr<MyClass> myclass)
{
  cout << "Hello " << myclass->me() << endl;
}


// Just a demo numproc
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


PYBIND11_MODULE(myngspy,m) {
  // import ngsolve such that python base classes are defined
  auto ngs = py::module::import("ngsolve");

  /*
    In python, classes are objects as well. Therefore we need to create all the class-objects
    we need in our module.
    py::class_ is a template class, which creates python classes.
    The first template parameter of is the class of the object we want to export, the order
    of the other template parameters doesn't matter, for consistency we declare next the holder
    type and then the base class(es).
    The holder type is either unique_ptr (default) or shared_ptr. In most cases we need shared_ptr,
    since we pass our objects to c++ code which may outlive our python objects.
    After that you specify already exported base classes for your object. Note that you do not
    specify not exported bases (since Python does not know anything about them).
    The constructor for our class object takes the module, a string as the class name in Python
    and another string for the docstring of the object. This string is displayed if you
    call the function 'help' on the object.
   */
  py::class_<MyBase, shared_ptr<MyBase>>
    (m, "MyBase", "My base class")
    .def(py::init<int>())
    .def("howmany", &MyBase::HowMany)
    ;

  py::class_<MyClass, shared_ptr<MyClass>, MyBase>
    (m, "MyClass", "My derived class")
    .def(py::init([] (string name)
         {
           return new MyClass(name);
         }), py::arg("name")="me")
    .def("me", [] (shared_ptr<MyClass> self)
         {
           return self->me();
         })
    ;

  /*
    you can export global functions taking shared_ptrs to exported objects easily, pybind
    will do the conversion from python objects to c++ ones for you.
  */
  m.def("Hello", &Hello);


  /*
    If you have old Numproces lying around you can use them from Python by exporting them.
    Note that a NumProc is not needed any more - you can either do all this in Python, or
    if you have some computationally expensive part, write a C++ function or class and export
    it to Python.
  */
  py::class_<NumProcPyDemo,shared_ptr<NumProcPyDemo>>
    (m, "NumProcPyDemo")
    .def ("Hello", &NumProcPyDemo::Hello)
    .def ("Sum", &NumProcPyDemo::Sum)
    ;

  /*
    Finally we export all the other classes and functions we created in this tutorial
   */
  ExportMyFESpace(m);
  ExportMyPreconditioner(m);
  ExportMyCoefficient(m);
  ExportMyIntegrator(m);

  ExportFESpace<MyHighOrderFESpace>(m, "MyHOFESpace");

  m.def("MyAssemble", &myassemble::MyAssemble);
  m.def("MyCoupling", &mycoupling::MyCoupling);

  // This adds documented flags to the docstring of objects (for example
  // FESpaces).
  ngs.attr("_add_flags_doc")(m);
}

static RegisterNumProc<NumProcPyDemo> npinit1("demopy");
