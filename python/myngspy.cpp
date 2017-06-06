#include <solve.hpp>
using namespace ngsolve;
#include <python_ngstd.hpp>
#include "../myElement.hpp"
#include "../myFESpace.hpp"
#include "../myCoefficient.hpp"
#include "../myDiffop.hpp"



class MyBase
{
protected:
  int count;
public:
  MyBase(int acount) : count(acount) { ; }
  int howmany() { return count; }
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


PYBIND11_PLUGIN(myngspy) {
  py::module m("myngspy", "myngspy doc-string");
  
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
    specify not exported bases (since Python does not have know anything about them).
    The constructor for our class object takes the module, a string as the class name in Python
    and another string for the docstring of the object. This string is displayed if you
    call the function 'help' on the object.
   */
  py::class_<MyBase, shared_ptr<MyBase>>
    (m, "MyBase", "My base class")
    .def(py::init<int>())
    .def("howmany", &MyBase::howmany)
    ;

  py::class_<MyClass, shared_ptr<MyClass>, MyBase>
    (m, "MyClass", "My derived class")
    .def("__init__", [] (MyClass* instance, string name)
         {
           new (instance) MyClass(name);
         },
         py::arg("name")="me")
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
  ExportMyElement(m);
  ExportMyFESpace(m);
  ExportMyCoefficient(m);
  ExportMyDiffOp(m);

  return m.ptr();
}




static RegisterNumProc<NumProcPyDemo> npinit1("demopy");
