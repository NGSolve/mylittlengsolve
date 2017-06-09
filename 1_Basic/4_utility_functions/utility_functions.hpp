#ifndef MY_UTILITY_FUNCTIONS_HPP
#define MY_UTILITY_FUNCTIONS_HPP

namespace myassemble
{
  shared_ptr<GridFunction> MyAssemble(shared_ptr<FESpace> fes,
                                      shared_ptr<BilinearFormIntegrator> bfi,
                                      shared_ptr<LinearFormIntegrator> lfi);
}

namespace mycoupling
{
  void MyCoupling(shared_ptr<GridFunction> gfu, shared_ptr<LinearForm> lfu);
}

#endif // MY_UTILITY_FUNCTIONS_HPP
