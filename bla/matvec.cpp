#include <bla.hpp>
using namespace ngbla;


int main()
{
  Vector<double> x(5);

  x = 3;
  x.Range(1,3) = 5;
  cout << "x = " << x << endl;
  
  Vector<double> y = 3*x;
  cout << "y = " << y << endl;

  Matrix<double> a(5,5);
  a = Identity(5);
  a.Rows(2,4).Cols(1,3) = 2;
  a.Col(3) = -1;
  cout << "a = " << endl << a << endl;
  y = a*x;
  cout << "y = " << y << endl;
}


