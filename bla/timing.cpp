#include <bla.hpp>
using namespace ngbla;


int main()
{
  cout << "# n   MFlops" << endl;
  
  for (size_t n = 4; n <= 2048; n *= 2)
    {
      cout << n << "   " << flush;
      Timer t("mat-mat "+ToString(n));
      Matrix a(n,n), b(n,n), c(n,n);
      a = 1; b = 2;

      size_t m = size_t(1e9 / (n*n*n)) + 1;

      t.Start();
      for (size_t j = 0; j < m; j++)
        c = a*b;
      t.Stop();
      t.AddFlops (m*n*n*n);
      cout << t.GetMFlops() << endl;
    }
}
