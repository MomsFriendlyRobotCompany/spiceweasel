
#include <cmath>
#include <iostream>
#include <spiceweasel.hpp>

using namespace std;
using namespace Eigen;

/*
Python Jupyter Result:
[[19.99999995 24.99999994]
 [ 4.99999999 -0.41614676]]
*/
VectorXd func(VectorXd x, double t) {
  VectorXd ret(2);
  ret(0) = x[0] * x[0] * x[1];
  ret(1) = 5.0 * x[0] + sin(x[1]);

  return ret;
}

int main() {
  Jacobian jac(func, 2);
  VectorXd x(2, 1);
  x << 5, 2;

  cout << jac.get(x, 0.0) << endl;

  return 0;
}