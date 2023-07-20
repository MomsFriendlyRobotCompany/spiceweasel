
#include <iostream>

#include <spiceweasel.hpp>

using namespace std;
using namespace Eigen;

/*
xdd = f/m
xd = v
*/
const VectorXd func(const VectorXd &x, double t) {
  VectorXd xx(2);
  xx(0) = 0.01 / 25;
  xx(1) = x[0];
  return xx;
}

int main() {
  size_t m  = 1; // measurements
  size_t n  = 2; // states
  double t  = 0.0;
  double dt = 0.1;
  MatrixXd C(m, n); // Output matrix
  MatrixXd Q(n, n); // Process noise covariance
  MatrixXd R(m, m); // Measurement noise covariance
  MatrixXd P(n, n); // Estimate error covariance
  VectorXd x(n), u(n);

  C << 1, 0;
  Q << 0.1, 0, 0, 0.1;
  R << 5;
  P << 0.1, 0, 0, 0.1;

  cout << C << endl;
  cout << Q << endl;
  cout << R << endl;
  cout << P << endl;
  cout << x.transpose() << endl;
  cout << u.transpose() << endl;

  ExtendedKalmanFilter ekf(func, C, Q, R, P, n);

  for (int i = 0; i < 100; ++i) {
    VectorXd z(1);
    x = ekf.update(z, u, t);
    t += dt;
    cout << t << ": " << x << endl;
  }
  return 0;
}