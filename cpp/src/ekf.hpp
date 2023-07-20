
#pragma once

using ekfFunc = std::function<cevd(cevd &, const double)>;

class Jacobian {
public:
  Jacobian(ekfFunc f, size_t n) : n(n), f(f), J(Eigen::MatrixXd(n, n)) {}

  Eigen::MatrixXd get(const Eigen::VectorXd &x, const double t) {
    Eigen::VectorXd xx = f(x, t);
    double dy          = 1e-8;

    for (int i = 0; i < n; ++i) {
      Eigen::VectorXd y = x;
      y(i) += dy;
      Eigen::VectorXd yy = f(y, t);
      J.col(i) << (yy - xx) / dy;
    }
    return J;
  }

protected:
  ekfFunc f;
  Eigen::MatrixXd J;
  const size_t n;
};

class ExtendedKalmanFilter {
public:
  ExtendedKalmanFilter(ekfFunc f, cemd &H, cemd &Q, cemd &R, cemd &P, size_t n)
      : H(H), Q(Q), R(R), P(P), x(n), I(Eigen::MatrixXd::Identity(n, n)),
        jac(Jacobian(f, n)), F(Eigen::MatrixXd(n, n)) {}

  Eigen::VectorXd update(cevd &z, cevd &u, double t) {
    F = jac.get(x, t);
    P = F * P * F.transpose() + Q;
    K = P * H.transpose() * (H * P * H.transpose() + R).inverse();
    x = F * x + K * (z - H * x);
    P = (I - K * H) * P;

    return x;
  }

protected:
  ekfFunc f;
  Jacobian jac;
  const Eigen::MatrixXd H, Q, R;
  Eigen::MatrixXd P, K, F;
  const Eigen::MatrixXd I;
  Eigen::VectorXd x;
};