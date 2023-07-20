

#pragma once

class KalmanFilter {
public:
  /**
   * Discrete Kalman filter
   *   A(n,n) - System dynamics matrix
   *   H(m,n) - Output matrix
   *   Q(n,n) - Process noise covariance
   *   R(m,m) - Measurement noise covariance
   *   P(n,n) - Estimate error covariance
   *   I(n,n) - Identity matrix
   */
  KalmanFilter(cemd &A, cemd &H, cemd &Q, cemd &R, cemd &P)
      : A(A), H(H), Q(Q), R(R), P(P), x(A.rows()),
        I(Eigen::MatrixXd::Identity(A.rows(), A.rows())) {}

  void init() { x.setZero(); }
  void init(cevd &x0) { x = x0; }

  cevd update(cevd &z) {
    P = A * P * A.transpose() + Q;
    K = P * H.transpose() * (H * P * H.transpose() + R).inverse();
    x = A * x + K * (z - H * x);
    P = (I - K * H) * P;

    return x;
  }

protected:
  const Eigen::MatrixXd A, H, Q, R;
  Eigen::MatrixXd P, K;
  const Eigen::MatrixXd I;
  Eigen::VectorXd x;
};
