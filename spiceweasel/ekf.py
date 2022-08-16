import numpy as np
from pyrk import RK4
from .jacobian import JacobianCenter


class EKF:
    """
    This is a continuous-discrete extended kalman filter
    """
    def __init__(self, func, dt):
        """
        func(t, x, u)
        x: state (n,1)
        z: measurement (m,1)
        R: measurement noise covariance (m,m)
        Q: process noise covariance (n,n)
        """
        self.func = func
        self.rk = RK4(func, dt)
        self.J = JacobianCenter(self.func)
        self.dt = dt

    def reset(self, n, m):
        """
        Resets the KF to correct dimensions as either an
        eye or zeros matrix.

        n: state size (Q)
        m: measurement size (R)
        """
        self.R = np.eye(m)
        self.Q = np.eye(n)
        self.P = np.eye(n)

        self.x = np.zeros(n)
        self.H = np.eye(m)
        self.I = np.eye(n)

    def predict(self, u):
        """
        Predicts the future state (x) and covarience
        matrix (P). The state predicution uses RK4 to
        integrate.
        """
        self.x = self.rk(self.dt, self.x, u)
        F = self.J(self.dt, self.x, u)
        self.P = F @ self.P @ F.T + self.Q

    def update(self, z):
        """
        Updates the current state estimate using the
        measurement z.
        """
        H = self.H
        I = self.I

        y = z - H.dot(self.x)
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K.dot(y)
        self.P = (I - K.dot(H)).dot(self.P)

        return self.x