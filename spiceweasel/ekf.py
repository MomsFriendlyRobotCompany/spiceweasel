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
        self.R = np.eye(m)
        self.Q = np.eye(n)
        self.P = np.eye(n)

        self.x = np.zeros(n)
        self.H = np.eye(m)
        self.I = np.eye(n)

    def predict(self, u):
        self.x = self.rk(self.dt, self.x, u)
        F = self.J(self.dt, self.x, u)
        self.P = F @ self.P @ F.T + self.Q

    def update(self, z):
        H = self.H
        I = self.I

        y = z - H.dot(self.x)
        # S = H.dot(self.P.dot(H.T)) + self.R
        S = H @ self.P @ H.T + self.R
        # K = self.P.dot(H.T.dot(np.linalg.inv(S)))
        K = self.P @ H.T @ np.linalg.inv(S)
        self.x = self.x + K.dot(y)
        self.P = (I - K.dot(H)).dot(self.P)

        # q = Quaternion(*self.x[6:])
        # if abs(1 - q.magnitude) > 1e-6:
        #     q.normalize
        #     self.x[6] = q.w
        #     self.x[7] = q.x
        #     self.x[8] = q.y
        #     self.x[9] = q.z

        # self.x = x
        # self.P = P

        return self.x