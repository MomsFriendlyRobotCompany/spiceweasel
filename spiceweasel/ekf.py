# -*- coding: utf-8 -*-
##############################################
# The MIT License (MIT)
# Copyright (c) 2022 Kevin Walchko
# see LICENSE for full details
##############################################
import numpy as np
from pyrk import RK4
from .jacobian import JacobianCenter


class EKF:
    """
    This is a continuous-discrete extended kalman filter
    Gleb, Table 6.1-1, p 188

    func(dt, x, u)

    x: state (n,1)
    z: measurement (m,1)
    P: covariance (n,n)
    R: measurement noise covariance (m,m)
    Q: process noise covariance (n,n)

    """
    def __init__(self, func, dt, n, m):
        """
        func(dt, x, u)
        dt: time step
        n: state size
        m: measurement size
        """
        self.func = func
        self.rk = RK4(func, dt)
        self.J = JacobianCenter(self.func)
        self.dt = dt
        self.m = m
        self.n = n
        self.reset()

    def reset(self):
        """
        Resets the KF to correct dimensions as either an
        eye or zeros matrix.
        """
        m = self.m
        n = self.n

        self.R = np.eye(m)
        self.Q = np.eye(n)
        self.P = np.eye(n)

        self.x = np.zeros(n)
        self.H = np.eye(m) # should this be n, even though m < n?
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
        # self.P = F @ self.P + self.P @ F.T + self.Q

    def update(self, z):
        """
        Updates the current state estimate using the
        measurement z.
        """
        H = self.H
        I = self.I

        S = H @ self.P @ H.T + self.R
        # change to: solve K*S = (P*H.T) -> Ax=B, faster, more stable?
        K = self.P @ H.T @ np.linalg.inv(S)

        y = z - H.dot(self.x)
        self.x = self.x + K.dot(y)

        self.P = (I - K.dot(H)).dot(self.P)

        return self.x