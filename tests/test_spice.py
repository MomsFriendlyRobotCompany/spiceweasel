import pytest
from spiceweasel import *


def test_ekf():
    def f(dt, x, u):
        return x

    try:
        ekf =  EKF(f, 0.1,2,2)
        ekf.reset()
        assert True
    except:
        assert False