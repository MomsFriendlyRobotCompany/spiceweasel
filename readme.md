![](https://github.com/MomsFriendlyRobotCompany/spiceweasel/raw/main/pics/elzar.png)

![GitHub](https://img.shields.io/github/license/MomsFriendlyRobotCompany/spiceweasel)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spiceweasel)
![PyPI](https://img.shields.io/pypi/v/spiceweasel)
![PyPI - Downloads](https://img.shields.io/pypi/dm/spiceweasel?color=aqua)

# Spice Weasel

>  I knocked it up a notch with my `spiceweasel`. Bam!

**under development**

The list of imported packages can be found in the [pyproject.toml][toml].

## Usage

Basic example:

```python
import numpy as np
from spiceweasel import EKF

def func(dt, x, u):
    """
    dt: time step
    x: state estimate
    u: control forces or other inputs
    """

    # differential equations
    return x

# create a kalman filter
ekf = EKF(func, dt, 2, 2)

# so reset puts R and Q to identify matrix, you should
# adjust them to your system
ekf.reset()
ekf.R *= [0.01,0.01,0.1] # measurement cov
ekf.Q *= [.05,.05,.1]    # process cov
ekf.x = np.array([1,-2]) # default sets this to zeros

# main filtering loop
for i in range(num):
    # ...
    ekf.predict(u)
    # ...
    y = ekf.update(meas)

# ...
```

# MIT License

**Copyright (c) 2022 Mom's Friendly Robot Company**

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

[toml]: https://github.com/MomsFriendlyRobotCompany/spiceweasel/blob/main/pyproject.toml
