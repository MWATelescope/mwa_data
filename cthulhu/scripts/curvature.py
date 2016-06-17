from scipy.interpolate import UnivariateSpline
import numpy as np

def curvature_splines(x, y=None, error=0.1):
    """Calculate the signed curvature of a 2D curve at each point
    using interpolating splines.

    Parameters
    ----------
    x,y: numpy.array(dtype=float) shape (n_points, )
         or
         y=None and
         x is a numpy.array(dtype=complex) shape (n_points, )

         In the second case the curve is represented as a np.array
         of complex numbers.

    error : float
        The admisible error when interpolating the splines

    Returns
    -------
    curvature: numpy.array shape (n_points, )

    Note: This is 2-3x slower (1.8 ms for 2000 points) than `curvature_gradient`
    but more accurate, especially at the borders.
    """

    # handle list of complex case
    if y is None:
        x, y = x.real, x.imag

    t = np.arange(x.shape[0])
    std = error * np.ones_like(x)

    fx = UnivariateSpline(t, x, k=4, w=1 / np.sqrt(std))
    fy = UnivariateSpline(t, y, k=4, w=1 / np.sqrt(std))

    xdash = fx.derivative(1)(t)
    xdashdash = fx.derivative(2)(t)
    ydash = fy.derivative(1)(t)
    ydashdash = fy.derivative(2)(t)
    curvature = (xdash* ydashdash - ydash* xdashdash) / np.power(xdash** 2 + ydash** 2, 3 / 2)
    return curvature
