#!/usr/bin/env python3
'''ErrorSurface class for containing resampled fitting errors'''

__author__ = "Juha Tiihonen"
__email__ = "tiihonen@iki.fi"
__license__ = "BSD-3-Clause"

from bisect import bisect
from scipy.interpolate import LinearNDInterpolator
from numpy import array, argsort, append, where, isscalar


class ErrorSurface():
    _E_mat = None  # Matrix of total errors
    _X_mat = None  # X-mesh
    _Y_mat = None  # Y-mesh
    X_res = None  # fractional X resolution
    Y_res = None  # fractional Y resolution
    verbosity = None  # Verbosity level of output:
    # 0 -> no output; 1 -> only critical errors; 2 -> all output

    def __init__(
        self,
        X_res=0.1,
        Y_res=0.1,
        verbosity=1
    ):
        # Initialize with zero error at the origin
        self._E_mat = array([[0.0]])
        self._X_mat = array([[0.0]])
        self._Y_mat = array([[0.0]])
        self.X_res = X_res
        self.Y_res = Y_res
        self.verbosity = verbosity
    # end def

    @property
    def Xs(self):
        return self._X_mat[0]
    # end def

    @property
    def Ys(self):
        return self._Y_mat[:, 0]
    # end def

    @property
    def X_mat(self):
        return self._X_mat
    # end def

    @property
    def Y_mat(self):
        return self._Y_mat
    # end def

    @property
    def E_mat(self):
        return self._E_mat
    # end def

    @property
    def T_mat(self):
        return self._X_mat >= self._Y_mat
    # end def

    def insert_row(self, y, row):
        if len(row) != len(self.Xs):
            raise ValueError(f"Cannot add row with len={len(row)} to data with len={len(self.Xs)}")
        # end if
        if not isscalar(y) or (y < 0.0):
            raise ValueError("Cannot add y < 0.0.")
        # end if
        X_mat = append(self._X_mat, [self.Xs], axis=0)
        Y_mat = append(self._Y_mat, [len(self.Xs) * [y]], axis=0)
        E_mat = append(self._E_mat, [row], axis=0)
        idx = argsort(Y_mat[:, 0])
        self._X_mat = X_mat[idx]
        self._Y_mat = Y_mat[idx]
        self._E_mat = E_mat[idx]
    # end def

    def insert_col(self, x, col):
        if len(col) != len(self.Ys):
            raise ValueError(f"Cannot add row with len={len(col)} to data with len={len(self.Ys)}")
        # end if
        if not isscalar(x) or (x < 0.0):
            raise ValueError("Cannot add x < 0.0.")
        # end if
        X_mat = append(self._X_mat, array([len(self.Ys) * [x]]).T, axis=1)
        Y_mat = append(self._Y_mat, array([self.Ys]).T, axis=1)
        E_mat = append(self._E_mat, array([col]).T, axis=1)
        idx = argsort(X_mat[0])
        self._X_mat = X_mat[:, idx]
        self._Y_mat = Y_mat[:, idx]
        self._E_mat = E_mat[:, idx]
    # end def

    def evaluate_target(self, epsilon):
        xi, yi = self._argmax_y(epsilon)
        # TODO: the performance can be further improved by interpolating to maximize epsilon
        return self.Xs[xi], self.Ys[yi]
    # end def

    def evaluate_surface(self, x, y):
        xi = bisect(self.Xs, x)
        yi = bisect(self.Ys, y)

        if xi < len(self.Xs):
            xi_prev = xi - 1
        else:
            if self.verbosity >= 2:
                print(f"  Requested x>={self.Xs[-1]}, reverting to x={self.Xs[-1]}")
            # end if
            xi, xi_prev = xi - 1, xi - 1
        # end if

        if yi < len(self.Ys):
            yi_prev = yi - 1
        else:
            if self.verbosity >= 2:
                print(f"  Requested y>={self.Ys[-1]}, reverting to y={self.Ys[-1]}")
            # end if
            yi, yi_prev = yi - 1, yi - 1
        # end if

        pts = [
            [self.Xs[xi], self.Ys[yi]],
            [self.Xs[xi_prev], self.Ys[yi]],
            [self.Xs[xi], self.Ys[yi_prev]],
            [self.Xs[xi_prev], self.Ys[yi_prev]]
        ]
        vals = [
            self.E_mat[yi, xi],
            self.E_mat[yi, xi_prev],
            self.E_mat[yi_prev, xi],
            self.E_mat[yi_prev, xi_prev]
        ]
        val_int = LinearNDInterpolator(pts, vals)
        val = val_int([x, y])
        return val
    # end def

    def request_points(self, epsilon):
        xi, yi = self._argmax_y(epsilon)

        x_vals = []
        y_vals = []
        if xi == 0:
            x_vals += self._treat_x_underflow()
        elif xi == len(self.Xs) - 1:
            x_vals += self._treat_x_overflow()
        else:
            x_vals += self._treat_x_res(xi)
        # end if

        if yi == 0:
            y_vals += self._treat_y_underflow()
        elif yi == len(self.Ys) - 1:
            y_vals += self._treat_y_overflow()
        else:
            y_vals += self._treat_y_res(yi)
        # end if

        return x_vals, y_vals
    # end def

    def _argmax_y(self, epsilon):
        """Return indices to the highest point in E matrix that is lower than epsilon"""
        if (epsilon <= 0.0):
            raise ValueError("Cannot optimize to epsilon <= 0.0")
        # end if
        xi, yi = 0, 0
        for i in range(len(self.E_mat), 0, -1):  # from high to low
            err = where((self.E_mat[i - 1] < epsilon) & (self.T_mat[i - 1]))[0]
            if len(err) > 0:
                yi = i - 1
                # pick xi from the middle of the available values
                xi = err[int(len(err) / 2)]
                break
            # end if
        # end for
        return xi, yi
    # end def

    # Fix x underflow by adding a new X value between the first and second
    def _treat_x_underflow(self):
        # The first x is zero
        X_this = self.Xs[0]
        X_right = self.Xs[1]
        # Add new x value to the right of 0
        X_new = (X_this + X_right) / 2
        X_diff = (X_right - X_new) / X_right
        # This cannot be False unless X_res >= 0.5, which would be foolish
        if X_diff > self.X_res:
            return [X_new]
        else:
            if self.verbosity >= 1:
                msg = f"  Persistent x-underflow. Could not add x={X_new}. "
                msg += "Check the data behind error surface: "
                msg += f"E(x={X_right}, y=0)={self.E_mat[1, 0]}."
                print(msg)
            # end if
            return []
        # end if
    # end def

    def _treat_x_overflow(self):
        X_this = self.Xs[-1]
        X_left = self.Xs[-2]
        # Add new x value to the left of x-max
        X_new = (X_this + X_left) / 2
        X_diff = (X_this - X_new) / X_this
        if X_diff > self.X_res:
            return [X_new]
        else:
            if self.verbosity >= 1:
                msg = f"  Persistent x-overflow. Did not add x={X_new} "
                msg += f"next to x-max={X_this} to maintain x-resolution={self.X_res}. "
                msg += f"Gain performance by adding data beyond x>{X_this}."
                print(msg)
            # end if
            return []
        # end if
    # end def

    def _treat_x_res(self, xi):
        X_this = self.Xs[xi]
        X_left = self.Xs[xi - 1]
        X_right = self.Xs[xi + 1]
        X_new_left = (X_this + X_left) / 2
        X_new_right = (X_this + X_right) / 2
        vals = []
        if (X_this - X_new_left) / X_this > self.X_res:
            # Add new X value to the left
            vals += [X_new_left]
        # end if
        if (X_right - X_this) / X_right > self.X_res:
            # Add new X value to the right
            vals += [X_new_right]
        # end if
        return vals
    # end def

    # Fix x underflow by adding a new sigma value between the first and second
    def _treat_y_underflow(self):
        # The first y is zero
        Y_this = self.Ys[0]
        Y_up = self.Ys[1]
        # Add new y value above zero
        Y_new = (Y_this + Y_up) / 2
        Y_diff = (Y_up - Y_new) / Y_up
        # This cannot be False unless Y_res >= 0.5, which would be foolish
        if Y_diff > self.Y_res:
            return [Y_new]
        else:
            if self.verbosity >= 1:
                msg = f"  Persistent y-underflow. Did not add y={Y_new} "
                msg += f"below y={Y_up} to maintain y-resolution={self.Y_res}. "
                print(msg)
            # end if
            return []
        # end if
    # end def

    # Fix y overflow by adding a new sigma value twice as high until max(Xs)
    def _treat_y_overflow(self):
        Y_this = self.Ys[-1]
        X_max = self.Xs[-1]
        Y_new = 2 * Y_this
        if Y_new > X_max:
            if self.verbosity >= 2:
                print(f"  Capping y-new to x-max={X_max}.")
            # end if
            Y_new = X_max
        # end if
        Y_diff = (Y_new - Y_this) / Y_new
        if Y_diff > self.Y_res:
            return [Y_new]
        else:
            if self.verbosity >= 1:
                msg = f"  Persistent y-overflow. Did not add y={Y_new} "
                msg += f"between y={Y_this} and x-max={X_max} "
                msg += f"to maintain y-resolution={self.Y_res}. "
                print(msg)
            # end if
            return []
        # end if
    # end def

    def _treat_y_res(self, yi):
        Y_this = self.Ys[yi]
        Y_down = self.Ys[yi - 1]
        Y_up = self.Ys[yi + 1]
        Y_new_up = (Y_this + Y_up) / 2
        Y_new_down = (Y_this + Y_down) / 2
        vals = []
        if (Y_new_up - Y_this) / Y_new_up > self.Y_res:
            # Add new Y value above
            vals += [Y_new_up]
        # end if
        if (Y_this - Y_new_down) / Y_this > self.Y_res:
            # Add new Y value below
            vals += [Y_new_down]
        # end if
        return vals
    # end def

# end class
