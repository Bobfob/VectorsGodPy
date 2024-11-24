'''
The Global Scope module contains many useful methods and functions. (By default, Global Scope is imported into the project automatically)
\nWarning! This module contains almost all functions and methods from math module!

Examples:
>>> # the lerp function (i.e. linear interpolation function)
>>> # linearly interpolates from a to b with a given weight t
>>> print(lerp(0.5, 1, 0.5))
>>> 0.75

>>>
>>> # the is_approx_zero function checks when some value
>>> # is approximately equals zero with given epsilon (i.e. coefficient of approximation)

>>>
>>> # if value is below the epsilon than this method returns true, otherwise returns false
>>> print(is_approx_zero(0.001, 1e-3))
>>> True
'''

from math import (
    cos, sin, tan, cosh, sinh, tanh, atan2,
    acos, asin, asinh, acosh, atanh, exp,
    log, log10, log2, floor, ceil, fmod,
    hypot, pow, sqrt, isinf, isnan, e, pi, tau
)

import time

# Class LCG for random number generator
class LCG:
  def __init__(self, seed):
    self.m = 2**32
    self.a = 1103515245
    self.c = 12345
    self.state = seed

  def next(self):
    self.state = (self.a * self.state + self.c) % self.m
    return self.state

# Definition of methods

def is_approx_zero(value: float, eps: float=1e-8) -> bool:
    '''Checks when some value is approximately equals zero with given epsilon (i.e. coefficient of approximation)'''
    return value < eps

def is_approx_equal(this: float, other: float, eps: float=1e-8) -> bool:
    '''Checks when some two values are approximately equals each other with given epsilon (i.e. coefficient of approximation)'''
    return (other - this) < eps

def lerp(a: float, b: float, t: float) -> float:
    '''linearly interpolates from a to b with a given weight t'''
    return a + (b - a) * t

def clamp(x: int | float, x_min: int | float, x_max: int | float) -> int | float:
    '''Clamps value x between x_min and x_max'''
    return min(max(x_min, x), x_max)

def clampi(x: int, x_min: int, x_max: int) -> int:
    '''Clamps value x between x_min and x_max (only integers)'''
    if isinstance(x, int) and isinstance(x_min, int) and isinstance(x_max, int):
        return min(max(x_min, x), x_max)
    else:
        raise TypeError("This function works only with integer values!")

def clampf(x: float, x_min: float, x_max: float) -> float:
    '''Clamps value x between x_min and x_max (only floats)'''
    if isinstance(x, float) and isinstance(x_min, float) and isinstance(x_max, float):
        return min(max(x_min, x), x_max)
    else:
        raise TypeError("This function only works with float values!")

def inverse_lerp(from_: float, to: float, weight: float) -> float:
    '''Returns an interpolation or extrapolation factor.'''
    if to == from_:
        return 0.0 # Avoid division by zero
    return (weight - from_) / (to - from_)

def is_finite(x) -> bool:
    '''Checks if a value can be converted to a finite number. Handles various input types.'''
    try:
        num = float(x)
        if isinf(num) or isnan(num):
            return False
        else:
            return True
    except ValueError:
        return False
    except TypeError:
        return False

def deg_to_rad(value: int | float) -> float:
    return value / 180 * pi

def rad_to_deg(value: float) -> float:
    return value / pi * 180

def lerp_angle(this: float, to: float, weight: float) -> float:
    '''Linearly interpolates between two angles (in radians).'''
    # Handle angle wrapping
    diff = to - this
    shortest_diff = atan2(sin(diff), cos(diff)) # handles angles > 180 degrees and < -180 degrees

    return this + shortest_diff * weight

def randf() -> float:
    '''Generates a pseudo-random float between 0.0 and 1.0 (inclusive) using LCG and nanosecond seed.'''
    nanoseconds = int(time.perf_counter_ns())
    generator = LCG(nanoseconds)
    return generator.next() / generator.m
