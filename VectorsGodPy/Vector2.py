from typing import overload, Union
from VectorsGodPy.globalScope import *


Scalar = Union[int, float]
MAX_INT32 = 2**31 - 1
MIN_INT32 = -2**31

class _Vector2i:
    '''This class is a dummy type! Don't use it!'''
    ...

class _Vector2:
    '''This class is a dummy type! Don't use it!'''
    ...

class Vector2i:
    # Overloaded constructors
    @overload
    def __init__(self) -> "Vector2i":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x : int) -> "Vector2i":
        '''Constructs Vector2i with x component'''
        ...
    @overload
    def __init__(self, x : int, y : int) -> "Vector2i":
        '''Constructs Vector2i with x and y components'''
        ...
    @overload
    def __init__(self, vector : "Vector2i") -> "Vector2i":
        '''Constructs Vector2i from other Vector2i'''
        ...
    @overload
    def __init__(self, vector : "Vector2") -> "Vector2i":
        '''Constructs Vector2i from other Vector2'''
        ...
    
    # Main init method
    def __init__(self, *args):
        if args == ():
            self.x = 0
            self.y = 0
        elif len(args) == 1 and isinstance(args[0], int):
            self.x = args[0]
            self.y = args[0]
        elif isinstance(args[0], int) and isinstance(args[1], int):
            self.x = args[0]
            self.y = args[1]
        elif isinstance(args[0], Vector2):
            self.x = int(args[0].x)
            self.y = int(args[0].y)
        elif isinstance(args[0], Vector2i):
            self.x = int(args[0].x)
            self.y = int(args[0].y)
        else:
            raise TypeError("Invalid Argument")
    
    # Constants
    def ZERO() -> "Vector2i":
        return Vector2i()
    def ONE() -> "Vector2i":
        return Vector2i(1)
    def M_INF() -> "Vector2i":
        return Vector2i(MAX_INT32)
    def INF() -> "Vector2":
        return Vector2i(MIN_INT32)
    def LEFT() -> "Vector2i":
        return Vector2i(-1, 0)
    def RIGHT() -> "Vector2i":
        return Vector2i(1, 0)
    def UP() -> "Vector2i":
        return Vector2i(0, -1)
    def DOWN() -> "Vector2i":
        return Vector2i(0, 1)
    
    # Representation
    def __repr__(self) -> str:
        return f"Vector2i({self.x}, {self.y})"
    
    # To string conversion
    def __str__(self) -> str:
        return f"({self.x}, {self.y})"
    
    # Operators dander methods
    def __lt__(self, other : "Vector2i") -> bool:
        return self.x < other.x and self.y < other.y
    
    def __gt__(self, other : "Vector2i") -> bool:
        return self.x > other.x and self.y > other.y
    
    def __lq__(self, other : "Vector2i") -> bool:
        return self.x <= other.x and self.y <= other.y
    
    def __gq__(self, other : "Vector2i") -> bool:
        return self.x >= other.x and self.y >= other.y
    
    def __eq__(self, other : "Vector2i") -> bool:
        return self.x == other.x and self.y == other.y
    
    def __ne__(self, other : "Vector2i") -> bool:
        return self.x != other.x and self.y != other.y
    
    def __neg__(self) -> "Vector2i":
        return Vector2(-self.x, -self.y)
    
    def __add__(self, other : "Vector2i") -> "Vector2i":
        if isinstance(other, Vector2i):
            return Vector2i(self.x + other.x, self.y + other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __sub__(self, other : "Vector2i") -> "Vector2i":
        if isinstance(other, Vector2i):
            return Vector2i(self.x - other.x, self.y - other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __mul__(self, other : int | float | _Vector2i) -> "Vector2i":
        if isinstance(other, int):
            return Vector2i(self.x * other, self.y * other)
        elif isinstance(other, Vector2i):
            return Vector2i(self.x * other.x, self.y * other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __rmul__(self, other : int | float | _Vector2i) -> "Vector2i":
        return self.__mul__(other)
    
    def __truediv__(self, other : int | float | _Vector2i) -> "Vector2i":
        if isinstance(other, int):
            return Vector2i(self.x // other, self.y // other)
        elif isinstance(other, Vector2i):
            return Vector2i(self.x // other.x, self.y // other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __floordiv__(self, other : int | _Vector2i) -> "Vector2i":
        if isinstance(other, int):
            return Vector2i(self.x // other, self.y // other)
        elif isinstance(other, Vector2i):
            return Vector2i(self.x // other.x, self.y // other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __mod__(self, other : int | _Vector2i) -> "Vector2i":
        if isinstance(other, int):
            return Vector2i(self.x % other, self.y % other)
        elif isinstance(other, Vector2i):
            return Vector2i(self.x % other.x, self.y % other.y)
        else:
            raise TypeError("Invalid Argument")

    def __getindex__(self, index : int) -> int:
        match (index):
            case 0:
                return self.x
            case 1:
                return self.y
            case _:
                raise IndexError("Index is out of range!")
    
    # Some useful dander/normal methods
    def __abs__(self) -> "Vector2i":
        return Vector2i(abs(self.x), abs(self.y))
    
    def is_approx_equal(self, other : "Vector2i", eps : float=1e-8) -> bool:
        '''Returns True if two vectors are approximately equal'''
        if isinstance(other, Vector2i):
            return (other - self) < eps
        else:
            raise TypeError("Invalid Argument")
    
    def is_approx_zero(self, eps : float=1e-8) -> bool:
        '''Returns True if vector is approximately equals zero'''
        return self.x < eps and self.y < eps
    
    def clamp(self, x_min : "Vector2i", x_max : "Vector2i") -> "Vector2i":
        '''Returns a new vector with all components clamped between the components of <b>min</b> and <b>max</b>'''
        if isinstance(x_min, Vector2i) and isinstance(x_max, Vector2i):
            return Vector2i(min(max(x_min.x, self.x), x_max.x), min(max(x_min.y, self.y), x_max.y))
        else:
            raise TypeError("Invalid Argument")
    
    def clampi(self, x_min : int, x_max : int) -> "Vector2":
        '''Returns a new vector with all components clamped between the <b>min</b> and <b>max</b>'''
        if isinstance(x_min, int) and isinstance(x_max, int):
            return Vector2(min(max(x_min, self.x), x_max), min(max(x_min, self.y), x_max))
        else:
            raise TypeError("Invalid argument")
    
    # Main methods for Vector2i in 2D space
    def length(self) -> float:
        '''Returns magnitude of <b>this</b> vector'''
        return (self.x ** 2 + self.y ** 2) ** 0.5
    
    def length_squared(self) -> int:
        '''Returns squared magnitude of <b>this</b> vector'''

    def distance_to(self, other : "Vector2i") -> float:
        '''Returns distance from <b>this</b> vector to <b>other</b> vector'''
        if isinstance(other, Vector2i):
            return (other - self).length()
        else:
            raise TypeError("Invalid Argument")
        
    def distance_to_squared(self, other : "Vector2i") -> int:
        '''Returns distance from <b>this</b> vector to <b>other</b> vector'''
        if isinstance(other, Vector2i):
            return (other - self).length_squared()
        else:
            raise TypeError("Invalid Argument")
    
    def sign(self) -> "Vector2i":
        '''Returns <b>signed</b> vector of <b>this</b> vector'''
        return self / abs(self)

class Vector2:
    # Overloaded constructors
    @overload
    def __init__(self) -> "Vector2":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x : Scalar) -> "Vector2":
        '''Constructs Vector2 with one component'''
        ...
    @overload
    def __init__(self, x : Scalar, y : Scalar) -> "Vector2":
        '''Constructs Vector2 with x and y components'''
        ...
    @overload
    def __init__(self, vector : "Vector2") -> "Vector2":
        '''Constructs Vector2 from other Vector2'''
        ...
    @overload
    def __init__(self, vector : "Vector2i") -> "Vector2":
        '''Constructs Vector2 from other Vector2i'''
        ...

    # Main init method
    def __init__(self, *args):
        if args == ():
            self.x = 0
            self.y = 0
        elif len(args) == 1 and isinstance(args[0], (int, float)):
            self.x = args[0]
            self.y = args[0]
        elif isinstance(args[0], (int, float)) and isinstance(args[1], (int, float)):
            self.x = args[0]
            self.y = args[1]
        elif len(args) == 1 and isinstance(args[0], Vector2):
            self.x = args[0].x
            self.y = args[0].y
            return
        elif len(args) == 1 and isinstance(args[0], Vector2i):
            self.x = args[0].x
            self.y = args[0].y
        else:
            raise TypeError("Invalid Argument")

    # Constants
    def ZERO() -> "Vector2":
        return Vector2()
    def ONE() -> "Vector2":
        return Vector2(1)
    def M_INF() -> "Vector2":
        return Vector2(float('-inf'))
    def INF() -> "Vector2":
        return Vector2(float('inf'))
    def LEFT() -> "Vector2":
        return Vector2(-1, 0)
    def RIGHT() -> "Vector2":
        return Vector2(1, 0)
    def UP() -> "Vector2":
        return Vector2(0, -1)
    def DOWN() -> "Vector2":
        return Vector2(0, 1)

    # Representation
    def __repr__(self) -> str:
        return f"Vector2({self.x}, {self.y})"
    
    # To string conversion
    def __str__(self) -> str:
        return f"({self.x}, {self.y})"
    
    # Operators dander methods
    def __lt__(self, other : "Vector2") -> bool:
        return self.x < other.x and self.y < other.y
    
    def __gt__(self, other : "Vector2") -> bool:
        return self.x > other.x and self.y > other.y
    
    def __lq__(self, other : "Vector2") -> bool:
        return self.x <= other.x and self.y <= other.y
    
    def __gq__(self, other : "Vector2") -> bool:
        return self.x >= other.x and self.y >= other.y
    
    def __eq__(self, other : "Vector2") -> bool:
        return self.x == other.x and self.y == other.y
    
    def __ne__(self, other : "Vector2") -> bool:
        return self.x != other.x and self.y != other.y
    
    def __neg__(self) -> "Vector2":
        return Vector2(-self.x, -self.y)
    
    def __add__(self, other : "Vector2") -> "Vector2":
        if isinstance(other, Vector2):
            return Vector2(self.x + other.x, self.y + other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __sub__(self, other : "Vector2") -> "Vector2":
        if isinstance(other, Vector2):
            return Vector2(self.x - other.x, self.y - other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __mul__(self, other : int | float | _Vector2) -> "Vector2":
        if isinstance(other, (int, float)):
            return Vector2(self.x * other, self.y * other)
        elif isinstance(other, Vector2):
            return Vector2(self.x * other.x, self.y * other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __rmul__(self, other) -> "Vector2":
        return self.__mul__(other)
    
    def __truediv__(self, other : int | float | _Vector2) -> "Vector2":
        if isinstance(other, (int, float)):
            return Vector2(self.x / other, self.y / other)
        elif isinstance(other, Vector2):
            return Vector2(self.x / other.x, self.y / other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __floordiv__(self, other : int | float | _Vector2) -> "Vector2":
        if isinstance(other, (int, float)):
            return Vector2(self.x // other, self.y // other)
        elif isinstance(other, Vector2):
            return Vector2(self.x // other.x, self.y // other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __mod__(self, other : int | _Vector2 | _Vector2i) -> "Vector2":
        if isinstance(other, (int, float)):
            return Vector2(self.x % other, self.y % other)
        elif isinstance(other, (Vector2, Vector2i)):
            return Vector2(self.x % other.x, self.y % other.y)
        else:
            raise TypeError("Invalid Argument")
    
    def __getindex__(self, index : int) -> int | float:
        match (index):
            case 0:
                return self.x
            case 1:
                return self.y
            case _:
                raise IndexError("Index is out of range!")
    
    # Some useful dander/normal methods
    def __abs__(self) -> "Vector2":
        return Vector2(abs(self.x), abs(self.y))
    
    def __floor__(self) -> "Vector2":
        return Vector2(floor(self.x), floor(self.y))
    
    def __round__(self, n_digits=None) -> "Vector2":
        return Vector2(round(self.x, n_digits), round(self.y, n_digits))
    
    def __ceil__(self) -> "Vector2":
        return Vector2(ceil(self.x), ceil(self.y))
    
    def __int__(self) -> "Vector2":
        return Vector2(int(self.x), int(self.y))
    
    def is_approx_equal(self, other : "Vector2", eps : float=1e-8) -> bool:
        '''Returns True if two vectors are approximately equal'''
        return (other - self) < eps
    
    def is_approx_zero(self, eps : float=1e-8) -> bool:
        '''Returns True if vector is approximately equals zero'''
        return self.x < eps and self.y < eps
    
    def clamp(self, x_min : "Vector2", x_max : "Vector2") -> "Vector2":
        '''Returns a new vector with all components clamped between the components of <b>min</b> and <b>max</b>'''
        return Vector2(min(max(x_min.x, self.x), x_max.x), min(max(x_min.y, self.y), x_max.y))
    
    def clampf(self, x_min : float, x_max : float) -> "Vector2":
        '''Returns a new vector with all components clamped between the <b>min</b> and <b>max</b>'''
        return Vector2(min(max(x_min, self.x), x_max), min(max(x_min, self.y), x_max))

    # Main methods for Vector2 in 2D space
    def length(self) -> float:
        '''Returns magnitude of <b>this</b> vector'''
        return (self.x ** 2 + self.y ** 2) ** 0.5
    
    def length_squared(self) -> int | float:
        '''Returns squared magnitude of <b>this</b> vector'''
        return self.x ** 2 + self.y ** 2
    
    def angle(self, as_degrees=False) -> float:
        '''Returns angle of <b>this</b> vector in radians (or degrees if needed)'''
        if not as_degrees:
            return atan2(self.y, self.x) % 2*pi
        else:
            return (atan2(self.y, self.x) * 180 / pi) % 360
        
    def normalized(self) -> "Vector2":
        '''Returns normalized vector'''
        return self / self.length()
    
    def rotate(self, angle : float, as_degrees=False) -> "Vector2":
        '''Returns rotated vector by some angle in radians (or degrees if needed)'''
        match (as_degrees):
            case False:
                if angle == pi/2 or angle == -3*pi/2 and not as_degrees:
                    return Vector2(-self.y, self.x)
                elif angle == pi or angle == -pi and not as_degrees:
                    return Vector2(-self.x, -self.y)
                elif angle == 3*pi/2 or angle == -pi/2 and not as_degrees:
                    return Vector2(self.y, -self.x)
                
                return Vector2(self.x * cos(angle) - self.y * sin(angle), self.x * sin(angle) + self.y * cos(angle))
            
            case True:
                if angle == 90 or angle == -270 and not as_degrees:
                    return Vector2(-self.y, self.x)
                elif angle == 180 or angle == -180 and not as_degrees:
                    return Vector2(-self.x, -self.y)
                elif angle == 270 or angle == -90 and not as_degrees:
                    return Vector2(self.y, -self.x)
                
                return Vector2(self.x * cos(angle*pi/180) - self.y * sin(angle*pi/180), self.x * sin(angle*pi/180) + self.y * cos(angle*pi/180))
    
    def angle_to(self, other : "Vector2", as_degrees=False) -> float:
        '''
        Returns angle to <b>other</b> vector in radians (or degrees if needed).
         <b>a.angle_to(b)</b> is equivalent to <b>(other - self).angle()</b>
        '''
        return (other - self).angle(as_degrees)
    
    def direction_to(self, other : "Vector2") -> "Vector2":
        '''
        Returns direction to <b>other</b> vector.
         <b>a.direction_to(b)</b> is equivalent to <b>(other - self).normalized()</b>
        '''
        return (other - self).normalized()
    
    def dot(self, other : "Vector2") -> float:
        '''Returns the result of the dot product of <b>this</b> vector and <b>other</b> vector'''
        return self.x * other.x + self.y * other.y
    
    def dist_to(self, other : "Vector2") -> float:
        '''Returns distance from <b>this</b> vector to <b>other</b> vector'''
        return (other - self).length()
    
    def dist_to_squared(self, other : "Vector2") -> int | float:
        '''Returns squared distance from <b>this</b> vector to <b>other</b> vector'''
        return (other - self).length_squared()
    
    def lerp(self, other : "Vector2", t : float) -> "Vector2":
        '''Returns the result of linear interpolation from <b>this</b> vector to <b>other</b> vector by given weight <b>t</b>'''
        return self + (other - self) * t
    
    def is_normalized(self, eps : float=1e-8) -> bool:
        '''Returns True if vector is normalized (uses approximation for more accurate result)'''
        return abs(self.length() - 1.0) < eps
    
    def sign(self) -> "Vector2":
        '''Returns <b>signed</b> vector of <b>this</b> vector'''
        return int(self) / int(abs(self))
    
    def aspect(self) -> "Vector2":
        '''Returns aspect ratio of <b>this</b> vector'''
        return self.x / self.y
    
    def bezier_derivative(self, control_1 : "Vector2", control_2 : "Vector2", end : "Vector2", t : float) -> "Vector2":
        '''Calculates the derivative of a Cubic Bezier curve at parameter t.'''
        return (
            3 * (1 - t) * (1 - t) * (control_1 - self)
            + 6 * (1 - t) * t * (control_2 - control_1)
            + 3 * t * t * (end - control_2)
        )
    
    def bezier_interpolate(self, control_1 : "Vector2", control_2 : "Vector2", end : "Vector2", t : float) -> "Vector2":
        '''Calculates the Cubic Bezier curve at parameter t'''
        q0 : Vector2 = self.lerp(control_1, t)
        q1 : Vector2 = control_1.lerp(control_2, t)
        q2 : Vector2 = control_2.lerp(end, t)

        r0 : Vector2 = q0.lerp(q1, t)
        r1 : Vector2 = q1.lerp(q2, t)

        return r0.lerp(r1, t)
