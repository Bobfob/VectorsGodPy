from VectorsGodPy.globalScope import *
from typing import overload, Union, TypeVar

Scalar = int | float
Scalar_or_Vector3 = Union[Scalar, "Vector3"]
Int_or_Vector3i = Union[int, "Vector3i"]
MaxInt32 = 2**31 - 1
MinInt32 = -2**31

class Vector3i:
    @overload
    def __init__(self) -> "Vector3i":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x: int) -> "Vector3i":
        '''Constructs Vector3i with x component'''
        ...
    @overload
    def __init__(self, x: int, y: int, z: int) -> "Vector3i":
        '''Constructs Vector3i with x, y and z components'''
        ...
    @overload
    def __init__(self, vector: "Vector3i") -> "Vector3i":
        '''Constructs Vector3i from other Vector3i'''
        ...
    @overload
    def __init__(self, vector: "Vector3") -> "Vector3i":
        '''Constructs Vector3i from other Vector3'''
        ...
    
    def __init__(self, *args):
        if len(args) == 0: # First constructor
            self.x = 0
            self.y = 0
            self.z = 0
        elif len(args) == 1 and isinstance(args[0], int): # Second constructor
            self.x = args[0]
            self.y = args[0]
            self.z = args[0]
        elif len(args) == 3 and isinstance(args[0], int) and isinstance(args[1], int) and isinstance(args[2], int): # Third constructor
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        elif len(args) == 1 and isinstance(args[0], Vector3i): # Fourth constructor
            self.x = args[0].x
            self.y = args[0].y
            self.z = args[0].z
        elif len(args) == 1 and isinstance(args[0], Vector3): # Fifth constructor
            self.x = int(args[0].x)
            self.y = int(args[0].y)
            self.z = int(args[0].z)
        else:
            raise TypeError("Invalid Argument")
    
    def INF() -> "Vector3i":
        return Vector3i(MaxInt32)
    def M_INF() -> "Vector3i":
        return Vector3i(MinInt32)
    
    def ZERO() -> "Vector3i":
        return Vector3i()
    def ONE() -> "Vector3i":
        return Vector3i(1)
    
    def FORWARD() -> "Vector3i":
        return Vector3i(0, 0, -1)
    def BACKWARD() -> "Vector3i":
        return Vector3i(0, 0, 1)
    def LEFT() -> "Vector3i":
        return Vector3i(-1, 0, 0)
    def RIGHT() -> "Vector3i":
        return Vector3i(1, 0, 0)
    def UP() -> "Vector3i":
        return Vector3i(0, 1, 0)
    def DOWN() -> "Vector3i":
        return Vector3i(0, -1, 0)
    
    def __str__(self) -> str:
        return f'({self.x}, {self.y}, {self.z})'
    
    def __repr__(self) -> str:
        return f'Vector3i({self.x}, {self.y}, {self.z})'

    def __add__(self, other : "Vector3i"):
        if isinstance(other, Vector3i):
            return Vector3i(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            raise TypeError("Invalid Argument")
    
    def __sub__(self, other: "Vector3i"):
        if isinstance(other, Vector3i):
            return Vector3i(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            raise TypeError("Invalid Argument")

    def __mul__(self, other: Int_or_Vector3i):
        if isinstance(other, int):
            return Vector3i(self.x * other, self.y * other, self.z * other)
        elif isinstance(other, Vector3i):
            return Vector3i(self.x * other.x, self.y * other.y, self.z * other.z)
        else:
            raise TypeError("Invalid Argument")
    
    def __rmul__(self, other: Int_or_Vector3i):
        return self.__mul__(other)
    
    def __truediv__(self, other: Int_or_Vector3i):
        if isinstance(other, int):
            return Vector3i(self.x // other, self.y // other, self.z // other)
        elif isinstance(other, Vector3i):
            return Vector3i(self.x // other.x, self.y // other.y, self.z // other.z)
        else:
            raise TypeError("Invalid Argument")
        
    def __floordiv__(self, other: Int_or_Vector3i):
        if isinstance(other, int):
            return Vector3i(self.x // other, self.y // other, self.z // other)
        elif isinstance(other, Vector3i):
            return Vector3i(self.x // other.x, self.y // other.y, self.z // other.z)
        else:
            raise TypeError("Invalid Argument")
    
    def __abs__(self) -> "Vector3i":
        return Vector3i(abs(self.x), abs(self.y), abs(self.z))
    
    def __min__(self, other: "Vector3i") -> "Vector3i":
        return Vector3i(min(self.x, other.x), min(self.y, other.y), min(self.z, other.z))

    def __max__(self, other: "Vector3i") -> "Vector3i":
        return Vector3i(max(self.x, other.x), max(self.y, other.y), max(self.z, other.z))

    def clamp(self, x_min: "Vector3i", x_max: "Vector3i") -> "Vector3i":
        if isinstance(x_min, Vector3i) and isinstance(x_max, Vector3i):
            return min(max(x_min, self), x_max)
    
    def distance_to(self, to: "Vector3i") -> float:
        if isinstance(to, Vector3i):
            return ((to.x - self.x)**2 + (to.y - self.y)**2 + (to.z - self.z)**2)**0.5
        else:
            raise TypeError("Invalid Argument")
    
    def distance_squared_to(self, to: "Vector3i") -> int:
        if isinstance(to, Vector3i):
            return ((to.x - self.x)**2 + (to.y - self.y)**2 + (to.z - self.z)**2)
        else:
            raise TypeError("Invalid Argument")
    
    def length(self) -> float:
        return (self.x*self.x + self.y*self.y + self.z*self.z)**0.5
    
    def length_squared(self) -> int:
        return self.x*self.x + self.y*self.y + self.z*self.z
    
    def __min__(self, other: "Vector3i") -> "Vector3i":
        if isinstance(other, Vector3i):
            return Vector3i(min(self.x, other.x), min(self.y, other.y), min(self.z, other.z))
        else:
            raise TypeError("Invalid Argument")
    
    def __max__(self, other: "Vector3i") -> "Vector3i":
        if isinstance(other, Vector3i):
            return Vector3i(max(self.x, other.x), max(self.y, other.y), max(self.z, other.z))
        else:
            raise TypeError("Invalid Argument")
        
    def sign(self) -> "Vector3i":
        return self / abs(self)

class Vector3:
    @overload
    def __init__(self) -> "Vector3":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x: Scalar) -> "Vector3":
        '''Constructs Vector3 with x component'''
        ...
    @overload
    def __init__(self, x: Scalar, y: Scalar, z: Scalar) -> "Vector3":
        '''Constructs Vector3 with x, y and z components'''
        ...
    @overload
    def __init__(self, vector: "Vector3") -> "Vector3":
        '''Constructs Vector3 from other Vector3'''
        ...
    @overload
    def __init__(self, vector: "Vector3i") -> "Vector3":
        '''Constructs Vector3 from other Vector3i'''
        ...
    
    def __init__(self, *args) -> "Vector3":
        if len(args) == 0: # First constructor
            self.x = 0
            self.y = 0
            self.z = 0
        elif len(args) == 1 and isinstance(args[0], (int, float)): # Second constructor
            self.x = args[0]
            self.y = args[0]
            self.z = args[0]
        elif len(args) == 3 and isinstance(args[0], (int, float)) and isinstance(args[1], (int, float)) and isinstance(args[2], (int, float)): # Third constructor
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
        elif len(args) == 1 and isinstance(args[0], Vector3): # Fourth constructor
            self.x = args[0].x
            self.y = args[0].y
            self.z = args[0].z
        elif len(args) == 1 and isinstance(args[0], Vector3i): # Fifth constructor
            self.x = float(args[0].x)
            self.y = float(args[0].y)
            self.z = float(args[0].z)
        else:
            raise TypeError("Invalid Argument")

    def INF() -> "Vector3":
        return Vector3(float("inf"))
    def M_INF() -> "Vector3":
        return Vector3(float('-inf'))

    def ZERO() -> "Vector3":
        return Vector3()
    def ONE() -> "Vector3":
        return Vector3(1)

    def FORWARD() -> "Vector3":
        return Vector3(0, 0, -1)
    def BACKWARD() -> "Vector3":
        return Vector3(0, 0, 1)
    def LEFT() -> "Vector3":
        return Vector3(-1, 0, 0)
    def RIGHT() -> "Vector3":
        return Vector3(1, 0, 0)
    def UP() -> "Vector3":
        return Vector3(0, 1, 0)
    def DOWN() -> "Vector3":
        return Vector3(0, -1, 0)

    # Main Methods

    def __str__(self) -> str:
        return f"({self.x}, {self.y}, {self.z})"

    def __repr__(self) -> str:
        return f"Vector3({self.x}, {self.y}, {self.z})"

    def __add__(self, other: "Vector3"):
        if isinstance(other, Vector3):
            return Vector3(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            raise TypeError("Invalid Argument")

    def __sub__(self, other: "Vector3"):
        if isinstance(other, Vector3):
            return Vector3(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            raise TypeError("Invalid Argument")

    def __mul__(self, other : Union[int, float, "Vector3"]):
        if isinstance(other, (int, float)):
            return Vector3(
                self.x * other, self.y * other, self.z * other
            )
        elif isinstance(other, Vector3):
            return Vector3(
                self.x * other.x, self.y * other.y, self.z * other.z
            )
        else:
            raise TypeError("Invalid Argument")

    def __rmul__(self, other : Union[int, float, "Vector3"]):
        return self.__mul__(other)

    def __floordiv__(self, other : Union[int, float, "Vector3"]):
        if isinstance(other, (int, float)):
            return Vector3(
                self.x // other, self.y // other, self.z // other
            )
        elif isinstance(other, Vector3):
            return Vector3(
                self.x // other.x, self.y // other.y, self.z // other.z
            )
        else:
            raise TypeError("Invalid Argument")

    def __truediv__(self, other : Union[int, float, "Vector3"]):
        if isinstance(other, (int, float)):
            return Vector3(
                self.x / other, self.y / other, self.z / other
            )
        elif isinstance(other, Vector3):
            return Vector3(
                self.x / other.x, self.y / other.y, self.z / other.z
            )
        else:
            raise TypeError("Invalid Argument")

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError("Index out of range!")

    def __neg__(self):
        return Vector3(-self.x, -self.y, -self.z)

    def __ne__(self, b: "Vector3"):
        return(self.x != b.x) and (self.y != b.y) and (self.z != b.z)

    def __eq__(self, b: "Vector3"):
        return (self.x == b.x) and (self.y == b.y) and (self.z == b.z)

    def __lt__(self, b: "Vector3"):
        return (self.x < b.x) and (self.y < b.y) and (self.z < b.z)

    def __gt__(self, b: "Vector3"):
        return (self.x > b.x) and (self.y > b.y) and (self.z > b.z)

    def __lq__(self, b: "Vector3"):
        return (self.x <= b.x) and (self.y <= b.y) and (self.z <= b.z)

    def __gq__(self, b: "Vector3"):
        return (self.x >= b.x) and (self.y >= b.y) and (self.z >= b.z)

    def __floor__(self):
        return Vector3(int(self.x), int(self.y), int(self.z))

    def __round__(self, ndigits=None):
        return Vector3(round(self.x, ndigits), round(self.y, ndigits), round(self.z, ndigits))

    def __ceil__(self):
        return Vector3(ceil(self.x), ceil(self.y), ceil(self.z))

    def __abs__(self) -> "Vector3":
        return Vector3(abs(self.x), abs(self.y), abs(self.z))

    def length(self) -> float:
        return (self.x**2 + self.y**2 + self.z**2) ** 0.5

    def length_squared(self) -> float:
        return self.x**2 + self.y**2 + self.z**2

    def normalize(self) -> "Vector3":
        if self.length() == 0:
            return Vector3(0, 0, 0)
        return Vector3(
            self.x / self.length(), self.y / self.length(), self.z / self.length()
        )

    def dot(self, vec: "Vector3") -> float:
        return self.x * vec.x + self.y * vec.y + self.z * vec.z

    def clamp(self, x_min: "Vector3", x_max: "Vector3") -> "Vector3":
        return Vector3(
            clamp(self.x, x_min.x, x_max.x),
            clamp(self.y, x_min.y, x_max.y),
            clamp(self.z, x_min.z, x_max.z)
        )

    def angle_to(self, vec: "Vector3") -> "Vector3":
        norm_self = self.normalize()
        norm_to = vec.normalize()

        dot = norm_self.dot(norm_to)

        dot = clamp(dot, -1, 1)

        return acos(dot)

    def bezier_derivative(
        self, control_1: "Vector3", control_2: "Vector3", end: "Vector3", t: float
    ) -> "Vector3":
        return (
            3 * (control_1 - self) * (1 - t) * (1 - t)
            + 6 * (control_2 - control_1) * (1 - t) * t
            + 3 * (end - control_2) * t * t
        )

    def bezier_interpolate(
        self, control_1: "Vector3", control_2: "Vector3", end: "Vector3", t: float
    ) -> "Vector3":
        u: float = 1 - t
        return (
            u * u * u * self
            + 3 * u * u * t * control_1
            + 3 * u * t * t * control_2
            + t * t * t * end
        )

    def bounce(self, n: "Vector3") -> "Vector3":
        return self - 2 * self.dot(n.normalize()) * n

    def cross(self, vec: "Vector3") -> "Vector3":
        cross_x = self.y * vec.z - self.z * vec.y
        cross_y = self.z * vec.x - self.x * vec.z
        cross_z = self.x * vec.y - self.y * vec.x

        return Vector3(cross_x, cross_y, cross_z)

    def direction_to(self, vec: "Vector3") -> "Vector3":
        return (vec - self).normalize()

    def inverse(self) -> "Vector3":
        return 1 / self

    def distance_to(self, vec: "Vector3") -> "Vector3":
        return (vec - self).length()

    def distance_to_squared(self, vec: "Vector3") -> "Vector3":
        return (vec - self).length_squared()

    def lerp(self, vec: "Vector3", t: float) -> "Vector3":
        return self + (vec - self) * t

    def min_axis_index(self) -> int:
        if (self.x == self.y and self.y == self.z and self.z == self.x) or (min((self.x, self.y, self.z)) == self.x):
            return 0
        elif min((self.x, self.y, self.z)) == self.y:
            return 1
        else:
            return 2

    def max_axis_index(self) -> int:
        if (self.x == self.y and self.y == self.z and self.z == self.x) or (max((self.x, self.y, self.z)) == self.x):
            return 0
        elif max((self.x, self.y, self.z)) == self.y:
            return 1
        else:
            return 2

    def __min__(self, vec: "Vector3") -> "Vector3":
        return Vector3(min(self.x, vec.x), min(self.y, vec.y), min(self.z, vec.z))

    def __max__(self, vec: "Vector3") -> "Vector3":
        return Vector3(max(self.x, vec.x), max(self.y, vec.y), max(self.z, vec.z))

    def rotated(self, axis: "Vector3", angle: float) -> "Vector3":
        axis = axis.normalize()
        ux, uy, uz = axis.x, axis.y, axis.z
        cos_angle = cos(angle)
        sin_angle = sin(angle)

        rotated_vector = Vector3(
            (cos_angle + ux * ux * (1 - cos_angle)) * self.x
            + (ux * uy * (1 - cos_angle) - uz * sin_angle) * self.y
            + (ux * uz * (1 - cos_angle) + uy * sin_angle) * self.z,
            (uy * ux * (1 - cos_angle) + uz * sin_angle) * self.x
            + (cos_angle + uy * uy * (1 - cos_angle)) * self.y
            + (uy * uz * (1 - cos_angle) - ux * sin_angle) * self.z,
            (uz * ux * (1 - cos_angle) - uy * sin_angle) * self.x
            + (uz * uy * (1 - cos_angle) + ux * sin_angle) * self.y
            + (cos_angle + uz * uz * (1 - cos_angle)) * self.z,
        )

        return rotated_vector

    def sign(self) -> "Vector3":
        """Get sign of all components"""
        return self / abs(self)
