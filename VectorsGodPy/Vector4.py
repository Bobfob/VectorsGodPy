from VectorsGodPy.globalScope import *
from typing import overload, Union

Scalar = int | float
Scalar_or_Vector3 = Union[Scalar, "Vector4"]
Int_or_Vector3i = Union[int, "Vector4i"]
MaxInt32 = 2**31 - 1
MinInt32 = -2**31

class Vector4i:
    @overload
    def __init__(self) -> "Vector4i":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x: int) -> "Vector4i":
        '''Constructs Vector4i with x component'''
        ...
    @overload
    def __init__(self, x: int, y: int, z: int, w: int) -> "Vector4i":
        '''Constructs Vector4i with x, y, z and w components'''
        ...
    @overload
    def __init__(self, vector: "Vector4i") -> "Vector4i":
        '''Constructs Vector4i from other Vector4i'''
        ...
    @overload
    def __init__(self, vector: "Vector4") -> "Vector4i":
        '''Constructs Vector4i from other Vector4'''
        ...
    
    def __init__(self, *args):
        if len(args) == 0:
            self.x = 0
            self.y = 0
            self.z = 0
            self.w = 0
        elif len(args) == 1 and isinstance(args[0], int):
            self.x = args[0]
            self.y = args[0]
            self.z = args[0]
            self.w = args[0]
        elif len(args) == 4 and isinstance(args[0], int) and isinstance(args[1], int) and isinstance(args[2], int) and isinstance(args[3], int):
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.w = args[3]
        elif len(args) == 1 and isinstance(args[0], Vector4i):
            self.x = args[0].x
            self.y = args[0].y
            self.z = args[0].z
            self.w = args[0].w
        elif len(args) == 1 and isinstance(args[0], Vector4):
            self.x = int(args[0].x)
            self.y = int(args[0].y)
            self.z = int(args[0].z)
            self.w = int(args[0].w)
        else:
            raise TypeError("Invalid Argument")
    
    @property
    def ZERO() -> "Vector4i":
        return Vector4i()
    @property
    def ONE() -> "Vector4i":
        return Vector4i(1)
    @property
    def M_INF() -> "Vector4i":
        return Vector4i(MinInt32)
    @property
    def INF() -> "Vector4i":
        return Vector4i(MaxInt32)

    def __str__(self) -> str:
        return f'({self.x}, {self.y}, {self.z}, {self.w})'
    
    def __repr__(self) -> str:
        return f'Vector4i({self.x}, {self.y}, {self.z}, {self.w})'

    def __add__(self, other: "Vector4i"):
        if isinstance(other, Vector4i):
            return Vector4i(self.x + other.x, self.y + other.y, self.z + other.z, self.w + other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __sub__(self, other: "Vector4i"):
        if isinstance(other, Vector4i):
            return Vector4i(self.x - other.x, self.y - other.y, self.z - other.z, self.w - other.w)

    @overload
    def __mul__(self, other: int) -> "Vector4i": ...
    @overload
    def __mul__(self, other: "Vector4i") -> "Vector4i": ...
    
    @overload
    def __rmul__(self, other: int) -> "Vector4i": ...
    @overload
    def __rmul__(self, other: "Vector4i") -> "Vector4i": ...

    @overload
    def __truediv__(self, other: int) -> "Vector4i": ...
    @overload
    def __truediv__(self, other: "Vector4i") -> "Vector4i": ...

    @overload
    def __floordiv__(self, other: int) -> "Vector4i": ...
    @overload
    def __floordiv__(self, other: "Vector4i") -> "Vector4i": ...

    @overload
    def __mod__(self, other: "Vector4i"): ...
    @overload
    def __mod__(self, other: int): ...
    
    def __mul__(self, other):
        if isinstance(other, int):
            return Vector4i(self.x * other, self.y * other, self.z * other, self.w * other)
        elif isinstance(other, Vector4i):
            return Vector4i(self.x * other.x, self.y * other.y, self.z * other.z, self.w * other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, int):
            return Vector4i(self.x // other, self.y // other, self.z // other, self.w // other)
        elif isinstance(other, Vector4i):
            return Vector4i(self.x // other.x, self.y // other.y, self.z // other.z, self.w // other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __floordiv__(self, other):
        return self.__truediv__(other)
    
    def __mod__(self, other):
        if isinstance(other, int):
            return Vector4i(self.x % other, self.y % other, self.z % other, self.w % other)
        elif isinstance(other, Vector4i):
            return Vector4i(self.x % other.x, self.y % other.y, self.z % other.z, self.w % other.w)
        else:
            raise TypeError("Invalid Argument")

    def __lt__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x < other.x and self.y < other.y and self.z < other.z and self.w < other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __gt__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x > other.x and self.y > other.y and self.z > other.z and self.w > other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __lq__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x <= other.x and self.y <= other.y and self.z <= other.z and self.w <= other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __gq__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x >= other.x and self.y >= other.y and self.z >= other.z and self.w >= other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __eq__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x == other.x and self.y == other.y and self.z == other.z and self.w == other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __ne__(self, other: "Vector4i") -> bool:
        if isinstance(other, Vector4i):
            return self.x != other.x and self.y != other.y and self.z != other.z and self.w != other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __neg__(self) -> "Vector4i":
        return Vector4i(
            -self.x,
            -self.y,
            -self.z,
            -self.w
        )
        
    def __abs__(self) -> "Vector4i":
        return Vector4i(abs(self.x), abs(self.y), abs(self.z), abs(self.w))
    
    def __getindex__(self, index: int) -> int:
        match (index):
            case 0:
                return self.x
            case 1:
                return self.y
            case 2:
                return self.z
            case 3:
                return self.w
            case _:
                raise IndexError("Index is out of range!")

    def __max__(self, other: "Vector4i") -> "Vector4i":
        if isinstance(other, Vector4i):
            return Vector4i(
                max(self.x, other.x),
                max(self.y, other.y),
                max(self.z, other.z),
                max(self.w, other.w)
            )
        else:
            raise TypeError("Invalid Argument")
    
    def __min__(self, other: "Vector4i") -> "Vector4i":
        if isinstance(other, Vector4i):
            return Vector4i(
                min(self.x, other.x),
                min(self.y, other.y),
                min(self.z, other.z),
                min(self.w, other.w)
            )
        else:
            raise TypeError("Invalid Argument")

    def clamp(self, x_min: "Vector4i", x_max: "Vector4i") -> "Vector4i":
        if isinstance(x_min, Vector4i) and isinstance(x_max, Vector4i):
            return min(max(x_min, self), x_max)
        else:
            raise TypeError("Invalid Argument")
    
    def clampi(self, x_min: int, x_max: int) -> "Vector4i":
        if isinstance(x_min, int) and isinstance(x_max, int):
            return Vector4(
                clampi(self.x, x_min, x_max),
                clampi(self.y, x_min, x_max),
                clampi(self.z, x_min, x_max),
                clampi(self.w, x_min, x_max)
            )
        else:
            raise TypeError("Invalid Argument")
    
    def distance_squared_to(self, other: "Vector4i") -> int:
        return (other.x - self.x)**2 + (other.y - self.y)**2 + (other.z - self.z)**2 + (other.w - self.w)**2
    
    def distance_to(self, other: "Vector4i") -> float:
        return self.distance_squared_to(other)**0.5

    def length(self) -> float:
        return (self.x**2 + self.y**2 + self.z**2 + self.w**2)**0.5
    
    def length_squared(self) -> int:
        return self.x**2 + self.y**2 + self.z**2 + self.w*self.w
    
    def sign(self) -> "Vector4i":
        return self / abs(self)

class Vector4:
    @overload
    def __init__(self) -> "Vector4":
        '''Default constructor with all components set to 0'''
        ...
    @overload
    def __init__(self, x: Scalar) -> "Vector4":
        '''Constructs Vector4 with x component'''
        ...
    @overload
    def __init__(self, x: Scalar, y: Scalar, z: Scalar, w: Scalar) -> "Vector4":
        '''Constructs Vector4 with x, y, z and w components'''
        ...
    @overload
    def __init__(self, vector: "Vector4") -> "Vector4":
        '''Constructs Vector4 from other Vector4'''
        ...
    @overload
    def __init__(self, vector: "Vector4i") -> "Vector4":
        '''Constructs Vector4 from other Vector4i'''
        ...
    
    def __init__(self, *args):
        if len(args) == 0:
            self.x = 0
            self.y = 0
            self.z = 0
            self.w = 0
        elif len(args) == 1 and isinstance(args[0], (int, float)):
            self.x = args[0]
            self.y = args[0]
            self.z = args[0]
            self.w = args[0]
        elif len(args) == 4 and isinstance(args[0], (int, float)) and isinstance(args[1], (int, float)) and isinstance(args[2], (int, float)) and isinstance(args[3], (int, float)):
            self.x = args[0]
            self.y = args[1]
            self.z = args[2]
            self.w = args[3]
        elif len(args) == 1 and isinstance(args[0], Vector4):
            self.x = args[0].x
            self.y = args[0].y
            self.z = args[0].z
            self.w = args[0].w
        elif len(args) == 1 and isinstance(args[0], Vector4i):
            self.x = args[0].x
            self.y = args[0].y
            self.z = args[0].z
            self.w = args[0].w
        else:
            raise TypeError("Invalid Argument")
    
    @property
    def ZERO() -> "Vector4":
        return Vector4()
    @property
    def ONE() -> "Vector4":
        return Vector4(1)
    @property
    def M_INF() -> "Vector4":
        return Vector4(float('-inf'))
    @property
    def INF() -> "Vector4":
        return Vector4(float('inf'))

    def __str__(self) -> str:
        return f'({self.x}, {self.y}, {self.z}, {self.w})'
    
    def __repr__(self) -> str:
        return f'Vector4i({self.x}, {self.y}, {self.z}, {self.w})'

    def __add__(self, other: "Vector4"):
        if isinstance(other, Vector4):
            return Vector4(self.x + other.x, self.y + other.y, self.z + other.z, self.w + other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __sub__(self, other: "Vector4"):
        if isinstance(other, Vector4):
            return Vector4(self.x - other.x, self.y - other.y, self.z - other.z, self.w - other.w)

    @overload
    def __mul__(self, other: Scalar) -> "Vector4": ...
    @overload
    def __mul__(self, other: "Vector4") -> "Vector4": ...
    
    @overload
    def __rmul__(self, other: Scalar) -> "Vector4": ...
    @overload
    def __rmul__(self, other: "Vector4") -> "Vector4": ...

    @overload
    def __truediv__(self, other: Scalar) -> "Vector4": ...
    @overload
    def __truediv__(self, other: "Vector4") -> "Vector4": ...

    @overload
    def __floordiv__(self, other: Scalar) -> "Vector4": ...
    @overload
    def __floordiv__(self, other: "Vector4") -> "Vector4": ...
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Vector4(self.x * other, self.y * other, self.z * other, self.w * other)
        elif isinstance(other, Vector4):
            return Vector4(self.x * other.x, self.y * other.y, self.z * other.z, self.w * other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Vector4(self.x / other, self.y / other, self.z / other, self.w / other)
        elif isinstance(other, Vector4):
            return Vector4(self.x / other.x, self.y / other.y, self.z / other.z, self.w / other.w)
        else:
            raise TypeError("Invalid Argument")
    
    def __floordiv__(self, other):
        return self.__truediv__(other)

    def __lt__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x < other.x and self.y < other.y and self.z < other.z and self.w < other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __gt__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x > other.x and self.y > other.y and self.z > other.z and self.w > other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __lq__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x <= other.x and self.y <= other.y and self.z <= other.z and self.w <= other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __gq__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x >= other.x and self.y >= other.y and self.z >= other.z and self.w >= other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __eq__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x == other.x and self.y == other.y and self.z == other.z and self.w == other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __ne__(self, other: "Vector4") -> bool:
        if isinstance(other, Vector4):
            return self.x != other.x and self.y != other.y and self.z != other.z and self.w != other.w
        else:
            raise TypeError("Invalid Argument")
    
    def __neg__(self) -> "Vector4":
        return Vector4(
            -self.x,
            -self.y,
            -self.z,
            -self.w
        )
        
    def __abs__(self) -> "Vector4":
        return Vector4(abs(self.x), abs(self.y), abs(self.z), abs(self.w))
    
    def __getindex__(self, index: int) -> Scalar:
        match (index):
            case 0:
                return self.x
            case 1:
                return self.y
            case 2:
                return self.z
            case 3:
                return self.w
            case _:
                raise IndexError("Index is out of range!")

    def __max__(self, other: "Vector4") -> "Vector4":
        if isinstance(other, Vector4):
            return Vector4(
                max(self.x, other.x),
                max(self.y, other.y),
                max(self.z, other.z),
                max(self.w, other.w)
            )
        else:
            raise TypeError("Invalid Argument")
    
    def __min__(self, other: "Vector4") -> "Vector4":
        if isinstance(other, Vector4):
            return Vector4(
                min(self.x, other.x),
                min(self.y, other.y),
                min(self.z, other.z),
                min(self.w, other.w)
            )
        else:
            raise TypeError("Invalid Argument")

    def clamp(self, x_min: "Vector4", x_max: "Vector4") -> "Vector4":
        if isinstance(x_min, Vector4) and isinstance(x_max, Vector4):
            return min(max(x_min, self), x_max)
        else:
            raise TypeError("Invalid Argument")
    
    def clampf(self, x_min: float, x_max: float) -> "Vector4":
        if isinstance(x_min, float) and isinstance(x_max, float):
            return Vector4(
                clampf(self.x, x_min, x_max),
                clampf(self.y, x_min, x_max),
                clampf(self.z, x_min, x_max),
                clampf(self.w, x_min, x_max)
            )
        else:
            raise TypeError("Invalid Argument")

    def length(self) -> Scalar:
        return (self.x**2 + self.y**2 + self.z**2 + self.w**2)**0.5
    
    def length_squared(self) -> Scalar:
        return (self.x**2 + self.y**2 + self.z**2 + self.w**2)

    def normalized(self) -> "Vector4":
        return self / self.length()

    def direction_to(self, other: "Vector4") -> "Vector4":
        if isinstance(other, Vector4):
            return (other - self).normalized()

    def distance_squared_to(self, other: "Vector4") -> Scalar:
        if isinstance(other, Vector4):
            return (other - self).length_squared()
        else:
            raise TypeError("Invalid Argument")
    
    def distance_to(self, other: "Vector4i") -> Scalar:
        return self.distance_squared_to(other)**0.5
    
    def dot(self, other: "Vector4") -> Scalar:
        return self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    
    def __floor__(self) -> "Vector4":
        return Vector4(
            floor(self.x),
            floor(self.y),
            floor(self.z),
            floor(self.w)
        )
    
    def __round__(self, ndigits=None) -> "Vector4":
        return Vector4(
            round(self.x, ndigits),
            round(self.y, ndigits),
            round(self.z, ndigits),
            round(self.w, ndigits)
        )
    
    def __ceil__(self) -> "Vector4":
        return Vector4(
            ceil(self.x),
            ceil(self.y),
            ceil(self.z),
            ceil(self.w)
        )
    
    def inverse(self) -> "Vector4":
        return Vector4(
            1.0 / self.x,
            1.0 / self.y,
            1.0 / self.z,
            1.0 / self.w
        )
    
    def is_approx_equal(self, other: "Vector4", eps: float = 1e-8) -> bool:
        return (
            is_approx_equal(self.x, other.x, eps) and
            is_approx_equal(self.y, other.y, eps) and
            is_approx_equal(self.z, other.z, eps) and
            is_approx_equal(self.w, other.w, eps)
        )
    
    def is_approx_zero(self, eps: float = 1e-8) -> bool:
        return (
            is_approx_zero(self.x, eps) and
            is_approx_zero(self.y, eps) and
            is_approx_zero(self.z, eps) and
            is_approx_zero(self.w, eps)
        )
    
    def is_finite(self) -> bool:
        return (
            is_finite(self.x) and
            is_finite(self.y) and
            is_finite(self.z) and
            is_finite(self.w)
        )
    
    def is_normalized(self, eps: float = 1e-8) -> bool:
        return abs(self.length() - 1.0) < eps
    
    def lerp(self, other: "Vector4", t: float) -> "Vector4":
        return Vector4(
            lerp(self.x, other.x, t),
            lerp(self.y, other.y, t),
            lerp(self.z, other.z, t),
            lerp(self.w, other.w, t)
        )
    
    def sign(self) -> "Vector4":
        return self / abs(self)
