import cmath
import math

ALPHA = cmath.exp(2/3*cmath.pi*1j)

class ThreePhase:
	def __init__(self, a: complex = 0, b: complex = 0, c: complex = 0):
		self._a = a
		self._b = b
		self._c = c
		self._theta = 0
		self._update_from_abc()

	def __add__(self, other):
		if not isinstance(other, type(self)):
			raise TypeError(f"Unsupported operand type(s) for +: '{type(self).__name__}' and '{type(other).__name__}'")
		return type(self)(a=self._a + other.a, b=self._b + other.b, c=self._c + other.c)

	def __sub__(self, other):
		if not isinstance(other, type(self)):
			raise TypeError(f"Unsupported operand type(s) for -: '{type(self).__name__}' and '{type(other).__name__}'")
		return type(self)(a=self._a - other.a, b=self._b - other.b, c=self._c - other.c)
	
	def __mul__(self, other):
		if not isinstance(other, type(self)):
			raise TypeError(f"Unsupported operand type(s) for *: '{type(self).__name__}' and '{type(other).__name__}'")
		return type(self)(a=self._a * other.a, b=self._b * other.b, c=self._c * other.c)	

	def __div__(self, other):
		if not isinstance(other, type(self)):
			raise TypeError(f"Unsupported operand type(s) for /: '{type(self).__name__}' and '{type(other).__name__}'")
		return type(self)(a=self._a * other.a, b=self._b * other.b, c=self._c * other.c)

	def __rlshift__(self, angle):
		if not (isinstance(angle, float) or isinstance(angle, int)):
			raise TypeError(f"Unsupported operand type(s) for <<: '{type(self).__name__}' and '{type(angle).__name__}'")
		_shift = cmath.exp(angle*1j)
		return type(self)(a=self._a * _shift, b=self._b * _shift, c=self._c * _shift)
	
	def __rrshift__(self, angle):
		if not (isinstance(angle, float) or isinstance(angle, int)):
			raise TypeError(f"Unsupported operand type(s) for >>: '{type(self).__name__}' and '{type(angle).__name__}'")
		_shift = cmath.exp(-angle*1j)
		return type(self)(a=self._a * _shift, b=self._b * _shift, c=self._c * _shift)
			
	def _update_fortescue_components(self):
		self._z = (self._a + self._b + self._c)/3
		self._p = (self._a + self._b*ALPHA + self._c*ALPHA**2)/3
		self._n = (self._a + self._b*ALPHA**2 + self._c*ALPHA)/3

	def _update_clarke_components(self):
		_k = 2/3
		self._x = _k * (self._a - self._b/2 - self._c/2)
		self._y = _k * (self._b * cmath.sqrt(3)/2 - self._c * cmath.sqrt(3)/2)
		#self._z = _k * (self._a/2 + self._b/2 + self._c/2)	

	def _update_park_components(self):
		self._d = cmath.cos(self._theta)*self._x + cmath.sin(self._theta)*self._y
		self._q = -cmath.sin(self._theta)*self._x + cmath.cos(self._theta)*self._y
		
	def _update_from_abc(self):
		self._update_fortescue_components()
		self._update_clarke_components()
		self._update_park_components()

	def _update_from_zpn(self):
		self._a = self._z + self._p + self._n
		self._b = self._z + self._p*ALPHA**2 + self._n*ALPHA
		self._c = self._z + self._p*ALPHA + self._n*ALPHA**2
		self._update_from_abc()

	def _update_from_xy(self):
		self._a = self._x + self._z
		self._b = -self._x/2 + self._y * cmath.sqrt(3)/2 + self._z
		self._c = -self._x - self._y * cmath.sqrt(3)/2 + self._z
		self._update_from_abc()

	def _update_from_dq(self):
		self._x = cmath.cos(self._theta)*self._d - cmath.sin(self._theta)*self._q
		self._y = cmath.sin(self._theta)*self._d + cmath.cos(self._theta)*self._q
		#self._z = -self._x - self._y * cmath.sqrt(3)/2 + self._z
		self._update_from_xy()

	@property
	def a(self):
		return self._a
	
	@a.setter
	def a(self, value: complex):
		self._a = value
		self._update_from_abc()

	@property
	def b(self):
		return self._b
	
	@b.setter
	def b(self, value: complex):
		self._b = value
		self._update_from_abc()

	@property
	def c(self):
		return self._c
	
	@c.setter
	def c(self, value: complex):
		self._c = value
		self._update_from_abc()

	@property
	def p(self):
		return self._p
	
	@p.setter
	def p(self, value: complex):
		self._p = value
		self._update_from_zpn()

	@property
	def n(self):
		return self._n
	
	@n.setter
	def n(self, value: complex):
		self._n = value
		self._update_from_zpn()

	@property
	def z(self):
		return self._z
	
	@z.setter
	def z(self, value: complex):
		self._z = value
		self._update_from_zpn()

	@property
	def x(self):
		return self._x
	
	@x.setter
	def x(self, value: complex):
		self._x = value
		self._update_from_xy()

	@property
	def y(self):
		return self._y
	
	@y.setter
	def y(self, value: complex):
		self._y = value
		self._update_from_xy()

	@property
	def d(self):
		return self._d
	
	@d.setter
	def d(self, value: complex):
		self._d = value
		self._update_from_dq()

	@property
	def q(self):
		return self._q
	
	@q.setter
	def q(self, value: complex):
		self._q = value
		self._update_from_dq()

	@property
	def theta(self):
		return self._theta
	
	@theta.setter
	def theta(self, value: complex):
		self._theta = value
		self._update_park_components()

	def to_dict(self, polar: bool = False, letters: bool = True) -> dict:
		if letters:
			_keys = ['A', 'B', 'C', 'P', 'N', 'Z', 'X', 'Y', 'D', 'Q']
		else:
			_keys = ['A', 'B', 'C', '+', '-', '0', 'X', 'Y', 'D', 'Q']
			
		_values = [self._a, self._b, self._c, self._p, self._n, self._z, self._x, self._y, self._d, self._q]

		if polar:
			for i, val in enumerate(_values):
				_mag, _ang = cmath.polar(_values[i])
				_ang = math.degrees(_ang)

				_values[i] = (_mag, _ang)

		return dict(zip(_keys, _values))
	
	def to_abc_dict(self, polar: bool = False):
		_keys = ['A', 'B', 'C']
		_dict = self.to_dict(polar=polar)
	
		return {k: _dict[k] for k in _keys if k in _dict}

	def to_dqz_dict(self, polar: bool = False, letters: bool = True) -> dict:
		if letters:
			_keys = ['D', 'Q', 'Z']
		else:
			_keys = ['D', 'Q', '0']

		_dict = self.to_dict(polar=polar, letters=letters)
	
		return {k: _dict[k] for k in _keys if k in _dict}
	
	def to_xyz_dict(self, polar: bool = False, letters: bool = True) -> dict:
		if letters:
			_keys = ['X', 'Y', 'Z']
		else:
			_keys = ['X', 'Y', '0']
			
		_dict = self.to_dict(polar=polar, letters=letters)
	
		return {k: _dict[k] for k in _keys if k in _dict}
	
	def to_pnz_dict(self, polar: bool = False, letters: bool = True) -> dict:
		if letters:
			_keys = ['P', 'N', 'Z']
		else:
			_keys = ['+', '-', '0']
			
		_dict = self.to_dict(polar=polar, letters=letters)
	
		return {k: _dict[k] for k in _keys if k in _dict}

	def conjugate(self):
		return type(self)(a=self._a.conjugate(), b=self._b.conjugate(), c=self._c.conjugate())


class WyeDelta:
	def __init__(self, z1: complex = 0, z2: complex = 0, z3: complex = 0):
		self._z1 = z1
		self._z2 = z2
		self._z3 = z3
		self._update_from_wye()

	def _update_from_wye(self):
		_zp = self._z1 * self._z2 + self._z2 * self._z3 + self._z1 * self._z3
		self._za = _zp/self._z1
		self._zb = _zp/self._z2
		self._zc = _zp/self._z3

	def _update_from_delta(self):
		_zsum = self._za + self._zb + self._zb
		self._z1 = (self._zb * self._zc)/_zsum
		self._z2 = (self._za * self._zc)/_zsum
		self._z3 = (self._za * self._zb)/_zsum

	@property
	def z1(self):
		return self._z1
	
	@z1.setter
	def z1(self, value: complex):
		self._z1 = value
		self._update_from_wye()

	@property
	def z2(self):
		return self._z2
	
	@z2.setter
	def z2(self, value: complex):
		self._z2 = value
		self._update_from_wye()

	@property
	def z3(self):
		return self._z3
	
	@z3.setter
	def z3(self, value: complex):
		self._z3 = value
		self._update_from_wye()

	@property
	def za(self):
		return self._za
	
	@za.setter
	def za(self, value: complex):
		self._za = value
		self._update_from_delta()

	@property
	def zb(self):
		return self._zb
	
	@zb.setter
	def zb(self, value: complex):
		self._zb = value
		self._update_from_delta()

	@property
	def zc(self):
		return self._zc
	
	@zc.setter
	def zc(self, value: complex):
		self._zc = value
		self._update_from_delta()
