from __future__ import division
from sympy.solvers import nsolve
from sympy import Symbol
import sympy

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
a = Symbol('a')

eq1 = a - 0.0071 +0.0071*sympy.cosh(-14057.73/y) + 100*sympy.sinh(-140.5773/y)/y
eq2 = -50/(z**2) + a
eq3 = -2354400/(x**2) +a
eq4 = x+y+z-1138000
print nsolve((eq1,eq2,eq3,eq4),(x,y,z,a),(5000.0,50000.0,50.0,100.0))

