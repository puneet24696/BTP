from __future__ import division
import sympy
from sympy import Symbol,nsolve
import mpmath
import numpy as np
import matplotlib.pyplot as plt

p = 0  # Weight of point load in Newton
l = 100 #total length of beam
#a_x = 10  #Position of point load
steps = 11
#a_x = [i for i in range(0, l+1, 1 + int(l / steps))]
a = 0
q = 10    # loading per unit length
iy = 1/12  # inertia of beam #CAUTION
e = 200 * 10**9 # youngs modulus
def_c = -1/(e*iy)

#v is shear force, m is moment , d is deflection
'''
*****************
Origin is at the leftmost position. 'x' is taken from the origin
'a' is the position of the point load taken from the origin
's_f' is the shear force at position 'x'.
'b_m' is the bending moment at position 'x'.
*****************
'''

'''
****
LEFT REGION
****
'''
#----data for graphs----
ymin = -0.008
ymax = 0.008
max_3 = []
l2_x = np.linspace(55,60,100)
#l2_x = [57.1,57.15]
for l2 in l2_x:
#    l2 = 50
    l3 = (l - l2)/2
    l1 = l3   #region 3
    ls = l1 + l2 #sum of first and second length
    l1_plot = np.linspace(l1,l1,10)
    ls_plot = np.linspace(ls,ls,10)
    column = np.linspace(0,ymin,10)
    person = np.linspace(0,100,steps)
    zeros = np.linspace(0,0,steps)

    # r1 and r2 are reactions
    a_plot = np.linspace(a,a,steps)
    r2 = q*l/2 - p*(l1-a)/l2
    r1 = q*l + p - r2
    c11 = Symbol('c11')
    c12 = Symbol('c12')
    c21 = Symbol('c21')
    c22 = Symbol('c22')
    c31 = Symbol('c31')
    c32 = Symbol('c32')
    c41 = Symbol('c41')
    c42 = Symbol('c42')

    f1 = c11 - c21 + p * a**2 / 2
    f2 = c11 * a + c12 + p * a**3 / 3 - c21 * a - c22
    f3 = p*(l1**3 /6 - a*l1**2 /2) + q * l1**4 /24 + c21 * l1 + c22
    f4 = c31 - c41 - r2* ls**2 /2
    f5 = c31 * ls + c32 + r2 *(ls**3 /6 - l1 * ls**2 /2 - l2 * ls**2 /2) - c41 *ls -c42
    f6 = c21 - l1**2 / 2 * r1 - c31
    f7 = l1**3 *r1 /3 + c31*l1 + c32 - c21 * l1 - c22
    f8 = q * ls**4 /24 + p* (ls**3 /6 - a* ls**2 /2) - r1* (ls**3 /6 - l1 * ls**2 /2) - r2 * (ls**3 /6 - l1*ls**2 /2 - l2*ls**2 /2) + c41 *ls + c42
    c11,c12,c21,c22,c31,c32,c41,c42 = nsolve((f1,f2,f3,f4,f5,f6,f7,f8),(c11,c12,c21,c22,c31,c32,c41,c42),(10,10,10,10,5,6,7,8))
    x_cordinate = np.linspace(0 , 100 , 101)
    v = []
    m = []
    d = []
    for x in x_cordinate:
	if x < a:               #region 1-a
#	    v.append(-q * x)
#    	    m.append(-q*x**2 /2)
	    d.append(def_c*(q*x**4 /24 + c11 *x + c12))
	elif x >= a and x < l1: #region 1-b
#	    v.append(-p -q*x)
#	    m.append(-p*(x-a)-q*x**2 /2)
	    d.append(def_c*(p*(x**3 /6 - a*x**2 /2)+q*x**4/ 24 +c21*x+c22))
	elif x >= l1 and x < ls: #region 2
#	    v.append(r1 - p - q*x)
#	    m.append(-p*(x-a)-q*x**2 /2 + r1*(x - l1))
	    d.append(def_c*(p*(x**3 /6 - a*x**2 /2) + q * x**4 /24 - r1*(x**3 /6 - l1 * x**2 /2) + c31*x + c32))
	else:                    #region 3
#	    v.append(-q*x -p +r1 +r2)
#	    m.append(-q* x**2 /2 -p*(x-a)+ r1*(x-l1) + r2*(x-l1-l2))
	    d.append(def_c*(q*x**4 /24 + p*(x**3 /6 - a*x**2 /2) - r1*(x**3 /6 - l1*x**2 /2) - r2*(x**3 /6 - l1*x**2 /2 - l2*x**2 /2) + c41*x + c42))
    print l2,"-", d[0],d[51],d[-1]
    max_3.append(max(d[0],d[51],d[-1]))
#    plt.plot(x_cordinate,d,l1_plot,column,'r',ls_plot,column,'r')
#    plt.axis([0,100,ymin,ymax])
#    plt.title("Deflection of beam")
#    plt.show()

print min(max_3), l2_x[max_3.index(min(max_3))]

