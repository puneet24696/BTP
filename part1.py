from __future__ import division
import sympy
from sympy import Symbol,nsolve
import mpmath
import numpy as np
import matplotlib.pyplot as plt

p = 6000  # Weight of point load in Newton
l = 100 #total length of beam
l1 = 30   #region 1
l2 = 40   #region 2
l3 = 30   #region 3
ls = l1 + l2 #sum of first and second length
#a_x = 10  #Position of point load
steps = 5
a_x = np.linspace(0,l,steps)
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
vmin = -10000
vmax = 10000
mmin = -200000
mmax = 200000
ymin = -0.07
ymax = 0.07
l1_plot = np.linspace(l1,l1,10)
ls_plot = np.linspace(ls,ls,10)
column = np.linspace(0,ymin,10)
person = np.linspace(0,100,steps)
zeros = np.linspace(0,0,steps)
for a in a_x:
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
    if a < l1:
        f1 = c11 - c21 + p * a**2 / 2
        f2 = c11 * a + c12 + p * a**3 / 3 - c21 * a - c22
        f3 = p*(l1**3 /6 - a*l1**2 /2) + q * l1**4 /24 + c21 * l1 + c22
        f4 = c31 - c41 - r2* ls**2 /2
        f5 = c31 * ls + c32 + r2 *(ls**3 /6 - l1 * ls**2 /2 - l2 * ls**2 /2) - c41 *ls -c42
        f6 = c21 - l1**2 / 2 * r1 - c31
        f7 = l1**3 *r1 /3 + c31*l1 + c32 - c21 * l1 - c22
        f8 = q * ls**4 /24 + p* (ls**3 /6 - a* ls**2 /2) - r1* (ls**3 /6 - l1 * ls**2 /2) - r2 * (ls**3 /6 - l1*ls**2 /2 - l2*ls**2 /2) + c41 *ls + c42
        c11,c12,c21,c22,c31,c32,c41,c42 = nsolve((f1,f2,f3,f4,f5,f6,f7,f8),(c11,c12,c21,c22,c31,c32,c41,c42),(10,10,10,10,5,6,7,8))

        x_cordinate = np.linspace(0 , 100 , 100)
        v = []
        m = []
        d = []
        for x in x_cordinate:
	    if x < a:               #region 1-a
	        v.append(-q * x)
    	        m.append(-q*x**2 /2)
	        d.append(def_c*(q*x**4 /24 + c11 *x + c12))
	    elif x >= a and x < l1: #region 1-b
	        v.append(-p -q*x)
	        m.append(-p*(x-a)-q*x**2 /2)
	        d.append(def_c*(p*(x**3 /6 - a*x**2 /2)+q*x**4/ 24 +c21*x+c22))
	    elif x >= l1 and x < ls: #region 2
	        v.append(r1 - p - q*x)
	        m.append(-p*(x-a)-q*x**2 /2 + r1*(x - l1))
	        d.append(def_c*(p*(x**3 /6 - a*x**2 /2) + q * x**4 /24 - r1*(x**3 /6 - l1 * x**2 /2) + c31*x + c32))
	    else:                    #region 3
	        v.append(-q*x -p +r1 +r2)
	        m.append(-q* x**2 /2 -p*(x-a)+ r1*(x-l1) + r2*(x-l1-l2))
	        d.append(def_c*(q*x**4 /24 + p*(x**3 /6 - a*x**2 /2) - r1*(x**3 /6 - l1*x**2 /2) - r2*(x**3 /6 - l1*x**2 /2 - l2*x**2 /2) + c41*x + c42))
        plt.figure(1)
        plt.subplot(411)
        plt.plot(a_plot,person)
        plt.axis([0,100,0,100])
        plt.subplot(412)
        plt.plot(x_cordinate,v,person,zeros)
        plt.axis([1,100,vmin,vmax])
        plt.subplot(413)
        plt.plot(x_cordinate,m,person,zeros)
        plt.axis([0,100,mmin,mmax])
        plt.subplot(414)
        plt.plot(x_cordinate,d,l1_plot,column,'r',ls_plot,column,'r')
        plt.axis([0,100,ymin,ymax])
        plt.show()

    elif a >= l1 and a < ls:
        f1 = c21 - c31 + p * a**2 / 2
        f2 = c21 * a + c22 + p * a**3 / 3 - c31 * a - c32
        f3 = q* l1**4 /24 + c11*l1 + c12
        f4 = c31 - c41 - r2* ls**2 /2
        f5 = c31 * ls + c32 + r2 *(ls**3 /6 - l1 * ls**2 /2 - l2 * ls**2 /2) - c41 *ls -c42
        f6 = c11 - l1**2 / 2 * r1 - c21
        f7 = l1**3 *r1 /3 + c21*l1 + c22 - c11 * l1 - c12
        f8 = q * ls**4 /24 + p* (ls**3 /6 - a* ls**2 /2) - r1* (ls**3 /6 - l1 * ls**2 /2) - r2 * (ls**3 /6 - l1*ls**2 /2 - l2*ls**2 /2) + c41 *ls + c42
        c11,c12,c21,c22,c31,c32,c41,c42 = nsolve((f1,f2,f3,f4,f5,f6,f7,f8),(c11,c12,c21,c22,c31,c32,c41,c42),(10,10,10,10,5,6,7,8))

        x_cordinate = np.linspace(0 , 100 , 100)
        v = []
        m = []
        d = []
        for x in x_cordinate:
	    if x < l1:               #region 1
	        v.append(-q * x)
	        m.append(-q*x**2 /2)
	        d.append(def_c*(q*x**4 /24 + c11 *x + c12))
	    elif x >= l1 and x < a: #region 2-a
	        v.append(r1 -q*x)
	        m.append(r1*(x-l1) - q*x**2 /2)
	        d.append(def_c*(-r1*(x**3 /6 - l1*x**2 /2)+q*x**4/ 24 +c21*x+c22))
	    elif x >= a and x < ls: #region 2-b
	        v.append(r1 - p - q*x)
	        m.append(-p*(x-a)-q*x**2 /2 + r1*(x - l1))
	        d.append(def_c*(p*(x**3 /6 - a*x**2 /2) + q * x**4 /24 - r1*(x**3 /6 - l1 * x**2 /2) + c31*x + c32))
	    else:                    #region 3
	        v.append(-q*x -p +r1 +r2)
	        m.append(-q* x**2 /2 -p*(x-a)+ r1*(x-l1) + r2*(x-l1-l2))
	        d.append(def_c*(q*x**4 /24 + p*(x**3 /6 - a*x**2 /2) - r1*(x**3 /6 - l1*x**2 /2) - r2*(x**3 /6 - l1*x**2 /2 - l2*x**2 /2) + c41*x + c42))

        plt.figure(1)
        plt.subplot(411)
        plt.plot(a_plot,person)
        plt.axis([0,100,0,100])
        plt.subplot(412)
        plt.plot(x_cordinate,v,person,zeros)
        plt.axis([1,100,vmin,vmax])
        plt.subplot(413)
        plt.plot(x_cordinate,m,person,zeros)
        plt.axis([0,100,mmin,mmax])
        plt.subplot(414)
        plt.plot(x_cordinate,d,l1_plot,column,'r',ls_plot,column,'r')
        plt.axis([0,100,ymin,ymax])
        plt.show()



#################if shatement 3 for region3 a>ls

    else:                          #  a > ls:
        f1 = c21 - c31 + r2 * ls**2 / 2
        f2 = c31 * a + c32 + p * a**3 / 3 - c41 * a - c42
        f3 = q* l1**4 /24 + c11*l1 + c12
        f4 = c31 - c41 + p* a**2 /2
        f5 = c21 * ls + c22 - r2 *ls**3 /3 - c31 *ls -c32
        f6 = c11 - l1**2 / 2 * r1 - c21
        f7 = l1**3 *r1 /3 + c21*l1 + c22 - c11 * l1 - c12
        f8 = q * ls**4 /24  - r1* (ls**3 /6 - l1 * ls**2 /2) + c21 *ls + c22
        c11,c12,c21,c22,c31,c32,c41,c42 = nsolve((f1,f2,f3,f4,f5,f6,f7,f8),(c11,c12,c21,c22,c31,c32,c41,c42),(10,10,10,10,5,6,7,8))

        x_cordinate = np.linspace(0 , 100 , 100)
        v = []
        m = []
        d = []
        for x in x_cordinate:
	    if x < l1:               #region 1
	        v.append(-q * x)
	        m.append(-q*x**2 /2)
	        d.append(def_c*(q*x**4 /24 + c11 *x + c12))
	    elif x >= l1 and x < ls: #region 2
	        v.append(r1 -q*x)
	        m.append(r1*(x-l1) - q*x**2 /2)
	        d.append(def_c*(-r1*(x**3 /6 - l1*x**2 /2)+q*x**4/ 24 +c21*x+c22))
	    elif x >= ls and x < a: #region 3-a
	        v.append(r1 + r2 - q*x)
	        m.append(r2*(x-ls) -q*x**2 /2 + r1*(x - l1))
	        d.append(def_c*(-r2*(x**3 /6 - ls*x**2 /2) + q * x**4 /24 - r1*(x**3 /6 - l1 * x**2 /2) + c31*x + c32))
	    else:                    #region 3-b
	        v.append(-q*x -p +r1 +r2)
	        m.append(-q* x**2 /2 -p*(x-a)+ r1*(x-l1) + r2*(x-l1-l2))
	        d.append(def_c*(q*x**4 /24 + p*(x**3 /6 - a*x**2 /2) - r1*(x**3 /6 - l1*x**2 /2) - r2*(x**3 /6 - l1*x**2 /2 - l2*x**2 /2) + c41*x + c42))

        plt.figure(1)
        plt.subplot(411)
        plt.plot(a_plot,person)
        plt.axis([0,100,0,100])
        plt.subplot(412)
        plt.plot(x_cordinate,v,person,zeros)
        plt.axis([1,100,vmin,vmax])
        plt.subplot(413)
        plt.plot(x_cordinate,m,person,zeros)
        plt.axis([0,100,mmin,mmax])
        plt.subplot(414)
        plt.plot(x_cordinate,d,l1_plot,column,'r',ls_plot,column,'r')
        plt.axis([0,100,ymin,ymax])
        plt.show()



