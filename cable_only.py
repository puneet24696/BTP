from __future__ import division
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import math

def get_x(n,l):
    x = []
    for i in range(n):
        x.append(i*l/(n-1))
    return x

#catenary start
def initial_slope_cat(w,T_0):
    return math.sinh(-w*l/(2*T_0))

def max_tension_cat(w,T_0):
    dy_dx = initial_slope_cat(w,T_0)
    return T_0* (1 + dy_dx**2)**0.5

def get_y_cat(x,w,T_0):
    y= []
    for i in range(len(x)):
        y.append((T_0 * (math.cosh(w* (x[i]- l/2)/T_0))/w - T_0 * (math.cosh(w* (-l/2)/T_0))/w))
    return y
#catenary end

### parabola start
def max_tension_para(w,d_s,l):
    return -(w* l* (16* d_s**2 + l**2 )**0.5)/(8*d_s)  #max tension at corner

def max_tension_para2(T_par0,d_s):
    dy_dx = -d_s*l/(2*T_par0)
    return T_par0* (1+ dy_dx**2)**0.5

def get_y_from_t(x,d_s,T_par0):
    y = []
    for i in range(len(x)):
        y.append(d_s * x[i]**2/(2*T_par0)- d_s*l*x[i]/(2*T_par0))
    return y

def get_y_para(x,d_s,l):
    y = []
    for i in range(len(x)):
        y.append(-4*d_s*x[i]**2/l**2 + 4*d_s*x[i]/l)
    return y
### parabola end

####### point start
def tension_point(x_load,d_p,T_P0):
    if d_p!=0:
        return T_P0* (d_p**2 + x_load**2)**0.5/x_load

def get_y_point(x,T_P0):
    y = x*(l-x)*W/(l*T_P0)
    return -1*y

def get_y_pointshape(x,x_load,d_p):
    y = []
    for i in range(len(x)):
        if d_p < 0:
            if x[i] < x_load:
                y.append(x[i] *d_p/x_load)
            elif x[i] == x_load:
                y.append(d_p)
            else:
                y.append((x[i]-l)*d_p/(x_load-l))
        else:
            y.append(0)
    return y

def elast_get_y(x,e):
    d = Symbol('d')
    d_elas = solve(((x**2 +d **2)**0.5 + ((l-x)**2 + d**2)**0.5 -l )/l - e ,d)
    print d_elas
    return d_elas

def get_e_with_t(T):
    return  T/(A*E)
###### point end



def plot_cable(x,y,fig_no,name,leg):
    plt.plot(x,y,label = leg)
    plt.plot([x[0],x[0]],[-10,10],"k")
    plt.plot([x[-1],x[-1]],[-10,10],"k")
    plt.plot([x[0],x[-1]],[0,0],"k")
    plt.xlabel("length of cable")
    plt.ylabel("deflection")
    plt.title(name)
    plt.legend()
    #plt.gca().set_aspect("equal","box")
    #plt.axis("equal")
    #plt.xlim([min(x),max(x)])
    #plt.ylim([-1*max(x),-1*min(x)])
    plt.savefig(fig_no)
    plt.clf()


###MAIN ## all dimensions in MKS
w = 14.33 # weight/length kg/m of cable 60mm dia, IWRC(independent wire rope core)
g = 9.81 #gravity
w = w*g
w_s = 24*g #suspended weight/unit length
mbl = (2268 * 10**3)/2 #min breaking limit(N) will give safety factor of 2
print "mbl ="
print mbl
W = 10000 #point load weight
A = 3.14* 0.03**2 #area
E = 158* 10**9 #youngs modulus
d_s = 2  #max deflection(centre) due to suspended weight
d_s = (-1)*d_s
#d_w = 1        # max deflection(centre) due to self weight
#d_w = -1*d_w
l = 200 # length of cable
n = 200 # number of points for plot.

x = []
x = get_x(n,l)
'''
#parabola
e = get_e_with_t(mbl)
print e
'''
strplot_t = [350000]
for T_par0 in strplot_t:
    t_para =  max_tension_para2(T_par0,d_s)
    y_para = get_y_from_t(x,w,T_par0)
#    plot_cable(x,y_para,"para"+str(T_par0),"Cable as Parabola","Max Deflection = "+str(abs(min(y_para))))
#    print min(y_para)

##### catenary
y_cat = []
#T_0 = (mbl - abs(t_pt + t_para))

#strplot_t = [1000000,800000,600000,400000,200000,100000]
strplot_t = [350000]
for T_0 in strplot_t:
    t_cat = max_tension_cat(w,T_0)
    print "t_cat"
    print t_cat
    y_cat = get_y_cat(x,w,T_0)
#    plot_cable(x,y_cat,"cat","Cable as Catenary","Max Deflection = "+str(abs(min(y_cat))))

###point

#strplot_t = [1000000,500000,100000,50000]
strplot_t = [1000000-T_0-T_par0]
y_pt = []
n_load = 1 # number of points where the load is to be applied.
#x_load = get_x(n_load,l)
#print x_load
x_load = [100]
for T_P0 in strplot_t:
    for i in range(n_load):
        d_p = get_y_point(x_load[i],T_P0) #deflection due to load
        y_pt = get_y_pointshape(x,x_load[i],d_p)
#        plot_cable(x,y_pt,"point"+str(8877),"Cable with point load","Max Deflection = "+str(abs(d_p)))
        t_pt = tension_point(x_load[i],d_p,T_P0)
        print "t_pt"
        print t_pt
	y_rope = []
	for i in range(len(x)):
    	    y_rope.append(y_cat[i]+y_pt[i])
	d_rope = abs(min(y_cat))+abs(min(y_pt)) + abs(min(y_para))
	print d_rope
	t_rope = t_pt + t_cat + t_para
	print t_rope

	plot_cable(x,y_rope,"suspension"+str(0),"Main cable of Suspension  bridge","Max deflection = "+str(d_rope))
#	print t_pt
	d_rope = abs(min(y_cat))+abs(min(y_pt)) + abs(min(y_para))
	print d_rope
	t_rope = t_pt + t_cat + t_para
	print t_rope


##rope
'''
y_rope = []
for i in range(len(x)):
    y_rope.append(y_cat[i]+y_pt[i])

plot_cable(x,y_rope,"rope"+str(x_load[0]),"Ropeway","Max deflection = "+str(abs(min(y_cat))+abs(min(y_pt))))
print t_pt
d_rope = abs(min(y_cat))+abs(min(y_pt))
print d_rope
t_rope = t_pt + t_cat
print t_rope
'''
'''
def cosT(tana):
    return 1/(1+tana**2)**0.5
##Deflection elastic
d_elas = [-4.6949,-6.1326,-6.8519,-7.075,-6.8519,-6.1326,-4.6949]
xload = [25,50,75,100,125,150,175]
for  i in range(len(d_elas)):
    y_elas = get_y_pointshape(x,xload[i],d_elas[i])
    tana = d_elas[i]/xload[i]
    cosa = cosT(tana)
    tanb = d_elas[i]/(l-xload[i])
    cosb = cosT(tanb)
    p_load = mbl*(cosa+cosb)
    print p_load
    #print y_elas
    plot_cable(x,y_elas,"elastic"+str(i),"Elastic Cable","Deflection = "+str(-d_elas[i]))

'''


