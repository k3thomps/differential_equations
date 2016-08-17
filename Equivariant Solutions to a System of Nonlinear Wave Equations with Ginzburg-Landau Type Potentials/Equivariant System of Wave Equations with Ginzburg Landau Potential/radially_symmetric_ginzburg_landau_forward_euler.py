# this is a numerical simulation of the 1d wave equation u_tt - u_xx = 0

import numpy as np
from scipy import exp
# from pylab import *
from matplotlib.pylab import *



# This is the plotting function
def save_graph(x,filename,f,limits,n):
    figure()
    axis(limits)
    
    line, = plot(x,x)

    for i in range(0,n):
        line.set_ydata(f[i])
        draw()
        savefig("graphs/" + '%s%04d' % (filename, i) + ".png")
	if i == 0*n/10:
		print "[>          ]"
	elif i == 1*n/10:
		print "[->         ]"
	elif i == 2*n/10:
		print "[-->        ]"
	elif i == 3*n/10:
		print "[--->       ]"
	elif i == 4*n/10:
		print "[---->      ]"
	elif i == 5*n/10:
		print "[----->     ]"
	elif i == 6*n/10:
		print "[------>    ]"
	elif i == 7*n/10:
		print "[------->   ]"
	elif i == 8*n/10:
		print "[-------->  ]"
	elif i == 9*n/10:
		print "[---------> ]"
	elif i == 10*n/10:
		print "[---------->]"


# domain = [a,b]
# let dr = |[a,b]|/m = (b-a)/m, where m is the number of divisions of the spatial axis
a = 100.
b = 101.
m = 500
dr = (b-a)/m
mid = (b+a)/2.

# epsilon parameter
epsilon = 1./200.

# we will discritize the interval [a,b] into m equally sized intervals and this will be the spatial lattice we will do the numerics on
pts = np.arange(a,b + dr/2.,dr)

# time steps
dt = 0.00005

# number of time steps
maxiter = 2500

# we set the initial data here. For convenience we use: u_0 = hyperbolic tangent and partial_t u_0 = 0
u0 = (exp((pts - mid)/epsilon) - exp(-(pts - mid)/epsilon))/(exp((pts - mid)/epsilon) + exp(-(pts - mid)/epsilon))
u1 = u0*1.
u2 = u0*1.

u = []
u.append(u0*1.)
u.append(u1*1.)

# we set a counting variable i and a common quantity that shows up (dt/dr)^2
i = 1
cfl = (dt/dr)**2

# forward method (i.e. for the next time step, a position only depends on the values of the 2 previous time steps and not on any of its neighbours

while (i <= maxiter):

	for j in range(1,m):
		u2[j] = 2.*u1[j] - u0[j] + (cfl/pts[j])*( (pts[j] + dr/2.)*u1[j+1] - 2.*pts[j]*u1[j] + (pts[j] - dr/2.)*u1[j-1] ) - ((dt/epsilon)**2.)*(u1[j]**2. - 1.)*u1[j]
	
	# Neumann boundary conditions
	u2[0] = 2.*u1[0] - u0[0] + (cfl/pts[0])*((pts[0] + dr/2.)*u1[0+1] - 2.*pts[0]*u1[0] + (pts[0] - dr/2.)*u1[0]) - ((dt/epsilon)**2.)*(u1[0]**2. - 1.)*u1[0]
	u2[m] = 2.*u1[m] - u0[m] + (cfl/pts[m])*((pts[m] + dr/2.)*u1[m] - 2.*pts[m]*u1[m] + (pts[m] - dr/2.)*u1[m-1]) - ((dt/epsilon)**2.)*(u1[m]**2. - 1.)*u1[m]

	# update
	u.append(u2*1.)
	u0 = u1*1.
	u1 = u2*1.
	i = i + 1

print "-----> Done simulation. Saving graphs. <-----"
save_graph(pts,"test",u,limits=[a,b,-2,2],n=maxiter)

