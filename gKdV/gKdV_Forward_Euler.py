# this is a numerical simulation of the gKdV using Forward Euler u_t + u_xxx + (u^p)_x = 0

import numpy as np
from scipy import exp
from matplotlib.pylab import *



# This is the plotting function
def save_graph(x,filename,f,limits,n):
    figure()
    axis(limits)
    
    line, = plot(x,x)

    for i in range(0,n):
	print i
        line.set_ydata(f[i])
        draw()
        savefig("graphs/" + '%s%04d' % (filename, i) + ".png")



# domain = [0,1]
# let h = |[0,1]|/m = 1/m, where m is the number of divisions of the spatial axis
a = -10.
b = 10.
m = 1000
dr = (b-a)/m

# we will discritize the interval [a,b] into m equally sized intervals and this will be the spatial lattice we will do the numerics on
pts = np.arange(a,b + dr/2.,dr)

# time steps
dt = 0.0001

# maximum number of iterations
i = 1
maxiter = 100

# parameter p appearing in gKdV
p = 3

# we set the initial data here. We just use a soliton for convenience.
u0 = ( (p+1)/( exp((p-1)/2*pts) + exp(-(p-1)/2*pts) ) )**(1/(p-1))
n = len(pts)
u1 = zeros(n)

u = []
u.append(u0*1.)

# forward method (i.e. for the next time step, a position only depends on the values of the 2 previous time steps and not on any of its neighbours

while (i <= maxiter):
	for j in range(1,n-1):
		u1[j] = u0[j] - dt/(2*dr**3)( u0[j+2] - 2*u0[j+1] + 2*u0[j-1] - u0[j - 2] ) - (dt/dr)*( u0[j+1]**p - u0[j-1]**p )

	u0[1] = u0[1] - dt/(2*dr**3)( u0[1+2] - 2*u0[1+1] + 2*u0[1-1] - u0[1 - 1] ) - (dt/dr)*( u0[1+1]**p - u0[1-1]**p )
	u0[0] = u0[0] - dt/(2*dr**3)( u0[0+2] - 2*u0[0+1] + 2*u0[0] - u0[0] ) - (dt/dr)*( u0[0+1]**p - u0[0]**p )

	u[n-1] = u0[n-1] - dt/(2*dr**3)( u0[n-1+1] - 2*u0[n-1+1] + 2*u0[n-1-1] - u0[n-1 - 2] ) - (dt/dr)*( u0[n-1+1]**p - u0[n-1-1]**p )
	u[n] = u0[n] - dt/(2*dr**3)( u0[n] - 2*u0[n] + 2*u0[n-1] - u0[n - 2] ) - (dt/dr)*( u0[n]**p - u0[n-1]**p )


	# update
	u.append(u2*1.)
	u0 = u1*1.
	i = i + 1

save_graph(pts,"test",u,limits=[-1,1,a,b],n=maxiter,)

