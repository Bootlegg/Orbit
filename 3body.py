import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.animation as animation

#units
AU = 1.5*10**11
M = 2*10**30
me = 6*10**24
mv = 4.87*10*24

re = 6371*10**3
rv = 6052*10**3
rs = 695700*10**3
G = 6.67*10**(-11)
dt = 300000 #seems about right


#iterations
nt = 100

x = np.linspace(-2*AU,2*AU)
y = np.linspace(-2*AU,2*AU)

#earth
xe = np.zeros((1+nt,1))
ye = np.zeros((1+nt,1))
ue = np.zeros((1+nt,1))
ve = np.zeros((1+nt,1))

#sun
xs = np.zeros((1+nt,1))
ys = np.zeros((1+nt,1))
us = np.zeros((1+nt,1))
vs = np.zeros((1+nt,1))

#venus
xv = np.zeros((1+nt,1))
yv = np.zeros((1+nt,1))
uv = np.zeros((1+nt,1))
vv = np.zeros((1+nt,1))

#starting conditions
xe[0,0] = AU
ve[0,0] = 29800

xv[0,0] = -108200000*10**3
vv[0,0] = -35000

for n in range(nt):
	#the equations
	#dxe/dt = ue
	#due/dt = Fx/m



	#ue[1] = ue[0]+dt*Fx/m
	#ve[1] = ve[0]+dt*Fy/m
	#earth
	ue[n+1,0] = ue[n,0]+dt*(-G*M/((xe[n,0]-xs[n,0])**2+(ye[n,0]-ys[n,0])**2)**(1.5))*(xe[n,0]-xs[n,0])\
								+dt*(-G*mv/((xe[n,0]-xv[n,0])**2+(ye[n,0]-yv[n,0])**2)**(1.5))*(xe[n,0]-xv[n,0])
	
	ve[n+1,0] = ve[n,0]+dt*(-G*M/((xe[n,0]-xs[n,0])**2+(ye[n,0]-ys[n,0])**2)**(1.5))*(ye[n,0]-ys[n,0])\
							+dt*(-G*mv/((xe[n,0]-xv[n,0])**2+(ye[n,0]-yv[n,0])**2)**(1.5))*(ye[n,0]-yv[n,0])
	#sun
	us[n+1,0] = us[n,0]+dt*(-G*me/((xe[n,0]-xs[n,0])**2+(ye[n,0]-ys[n,0])**2)**(1.5))*(xs[n,0]-xe[n,0])\
							+dt*(-G*mv/((xv[n,0]-xs[n,0])**2+(yv[n,0]-ys[n,0])**2)**(1.5))*(xs[n,0]-xv[n,0])
								
	
	vs[n+1,0] = vs[n,0]+dt*(-G*me/((xe[n,0]-xs[n,0])**2+(ye[n,0]-ys[n,0])**2)**(1.5))*(ys[n,0]-ye[n,0])\
							+dt*(-G*mv/((xv[n,0]-xs[n,0])**2+(yv[n,0]-ys[n,0])**2)**(1.5))*(ys[n,0]-yv[n,0])
	
	#venus
	uv[n+1,0] = uv[n,0]+dt*(-G*M/((xv[n,0]-xs[n,0])**2+(yv[n,0]-ys[n,0])**2)**(1.5))*(xv[n,0]-xs[n,0])\
						+dt*(-G*me/((xv[n,0]-xe[n,0])**2+(yv[n,0]-ye[n,0])**2)**(1.5))*(xv[n,0]-xe[n,0])
						
						
	vv[n+1,0] = vv[n,0]+dt*(-G*M/((xv[n,0]-xs[n,0])**2+(yv[n,0]-ys[n,0])**2)**(1.5))*(yv[n,0]-ys[n,0])\
						+dt*(-G*me/((xv[n,0]-xe[n,0])**2+(yv[n,0]-ye[n,0])**2)**(1.5))*(yv[n,0]-ye[n,0])
	
	#discretized
	#forward, explicit
	#earth
	xe[n+1,0] = xe[n,0]+dt*ue[n+1,0]
	ye[n+1,0] = ye[n,0]+dt*ve[n+1,0]
	#sun
	xs[n+1,0] = xs[n,0]+dt*us[n+1,0]
	ys[n+1,0] = ys[n,0]+dt*vs[n+1,0]
	#venus
	xv[n+1,0] = xv[n,0]+dt*uv[n+1,0]
	yv[n+1,0] = yv[n,0]+dt*vv[n+1,0]
	
	
plt.hold(True)
fig = plt.figure()
ax = fig.gca()
ax.set_xlim(-2*AU,2*AU)
ax.set_ylim(-2*AU,2*AU)
ax.set_title('Kappa')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')
ax.scatter(xe[0,0], ye[0,0],s=30)
ax.scatter(xv[0,0], yv[0,0],s=30)
ax.scatter(xs[0,0], ys[0,0],s=100, c='yellow')



def animate(i): #i increment with 1 each step
	ax.clear()
	ax.set_xlim(-2*AU,2*AU)
	ax.set_ylim(-2*AU,2*AU)
	ax.set_title('Kappa')
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	ax.scatter(xe[i,0], ye[i,0],s=30)
	ax.scatter(xv[i,0], yv[i,0],s=30)
	ax.scatter(xs[i,0], ys[i,0],s=100,c='yellow')
	if i == 50:
		plt.savefig('Orbit.png', bbox_inches='tight')
	return None


anim = animation.FuncAnimation(fig, animate, frames=nt, interval=100)

plt.show()


