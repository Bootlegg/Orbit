import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.animation as animation



#units
AU = 1.5*10**11
M = 2*10**30
me = 6*10**24
mv = 4.87*10*24
mm = 3.3*10**23
mma = 6.417*10**23
mj = 1.898*10**27


Np = 6

m = np.zeros(Np)
m[0] = M
m[1] = mm
m[2] = mv
m[3] = me
m[4] = mma
m[5] = mj

re = 6371*10**3
rv = 6052*10**3
rs = 695700*10**3
rm = 2440*10**3
rmma = 3389*10**3
rj = 69.911*10**3

G = 6.67*10**(-11)
dt = 300000 #seems about right


#iterations
nt = 1000

x = np.linspace(-2*AU,2*AU)
y = np.linspace(-2*AU,2*AU)

xp = np.zeros((1+nt,Np))
yp = np.zeros((1+nt,Np))
axp = np.zeros((1+nt,Np))
ayp = np.zeros((1+nt,Np))
up = np.zeros((1+nt,Np))
vp = np.zeros((1+nt,Np))

#Earth
xp[0,3] = AU
vp[0,3] = 29800

#Venus
xp[0,2] = -108200000*10**3
vp[0,2] = -35000

#mercury
yp[0,1] = -46001200*10**3
up[0,1] = 47.362*10**3

#mars
yp[0,4] = 1.3814*AU
up[0,4] = -24.077*10**3

#jupiter
yp[0,5] = 778500000*10**3
up[0,5] = -13.07*10**3


for n in range(nt):
	#the equations
	#dxe/dt = ue
	#due/dt = Fx/m


	#ue[1] = ue[0]+dt*Fx/m
	#ve[1] = ve[0]+dt*Fy/m
	
	#Make acceleration
	for i in range(Np):
		for j in range(Np):
			if j != i:
				axp[n+1,i] = axp[n+1,i]+(-G*m[j]/((xp[n,i]-xp[n,j])**2+(yp[n,i]-yp[n,j])**2)**(1.5))*(xp[n,i]-xp[n,j])
				ayp[n+1,i] = ayp[n+1,i]+(-G*m[j]/((xp[n,i]-xp[n,j])**2+(yp[n,i]-yp[n,j])**2)**(1.5))*(yp[n,i]-yp[n,j])
				
				
	#for i in range(Np):
	#	up[n+1,i] = up[n,i] + dt*axp[n+1,i]
	#	vp[n+1,i] = vp[n,i] + dt*ayp[n+1,i]
	up[n+1,:] = up[n,:] + dt*axp[n+1,:]
	vp[n+1,:] = vp[n,:] + dt*ayp[n+1,:]	
	
	xp[n+1,:] = xp[n,:]+dt*up[n+1,:]
	yp[n+1,:] = yp[n,:]+dt*vp[n+1,:]
	
	
	#for i in range(Np):
		#Euler Forward Method
		#xp[n+1,i] = xp[n,i]+dt*up[n+1,i]
		#yp[n+1,i] = yp[n,i]+dt*vp[n+1,i]

		
		#Verlet integration??
		#xp[n+1,i] = 2*xp[n,i]-xp[n-1,i]+dt**2*axp[n,i]
		#yp[n+1,i] = 2*yp[n,i]-yp[n-1,i]+dt**2*ayp[n,i]
		
		#Average of timesteps
		#xp[n+1,i] = xp[n,i]+dt*(0.5*up[n+1,i]+0.5*up[n,i])
		#yp[n+1,i] = yp[n,i]+dt*(0.5*vp[n+1,i]+0.5*vp[n,i])
	
	
plt.hold(True)
fig = plt.figure()
ax = fig.gca()
ax.set_xlim(-2*AU,2*AU)
ax.set_ylim(-2*AU,2*AU)
ax.set_title('Kappa')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')


def animate(i): #i increment with 1 each step
	ax.clear()
	ax.set_xlim(-6*AU,6*AU)
	ax.set_ylim(-6*AU,6*AU)
	ax.set_title('Kappa')
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	
	colors = ['yellow','black','black','blue','red','black']
	for j in range(Np):
		ax.plot(xp[0:i,j], yp[0:i,j],c=colors[j])
		ax.scatter(xp[i,j], yp[i,j],s=30,c=colors[j])
	
	if i == 50:
		plt.savefig('Orbit.png', bbox_inches='tight')
	return None


anim = animation.FuncAnimation(fig, animate, frames=nt, interval=100)



plt.show()


