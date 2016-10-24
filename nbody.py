import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.animation as animation
import add
import accel
import velocity



class Planet:
	'Planets and objects including Sun'
	PlanetCount = 0
	planets = []
	m = []

	def __init__(self,mass,x0,y0,u0,v0,radius):
		#x0,y0,u0,v0 are initial conditions for motion
		self.mass = mass
		self.radius = radius
		self.x0 = x0
		self.y0 = y0
		self.u0 = u0
		self.v0 = v0
		Planet.PlanetCount += 1
		Planet.planets.append(self)
		Planet.m.append(self.mass)




#instance Planet(mass,x,y,u,vradius)
#x,y,u,v are initial conditions for the motion and position
Sun = Planet(2*10**30,0,0,0,0,695700*10**3)

Mercury = Planet(3.3*10**23,0,-46001200*10**3,47.362*10**3,0,2440*10**3)

Venus = Planet(4.87*10**24,-108200000*10**3,0,0,-35000,6052*10**3)

Earth = Planet(6*10**24,1.5*10**11,0,0,29800,6371*10**3)

Mars = Planet(6.417*10**23,0,1.3814*1.5*10**11,-24.077*10**3,0,3389*10**3)

Jupiter = Planet(1.898*10**27,0,778500000*10**3,-13.07*10**3,0,69.911*10**3)


print("Number of planets = {}".format(Planet.PlanetCount))
print("Planet masses")
for i in range(Planet.PlanetCount):
	print(Planet.planets[i].mass)

m = np.array(Planet.m)


#print(P0.name)

# planets[0] = Sun
# planets[1] = mercury
# for i in range(Planet.PlanetCount):
# 	planets[i] = 


#units
AU = 1.5*10**11
Np = Planet.PlanetCount

G = 6.67*10**(-11)
dt = add.num(56)[0]#300000 #seems about right
#print(add.num()[0])
#print(add.num()[1])
print("Done with add, now for accel")

#iterations
nt = 1000

x = np.zeros((1+nt,Np))
y = np.zeros((1+nt,Np))
ax = np.zeros((1+nt,Np))
ay = np.zeros((1+nt,Np))
u = np.zeros((1+nt,Np))
v = np.zeros((1+nt,Np))
T = np.zeros(1+nt)
V = np.zeros(1+nt)

Tavg = np.zeros(1+nt)
Vavg = np.zeros(1+nt)

for i in range(Np):
	x[0,i] = Planet.planets[i].x0
	y[0,i] = Planet.planets[i].y0
	u[0,i] = Planet.planets[i].u0
	v[0,i] = Planet.planets[i].v0

#Initial kinetic and potential energy
T[0] = np.array([0.5*m[i]*(u[0,i]**2+v[0,i]**2) for i in range(Np)]).sum()

for i in range(Np):
	for j in range(i+1,Np):
		V[0] = V[0] + (-G*m[i]*m[j]/np.sqrt((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2))

#Basic stormer-verlet integration, first timestep
#Accelerations[0]
for i in range(Np):
		for j in range(Np):
			if j != i:
				ax[0,i] = ax[0,i]+(-G*m[j]/((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2)**(1.5))*(x[0,i]-x[0,j])
				ay[0,i] = ay[0,i]+(-G*m[j]/((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2)**(1.5))*(y[0,i]-y[0,j])
#positions step 1
x[1,:] = x[0,:] + u[0,:]*dt+0.5*ax[0,:]*dt**2 
y[1,:] = y[0,:] + v[0,:]*dt+0.5*ay[0,:]*dt**2 


print("x")
print(np.array(accel.accel(G,Np,m,x[0,:],y[0,:]))[0])
print(ax[0,:])

print("y")
print(np.array(accel.accel(G,Np,m,x[0,:],y[0,:]))[1])
print(ay[0,:])

for n in range(1,nt):
	#the equations
	#dxe/dt = ue
	#due/dt = Fx/m


	#ue[1] = ue[0]+dt*Fx/m
	#ve[1] = ve[0]+dt*Fy/m
	
	#Acceleration
	for i in range(Np):
		for j in range(Np):
			if j != i: #Here, either do n or n+1
				ax[n,i] = ax[n,i]+(-G*m[j]/((x[n,i]-x[n,j])**2+(y[n,i]-y[n,j])**2)**(1.5))*(x[n,i]-x[n,j])
				ay[n,i] = ay[n,i]+(-G*m[j]/((x[n,i]-x[n,j])**2+(y[n,i]-y[n,j])**2)**(1.5))*(y[n,i]-y[n,j])
				
				
	#update velocity
	#velocity either n or n+1, n+1 works with Euler Forward
	u[n,:] = u[n-1,:] + dt*ax[n,:]
	v[n,:] = v[n-1,:] + dt*ay[n,:]	
	#update positions
	#x[n+1,:] = x[n,:]+dt*u[n+1,:]
	#y[n+1,:] = y[n,:]+dt*v[n+1,:]

	#update positions verlet
	x[n+1,:] = 2*x[n,:]-x[n-1,:]+ax[n,:]*dt**2
	y[n+1,:] = 2*y[n,:]-y[n-1,:]+ay[n,:]*dt**2
	
	#Kinetic energy either do n+1 or n...
	T[n] = np.array([0.5*m[i]*(u[n,i]**2+v[n,i]**2) for i in range(Np)]).sum()
	#Potential energy
	for i in range(Np):
		for j in range(i+1,Np):
			V[n] = V[n] + (-G*m[i]*m[j]/np.sqrt((x[n,i]-x[n,j])**2\
							+ (y[n,i]-y[n,j])**2))

	#time average kinetic energy
	Tavg[n] = np.array(T[:n]).sum()*dt/((n)*dt)

	#time average potential energy
	Vavg[n] = np.array(V[:n]).sum()*dt/((n)*dt)


	##Different schemes
		#Average of timesteps
		#x[n+1,i] = x[n,i]+dt*(0.5*u[n+1,i]+0.5*u[n,i])
		#y[n+1,i] = y[n,i]+dt*(0.5*v[n+1,i]+0.5*v[n,i])
# fig1 = plt.figure(1)
# ax1 = fig1.gca()
# #ax1.plot(T,c='red')
# #ax1.plot(V,c='blue')
# ax1.plot(Tavg[1:-2],c='red')
# ax1.plot(-0.5*Vavg[1:-2],c='blue')
# ax1.set_title('Virial Theorem')
# ax1.set_xlabel('timestep n')
# ax1.set_ylabel('Energy[J]')
# ax1.legend(['T', '-1/2 V'])
# fig1.savefig('Virial.png', bbox_inches='tight')

# plt.hold(True)
# fig = plt.figure(2)
# ax = fig.gca()
# ax.set_xlim(-2*AU,2*AU)
# ax.set_ylim(-2*AU,2*AU)
# ax.set_title('Kappa')
# ax.set_ylabel('y [m]')
# ax.set_xlabel('x [m]')


# def animate(i): #i increment with 1 each step
# 	ax.clear()
# 	ax.set_xlim(-6*AU,6*AU)
# 	ax.set_ylim(-6*AU,6*AU)
# 	ax.set_title('Kappa')
# 	ax.set_ylabel('y [m]')
# 	ax.set_xlabel('x [m]')
	
# 	colors = ['yellow','black','black','blue','red','black']
# 	for j in range(Np):
# 		ax.plot(x[0:i,j], y[0:i,j],c=colors[j])
# 		ax.scatter(x[i,j], y[i,j],s=30,c=colors[j])
	
# 	if i == 50:
# 		fig.savefig('Orbit.png', bbox_inches='tight')
# 	return None


#anim = animation.FuncAnimation(fig, animate, frames=nt, interval=100)



#plt.show()


