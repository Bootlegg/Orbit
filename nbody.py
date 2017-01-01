import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.animation as animation
import mainr



class Planet:
	'Planets and objects including Sun'
	PlanetCount = 0
	planets = []


	def __init__(self,
				mass, 	#mass of planet
				x0, 	# initial x starting position
				y0,		# initial y starting position
				u0,		# initial u velocity
				v0,		# initial v velocity
				radius	#radius of planet, not used thus far
				):
	
		self.mass = mass
		self.radius = radius
		self.x0 = x0
		self.y0 = y0
		self.u0 = u0
		self.v0 = v0

		Planet.PlanetCount += 1
		Planet.planets.append(self)



	def posvos(nt):
		#Appending to arrays i numpy er ikke godt, og kræver re-allocation af hele array og 
		#temporarily doubles memory requirement... så bedst at have size of array til at starte med
		
		#Position
		Planet.xp = np.zeros((nt+1,Planet.PlanetCount))
		Planet.yp = np.zeros((nt+1,Planet.PlanetCount))
		
		#Velocity
		Planet.up = np.zeros((nt+1,Planet.PlanetCount))
		Planet.vp = np.zeros((nt+1,Planet.PlanetCount))

		#Acceleration
		Planet.axp = np.zeros((nt+1,Planet.PlanetCount))
		Planet.ayp = np.zeros((nt+1,Planet.PlanetCount))

		#Energies
		Planet.Tkin = np.zeros((1+nt,Planet.PlanetCount))
		Planet.Vpot = np.zeros((1+nt,Planet.PlanetCount))
		Planet.Tavg = np.zeros((1+nt,Planet.PlanetCount))
		Planet.Vavg = np.zeros((1+nt,Planet.PlanetCount))


		Planet.m = np.zeros(Planet.PlanetCount)

	def setposmass():
		for i in range(Planet.PlanetCount):
			Planet.xp[0,i] = Planet.planets[i].x0
			Planet.yp[0,i] = Planet.planets[i].y0
			Planet.up[0,i] = Planet.planets[i].u0
			Planet.vp[0,i] = Planet.planets[i].v0
			Planet.m[i] = Planet.planets[i].mass

def CallFortran():
	(Planet.axp,Planet.ayp,
	Planet.up,Planet.vp,
	Planet.xp,Planet.yp,
	Planet.Tavg,Planet.Vavg,
	Planet.Tkin,Planet.Vpot) = np.array(mainr.mainr(
									Planet.axp,Planet.ayp,
									Planet.up,Planet.vp,
									Planet.xp,Planet.yp,
									Planet.m,
									Planet.Tavg,Planet.Vavg,
									Planet.Tkin,Planet.Vpot,
									dt,G,nt,Np))


#instance Planet(mass,x,y,u,vradius)
#x,y,u,v are initial conditions for the motion and position
Sun = Planet(2.0*10**30,0,0,0,0,695700*10**3)

Mercury = Planet(3.3*10**23,0,-46001200*10**3,47.362*10**3,0,2440*10**3)

Venus = Planet(4.87*10**24,-108200000*10**3,0,0,-35000,6052*10**3)

Earth = Planet(6.0*10**24,1.5*10**11,0,0,29800,6371*10**3)

Mars = Planet(6.417*10**23,0,1.3814*1.5*10**11,-24.077*10**3,0,3389*10**3)

Jupiter = Planet(1.898*10**27,0,778500000*10**3,-13.07*10**3,0,69.911*10**3)


print("Number of planets = {}".format(Planet.PlanetCount))
print("Planet masses")
for i in range(Planet.PlanetCount):
	print(Planet.planets[i].mass)



#units
AU = 1.5*10**11
Np = Planet.PlanetCount

G = 6.67*10**(-11)
dt = 300000 #300000 #seems about right

#iterations
nt = 1000


Planet.posvos(nt)
Planet.setposmass()
CallFortran()

def TVplot():
	fig1 = plt.figure(1)
	ax1 = fig1.gca()
	#ax1.plot(T,c='red')
	#ax1.plot(V,c='blue')
	ax1.plot(Planet.Tavg[1:-2,0],c='red')
	ax1.plot(-0.5*Planet.Vavg[1:-2,0],c='blue')
	ax1.set_title('Virial Theorem')
	ax1.set_xlabel('timestep n')
	ax1.set_ylabel('Energy[J]')
	ax1.legend(['T', '-1/2 V'])
	fig1.savefig('figs/Virial.png', bbox_inches='tight')

TVplot()

plt.hold(True)
fig = plt.figure(2)
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
		ax.plot(Planet.xp[0:i,j], Planet.yp[0:i,j],c=colors[j])
		ax.scatter(Planet.xp[i,j], Planet.yp[i,j],s=30,c=colors[j])
	
	if i == 50:
		fig.savefig('figs/Orbit.png', bbox_inches='tight')
	return None


anim = animation.FuncAnimation(fig, animate, frames=nt, interval=100)

plt.show()


