import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.animation as animation
from matplotlib.ticker import FormatStrFormatter


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
Sun = Planet(2*10**30,0,-695700*10**3/2,347850/(365*24*3.6),0,695700*10**3)

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

tfaktor = 1
dt = tfaktor*300000#300000 #seems about right
year = 365.20*24*60*60

timeratio = dt/year


#print(add.num()[0])
#print(add.num()[1])

#====================================
#Print critical time #1.5*10**11 778500000*10**3
print("Largest critical time, using R = Mercury distance to sun")
tcr = np.sqrt(G*np.sum(m)/(46001200*10**3)**3)
print(tcr)
print("Timesteps to reach dynamic equilibrium")
print(tcr/dt)


print("Smallest critical time, using R = Jupiter distance to sun")
tcr = np.sqrt(G*np.sum(m)/(778500000*10**3)**3)
print(tcr)
print("Timesteps to reach dynamic equilibrium")
print(tcr/dt)


print("Done with add, now for accel")

#iterations
nt = 10000


#========================================================================
#Scaled nbody units
M0 = np.sum(m)
R0 = 1.5*10**11
T0 = np.sqrt(R0**3/(G*M0))
V0 = R0/T0
A0 = R0/(T0*T0)


ms = m/M0
dts = dt/T0
Gs = G*M0*T0*T0/R0**3

xs = np.zeros((1+nt,Np))/R0
ys = np.zeros((1+nt,Np))/R0
axs = np.zeros(Np)/A0
ays = np.zeros(Np)/A0
axs1 = np.zeros(Np)/A0
ays1 = np.zeros(Np)/A0
us = np.zeros((1+nt,Np))/V0
vs = np.zeros((1+nt,Np))/V0
Ts = np.zeros(1+nt)
Vs = np.zeros(1+nt)

Lzs = np.zeros(1+nt) #Angular momentum, only in z direction
Pxs = np.zeros(1+nt) #linear momentum
Pys = np.zeros(1+nt)

Tavgs = np.zeros(1+nt)
Vavgs = np.zeros(1+nt)

es = np.zeros(1+nt) #eccentricity 
e = np.zeros(1+nt)

for i in range(Np):
	#x[0,i] = Planet.planets[i].x0
	#y[0,i] = Planet.planets[i].y0
	#u[0,i] = Planet.planets[i].u0
	#v[0,i] = Planet.planets[i].v0
	
	xs[0,i] = Planet.planets[i].x0/R0
	ys[0,i] = Planet.planets[i].y0/R0
	us[0,i] = Planet.planets[i].u0/V0
	vs[0,i] = Planet.planets[i].v0/V0


#====================================================================
print("Time unit")
print(2*np.sqrt(2))
print(T0)




#=============================================================
#Initial linear momentum
#for i in range(Np):
#	print(u[0,i]*m[i])
#	Px[0] += u[0,i]*m[i]
#	Py[0] += v[0,i]*m[i]
print("Initial Linear Momentum Pxs")
for i in range(Np):
	print(us[0,i]*ms[i])
	Pxs[0] += us[0,i]*ms[i]
	Pys[0] += vs[0,i]*ms[i]

#=============================================================
#Initial angular momentum
print("Initial Angular Momentum")
for i in range(Np):
	#Lz[0] = m[i]*(x[0,i]*v[0,i]-y[0,i]*u[0,i])
	Lzs[0] = ms[i]*(xs[0,i]*vs[0,i]-ys[0,i]*us[0,i])
	print(Lzs[0])
	
	
#============================================================
#Initial eccentricity
#mu = 2*10**30*6*10**24/(2*10**30+6*10**24)

#mu = 12*10**30/(2*10**6+6)

mu = ms[0]*ms[3]/(ms[0]+ms[3])
#e[0] = np.sqrt(1+2*Lzs[0]**2*(Ts[0]+Vs[0])/(mu*G**2))
Learth = ms[3]*(xs[0,3]*vs[0,3]-ys[0,3]*us[0,3])
Eearth = 0.5*ms[3]*(us[0,3]**2+vs[0,3]**2) - ms[0]*ms[3]/np.sqrt((xs[0,3]-xs[0,0])**2+(ys[0,3]-ys[0,0])**2)
e[0] = np.sqrt(1+2*Learth**2*Eearth/mu)

#=========================================================================
#Timevariable
dtmin = dts




#=============================================================
#Start main integraton loop!

def Leapfrog(axs,ays,us,vs,Ts,Vs,Tavgs,Vavgs,dts):
	###
	#Second order leap frog, symplectic
	###

	for n in range(0,nt):
		#Calculate a0x,a0y
		axs[:] = 0
		ays[:] = 0
		#Acceleration
		for i in range(Np):
			for j in range(Np):
				if j != i: #Here, either do n or n+1
					axs[i] += -ms[j]/((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)**(1.5)*(xs[n,i]-xs[n,j])
					ays[i] += -ms[j]/((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)**(1.5)*(ys[n,i]-ys[n,j])
		
		#Position
		xs[n+1,:] = xs[n,:] + us[n,:]*dts + 0.5*axs*dts*dts
		ys[n+1,:] = ys[n,:] + vs[n,:]*dts + 0.5*ays*dts*dts
		
		#Velocity requires acceleration at x[n+1], so first we calculate the new positions
		#Calculate a1x,a1y
		axs1[:] = 0
		ays1[:] = 0
		#Acceleration
		for i in range(Np):
			for j in range(Np):
				if j != i: #Here, either do n or n+1
					axs1[i] += -ms[j]/((xs[n+1,i]-xs[n+1,j])**2+(ys[n+1,i]-ys[n+1,j])**2)**(1.5)*(xs[n+1,i]-xs[n+1,j])
					ays1[i] += -ms[j]/((xs[n+1,i]-xs[n+1,j])**2+(ys[n+1,i]-ys[n+1,j])**2)**(1.5)*(ys[n+1,i]-ys[n+1,j])
		
		#Velocity
		us[n+1,:] = us[n,:]+0.5*(axs1+axs)*dts
		vs[n+1,:] = vs[n,:]+0.5*(ays1+ays)*dts
		
		
		#==============================================================
		#Kinetic and potential energy, nbody units
		Ts[n] = np.array([0.5*ms[i]*(us[n,i]**2+vs[n,i]**2) for i in range(Np)]).sum()
		#Potential energy
		for i in range(Np):
			for j in range(i+1,Np):
				Vs[n] += (-ms[i]*ms[j]/np.sqrt((xs[n,i]-xs[n,j])**2\
								+ (ys[n,i]-ys[n,j])**2))
		
		
		#time average kinetic energy, removed dts/dts faktor
		Tavgs[n] = np.array(Ts[:n]).sum()/(n)

		#time average potential energy, removed dts/dts faktor
		Vavgs[n] = np.array(Vs[:n]).sum()/(n)	

		#=============================================================
		#linear momentum
		for i in range(Np):
			#print(u[n,i]*m[i])
			#Px[n] += u[n,i]*m[i]
			#Py[n] += v[n,i]*m[i]
			
			Pxs[n] += us[n,i]*ms[i]
			Pys[n] += vs[n,i]*ms[i]
		
		#=============================================================
		#Initial angular momentum
		for i in range(Np):
			#Lz[n] += m[i]*(x[n,i]*v[n,i]-y[n,i]*u[n,i])
			
			Lzs[n] += ms[i]*(xs[n,i]*vs[n,i]-ys[n,i]*us[n,i])
			
			
		#============================================================
		#eccentricity of earth, each timestep
		#e[n] = np.sqrt(1+2*Lz[n]**2*(T[n]+V[n])/(mu*G**2))
		#e[n] = np.sqrt(1+2*Lzs[n]**2*(Ts[n]+Vs[n])/mu)
		
		
		Learth = ms[3]*(xs[n,3]*vs[n,3]-ys[n,3]*us[n,3])
		Eearth = 0.5*ms[3]*(us[n,3]**2+vs[n,3]**2) - ms[0]*ms[3]/np.sqrt((xs[n,3]-xs[n,0])**2+(ys[n,3]-ys[n,0])**2)
		e[n] = np.sqrt(1+2*Learth**2*Eearth/mu)
		
		
		
		#=============================================================
		#Theoretical timestep
		#Can put this with acceleration
		dtmin = dts
		for i in range(Np):
			for j in range(Np):
				if j != i:
					dtij = (np.sqrt((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)
					/np.sqrt((us[n,i]-us[n,j])**2+(vs[n,i]-vs[n,j])**2))
					if dtij < dtmin:
						#print("Smaller timestep found")
						dtmin = dtij
		dts = dtmin

	
	
def StormerVerlet(axs,ays,us,vs,Ts,Vs,Tavgs,Vavgs,dts,dtmin):

	#=====================================================================
	#Basic stormer-verlet integration, first timestep
	#Accelerations[0]
	for i in range(Np):
			for j in range(Np):
				if j != i:
					#ax[0,i] = ax[0,i]+(-G*m[j]/((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2)**(1.5))*(x[0,i]-x[0,j])
					#ay[0,i] = ay[0,i]+(-G*m[j]/((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2)**(1.5))*(y[0,i]-y[0,j])
					
					axs[i] += (-ms[j]/((xs[0,i]-xs[0,j])**2+(ys[0,i]-ys[0,j])**2)**(1.5))*(xs[0,i]-xs[0,j])
					ays[i] += (-ms[j]/((xs[0,i]-xs[0,j])**2+(ys[0,i]-ys[0,j])**2)**(1.5))*(ys[0,i]-ys[0,j])
					
	#positions step 1
	#x[1,:] = x[0,:] + u[0,:]*dt + 0.5*ax[0,:]*dt**2 
	#y[1,:] = y[0,:] + v[0,:]*dt + 0.5*ay[0,:]*dt**2 



	xs[1,:] = xs[0,:] + us[0,:]*dts+0.5*axs[:]*dts**2 
	ys[1,:] = ys[0,:] + vs[0,:]*dts+0.5*ays[:]*dts**2 

	#=============================================================
	#Initial kinetic and potential energy
	#T[0] = np.array([0.5*m[i]*(u[0,i]**2+v[0,i]**2) for i in range(Np)]).sum()

	#for i in range(Np):
	#	for j in range(i+1,Np):
	#		V[0] += (-G*m[i]*m[j]/np.sqrt((x[0,i]-x[0,j])**2+(y[0,i]-y[0,j])**2))
			
	
	Ts[0] = np.array([0.5*ms[i]*(us[0,i]**2+vs[0,i]**2) for i in range(Np)]).sum()
	print("1/|rsi-rsj| terms in potential energy")
	for i in range(Np):
		for j in range(i+1,Np):
			Vs[0] += (-ms[i]*ms[j]/np.sqrt((xs[0,i]-xs[0,j])**2+(ys[0,i]-ys[0,j])**2))
			print(1/np.sqrt((xs[0,i]-xs[0,j])**2+(ys[0,i]-ys[0,j])**2))
	print("Initial Potential Energy Vs")
	print(Vs[0])
		
	
	for n in range(1,nt):
		#the equations
		#dxe/dt = ue
		#due/dt = Fx/m
		
		
		#ue[1] = ue[0]+dt*Fx/m
		#ve[1] = ve[0]+dt*Fy/m
		axs[:] = 0
		ays[:] = 0
		#Acceleration
		for i in range(Np):
			for j in range(Np):
				if j != i: #Here, either do n or n+1
					axs[i] += -ms[j]/((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)**(1.5)*(xs[n,i]-xs[n,j])
					ays[i] += -ms[j]/((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)**(1.5)*(ys[n,i]-ys[n,j])
		#update velocity
		#velocity either n or n+1, n+1 works with Euler Forward

		us[n,:] = us[n-1,:] + dts*axs[:]
		vs[n,:] = vs[n-1,:] + dts*ays[:]	
		
		#update positions
		#x[n+1,:] = x[n,:]+dt*u[n+1,:]
		#y[n+1,:] = y[n,:]+dt*v[n+1,:]

		#update positions verlet

		
		xs[n+1,:] = 2*xs[n,:]-xs[n-1,:]+axs[:]*dts**2
		ys[n+1,:] = 2*ys[n,:]-ys[n-1,:]+ays[:]*dts**2
		
					
							
		#==============================================================
		#Kinetic and potential energy, nbody units
		Ts[n] = np.array([0.5*ms[i]*(us[n,i]**2+vs[n,i]**2) for i in range(Np)]).sum()
		#Potential energy
		for i in range(Np):
			for j in range(i+1,Np):
				Vs[n] += (-ms[i]*ms[j]/np.sqrt((xs[n,i]-xs[n,j])**2\
								+ (ys[n,i]-ys[n,j])**2))
		
		
		#time average kinetic energy, removed dts/dts faktor
		Tavgs[n] = np.array(Ts[:n]).sum()/(n)

		#time average potential energy, removed dts/dts faktor
		Vavgs[n] = np.array(Vs[:n]).sum()/(n)	

		#=============================================================
		#linear momentum
		for i in range(Np):
			#print(u[n,i]*m[i])
			#Px[n] += u[n,i]*m[i]
			#Py[n] += v[n,i]*m[i]
			
			Pxs[n] += us[n,i]*ms[i]
			Pys[n] += vs[n,i]*ms[i]
		
		#=============================================================
		#Initial angular momentum
		for i in range(Np):
			#Lz[n] += m[i]*(x[n,i]*v[n,i]-y[n,i]*u[n,i])
			
			Lzs[n] += ms[i]*(xs[n,i]*vs[n,i]-ys[n,i]*us[n,i])
			
			
		#============================================================
		#eccentricity of earth, each timestep
		#e[n] = np.sqrt(1+2*Lz[n]**2*(T[n]+V[n])/(mu*G**2))
		#e[n] = np.sqrt(1+2*Lzs[n]**2*(Ts[n]+Vs[n])/mu)
		
		
		Learth = ms[3]*(xs[n,3]*vs[n,3]-ys[n,3]*us[n,3])
		Eearth = 0.5*ms[3]*(us[n,3]**2+vs[n,3]**2) - ms[0]*ms[3]/np.sqrt((xs[n,3]-xs[n,0])**2+(ys[n,3]-ys[n,0])**2)
		e[n] = np.sqrt(1+2*Learth**2*Eearth/mu)
		
		
		
		#=============================================================
		#Theoretical timestep
		#Can put this with acceleration
		dtmin = dts
		for i in range(Np):
			for j in range(Np):
				if j != i:
					dtij = (np.sqrt((xs[n,i]-xs[n,j])**2+(ys[n,i]-ys[n,j])**2)
					/np.sqrt((us[n,i]-us[n,j])**2+(vs[n,i]-vs[n,j])**2))
					if dtij < dtmin:
						#print("Smaller timestep found")
						dtmin = dtij
		dts = dtmin
		##Different schemes
			#Average of timesteps
			#x[n+1,i] = x[n,i]+dt*(0.5*u[n+1,i]+0.5*u[n,i])
			#y[n+1,i] = y[n,i]+dt*(0.5*v[n+1,i]+0.5*v[n,i])

StormerVerlet(axs,ays,us,vs,Ts,Vs,Tavgs,Vavgs,dts,dtmin)
#Leapfrog(axs,ays,us,vs,Ts,Vs,Tavgs,Vavgs,dts)

#print(Px[0:20])
#print(Px[-40:-2])
print("timesteps")
print(dts)
print(dtmin)
#=============================================================
#Virial plot
fig1 = plt.figure(1)
#ax1 = fig1.gca()
ax1 = fig1.add_subplot(121)
#ax1.plot(T,c='red')
#ax1.plot(V,c='blue')
ax1.plot(Tavgs[1:-2]/10**(-5),c='red')
ax1.plot(-0.5*Vavgs[1:-2]/10**(-5),c='blue')
ax1.set_title('Virial Theorem')
ax1.set_xlabel('timestep n')
ax1.set_ylabel('Energy / e-5')
ax1.legend(['T', '-1/2 V'])
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))



Etotal = Ts[:-2]+Vs[:-2]
#ticks = np.linspace(min(Etotal),max(Etotal),5)

#Total energy plot
#fig2 = plt.figure(2)
#ax2 = fig2.gca()
ax2 = fig1.add_subplot(122)
#ax2.plot(Etotal/10**(-5),color="black")
ax2.plot(Etotal,color="black")
ax2.set_xlabel("timstep n")
ax2.set_ylabel("Energy/e-5")
ax2.set_title("Total energy")
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax2.yaxis.set_label_position("right")
ax2.yaxis.tick_right()

fig1.savefig('Virial.png', bbox_inches='tight')

print("Kinetic energy")
print(Tavgs[1:15])
print("Potential energy")
print(Vavgs[1:15])

print("Potential energy from nbody units")
print(-Gs*M0**2/R0)


print("Gravitational Nbody constant, should be Gs = 1?")
print(Gs)

#=============================================================
#Total angular momentum
fig3 = plt.figure(3)
ax3 = fig3.gca()
#ax3 = fig3.add_subplot(121)
ax3.plot(Lzs[1:-2]/10**(-3),color="black")
#ax3.plot(Py[:-2],color="red")
ax3.set_xlabel("timestep n")
ax3.set_ylabel("Angular Momentum/e-3")
ax3.set_title("Total Angular Momentum")
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#
fig3.savefig('AngularMomentum.png', bbox_inches='tight')

#Linear momentum
#fig3 = plt.figure(4)
#ax3 = fig3.gca()
# ax3 = fig3.add_subplot(122)
# ax3.plot(Pxs[1:-2]/10**(-5),color="black")#
# #ax3.plot(Pys[:-2],color="red")
# ax3.set_xlabel("timestep n")
# ax3.set_ylabel("Momentum/e-2")
# ax3.set_title("Total Linear Momentum")
# ax3.axis([0,nt,min(Pxs[1:-2]),max(Pxs[1:-2])])
# ax3.yaxis.set_label_position("right")
# ax3.yaxis.tick_right()
#ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#===============================================================
#analyze sun motion? or maybe earth? let's start with earth
#should be same periods no matter what we do with it, right? Maybe?
plt.figure(6)
plt.plot(e[:-2],color="black")
plt.xlabel("timestep n")
plt.ylabel("e eccentricity")
plt.title("Eccentricity of Earth-Sun")
plt.axis([0,nt,min(e[:-2]),max(e[:-2])])

print("Eccentricities of Earth-Sun")
print(e[-20:-2])

#==============================================================
#Sun velocity... måske lidt bedre end position, fordi, velocity får man fra rødforskydning etc når solen kommer tættere/nærmere
#Btw, vi har måske ikke hele period fra Jupiter med?!
plt.figure(7)
plt.plot(xs[:-2,0],color="black")
plt.plot(ys[:-2,0],color="black")
plt.xlabel("timestep n")
#plt.ylabel("v(t) m/s")
plt.title("Periodic data from Sun")


#DFT af sun data v[:-2,0]

#Jeg mangler vel lidt med dt = 3000000 etc? Der skal lidt mere info ind i denne fourer transform i think


def DFT():
	print("Doing DFT")
	k = 10 #siger vi lige
	Xk = 0
	#Vdat1 = vs[:-2,0]
	Vdat1 = ys[:-2,0]
	#Vdat1 = Vdat1[]
	Ndat = len(Vdat1)
	vdat = np.array(Vdat1,dtype=complex)
	Xkarray = np.zeros(Ndat,dtype=complex)

	pi2= 2.0*np.pi

	for k in range(Ndat):
		for n in range(Ndat):
			Xkarray[k] += vdat[n]*np.exp(-1j*pi2*k*n/Ndat)
	#print(Xk)
	Xkarray *= 1/Ndat
	print(Xkarray[0:10])

	Xkarray1 = np.fft.fft(Vdat1)

	plt.figure(8)
	plt.xlabel("k")
	plt.ylabel("Xk")
	plt.title("Fourier transform")
	#plt.plot(Xkarray1.imag,color="red")
	#plt.plot(Xkarray1.real,color="blue")
	#plt.axis([0,100,0,10])
	plt.plot(np.absolute(Xkarray),color="black")


#DFT()
	
def PlotOrbit():
	plt.figure(9)
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("Orbit")
	plt.plot(xs[0:-2,3], ys[0:-2,3],c="blue")
	plt.plot(xs[0:-2,0], ys[0:-2,0],c="black")
	
PlotOrbit()


fig = plt.figure(5)
ax = fig.gca()
#ax.set_xlim(-2*R0,2*R0)
#ax.set_ylim(-2*R0,2*R0)
ax.set_title('N-body units')
ax.set_ylabel('ys')
ax.set_xlabel('xs')


def animateUnits(i): #i increment with 1 each step
	ax.clear()
	ax.set_xlim(-6,6)
	ax.set_ylim(-6,6)
	ax.set_title('Kappa n={}, yrs={:0.1f}'.format(i,timeratio*i))
	ax.set_ylabel('ys')
	ax.set_xlabel('xs')
	
	colors = ['yellow','black','black','blue','red','black']
	for j in range(Np):
		ax.plot(xs[0:i,j], ys[0:i,j],c=colors[j])
		ax.scatter(xs[i,j], ys[i,j],s=30,c=colors[j])
	
	#if i == 50:
	#	fig.savefig('Orbit.png', bbox_inches='tight')
	return None

	
animUnits = animation.FuncAnimation(fig, animateUnits, frames=nt, interval=100)



plt.show()


