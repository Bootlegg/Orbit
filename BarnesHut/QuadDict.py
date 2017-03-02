"""
Quadtree

2017
Dan Krog
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation




Dict = {}


##############################
#Sæt plot dict og print dict her?
def PlotDict(level,rectdict):
	global ListPlot
	#print(level)
	
	#if level == maxlevel-1: #we're at 2 now, while the class is at level 1... that's because
	#the class will append the rect to the NEXT level, which is 2, thus we need to be at 2 for plotting
	#print(len(rectdict.keys()))

	level += 1
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			for i in range(4):
			
				#Here should maybe append *rectdict[Areas[j]]["A"][i] or rectdict[Areas[j]]["A"][i]
				#to a list, if it isn't already in
				if rectdict[area]["A"][i] not in ListPlot:
					ListPlot.append(rectdict[area]["A"][i])
					plt.plot(*rectdict[area]["A"][i], color = "black")
				#else:
					#pass
			

	#else: #Fjern den her
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			#Den kalder rectdict HVIS rectdict["A1"] etc findes, ellers lader den vær.
			#Så når man kommer dybt nok er der ikke flere "A1","A2" etc, så stopper den loop.

			#Et problem med den her er at alde sub-calls vil jo add til level += 1.
			#Så hvis jeg bruger while level < 2 fx... så bliver det buggy
			
			PlotDict(level,rectdict[area])

def PrintDict():
	global levuls, totalparticles
	#Lets print the dict
	print("PRINTING DICT")
	totalparticles = 0
	level = 0
	levuls = [0 for i in range(maxlevel)]

	levelarea = []
	PrintNumChildren(Dict,level,levelarea)
	print("Total number of counted particles")
	print(totalparticles)
	print("Total len of points")
	print(len(points))
	print("LET'S LOOK AT LEVELS")
	print(levuls)


def PrintNumChildren(rectdict,level,levela):
	global totalparticles, levuls
	#Need to add to this a level, maybe this is something like level = level, like for functional programming

	#if level == 5:
	for i in range(maxlevel):
		if level == i:
			levuls[i] += 1
	print("Current box:",levela)
	print(rectdict.keys())

	children = 0

	if "Partics" in rectdict:
		children = len(rectdict["Partics"])
	
	totalparticles = totalparticles + children
	


	print("level {} and {} children".format(level,children))
	if children != 0:
		print(rectdict["Partics"])
	
	

	
	#level += 1

	#Den første her siger bare, HVIS der er en deeper level, så laver vi level += 1
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			level +=1
			break
	#Her der ved vi allerede der er en deeper level, så vi prøver bare hver area...
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict: 
			levela.append(area)
			#print("Current box:",area)
			PrintNumChildren(rectdict[area],level,levela)
			levela.remove(area)




def NewSplit(level,lims,rectdict):
	#dictrect is a dictionary atm,
	#empty dict {}

	Lx0, Lx, Ly0, Ly = lims
	area1 = [	([Lx0,Lx0 + (Lx-Lx0)/2]				,[Ly0,Ly0]),
				([Lx0,Lx0 + (Lx-Lx0)/2]				,[Ly0 + (Ly-Ly0)/2,Ly0 + (Ly-Ly0)/2]),
				([Lx0,Lx0]							,[Ly0,Ly0 + (Ly-Ly0)/2]),
				([Lx0 + (Lx-Lx0)/2,Lx0 + (Lx-Lx0)/2],[Ly0,Ly0 + (Ly-Ly0)/2])]

	#A1coords = [Lx0,Lx0+(Lx-Lx0)/2,Ly0,Ly0 + (Ly-Ly0)/2]
	
	area2 = [	([Lx0 + (Lx-Lx0)/2,Lx]				,[Ly0,Ly0]),
				([Lx0 + (Lx-Lx0)/2,Lx]				,[Ly0 + (Ly-Ly0)/2,Ly0 + (Ly-Ly0)/2]),
				([Lx0 + (Lx-Lx0)/2,Lx0 + (Lx-Lx0)/2],[Ly0,Ly0 + (Ly-Ly0)/2]),
				([Lx,Lx]							,[Ly0,Ly0 + (Ly-Ly0)/2])]

	#A2coords = [Lx0+(Lx-Lx0)/2,Lx,Ly0,Ly0 + (Ly-Ly0)/2]

	area3 = [	([Lx0,Lx0 + (Lx-Lx0)/2]				,[Ly,Ly]),
				([Lx0,Lx0 + (Lx-Lx0)/2]				,[Ly0 + (Ly-Ly0)/2,Ly0 + (Ly-Ly0)/2]),
				([Lx0,Lx0]							,[Ly0 + (Ly-Ly0)/2,Ly]),
				([Lx0 + (Lx-Lx0)/2,Lx0 + (Lx-Lx0)/2],[Ly0 + (Ly-Ly0)/2,Ly])]

	#A3coords = [Lx0,Lx0+(Lx-Lx0)/2,Ly0 + (Ly-Ly0)/2,Ly]

	area4 = [	([Lx0 + (Lx-Lx0)/2,Lx]				,[Ly,Ly]),
				([Lx0 + (Lx-Lx0)/2,Lx]				,[Ly0 + (Ly-Ly0)/2,Ly0 + (Ly-Ly0)/2]),
				([Lx0 + (Lx-Lx0)/2,Lx0 + (Lx-Lx0)/2],[Ly0 + (Ly-Ly0)/2,Ly]),
				([Lx,Lx]							,[Ly0 + (Ly-Ly0)/2,Ly])]

	#A4coords = [Lx0+(Lx-Lx0)/2,Lx,Ly0 + (Ly-Ly0)/2,Ly]			


	#Tror at lave A1 kan cancel en previous mere detailed A1, så kun hvis A1 allerede ikke eksistere.
	#Her tror jeg faktisk IKKE jeg skal append... det gør pointcoords for mig!
	#Så jeg prøver at outcomment disse point appends


	if "Partics" in rectdict:
		children = len(rectdict["Partics"])
		

		if children >= 2:

			if "A1" not in rectdict:
				rectdict["A1"] = {}

			if "A2" not in rectdict:
				rectdict["A2"] = {}

			if "A3" not in rectdict:
				rectdict["A3"] = {}

			if "A4" not in rectdict:
				rectdict["A4"] = {}

			rectdict["A1"]["A"] = area1
			rectdict["A2"]["A"] = area2
			rectdict["A3"]["A"] = area3
			rectdict["A4"]["A"] = area4


			rectdict["A1"]["r"] = (Lx-Lx0)/2
			rectdict["A2"]["r"] = (Lx-Lx0)/2
			rectdict["A3"]["r"] = (Lx-Lx0)/2
			rectdict["A4"]["r"] = (Lx-Lx0)/2
			
			#Doesnt seem to be a maxlevel problem!
			if level < maxlevel: 
				

				level += 1
				#Doesn't have to be append, really
				#ListParticles = []
				#for particle in rectdict["Particles"]:

				#	ListParticles.append(particle)

				#Removing all the particles from dict by setting list to empty list = []
				#rectdict["Particles"] = []
				#for particle in ListParticles:
				#	PointCoords(level,particle,lims,rectdict)
					#rectdict["Particles"].remove(particle)


				##################################################################
				#DO THE PARTICS DICT

				DictPartics = {}
				num = 0
				#for partic in rectdict["Partics"]:
				#	DictPartics[partic] = rectdict["Partics"][partic]
				DictPartics.update(rectdict["Partics"])

				#Gøre den empty, men hvad med, tbh, at fjerne "Partics" key completely?
				#Der er mange måder...
				#rectdict["Partics"].clear()
				rectdict.pop("Partics")
				#rectdict.pop('Partics', None)


				for partic in DictPartics:
					PointCoords(level,DictPartics[partic],partic,lims,rectdict)
					#pass

def PointCoords(level,point,pointkey,lims,rectdict):
	#Første gang jeg caller pointcoords tager den lims = [0,10,0,10]
	#Anden gang jeg caller pointcoords,hvis den lander i A1, tager den lims = [0,5,0,5]


	#JEg skal lave particles om til dictionaries på en måde,
	#Så det er mere overskueligt med timestepping etc...




	#if level <= maxlevel:
		#xdot,ydot = point
	xdot,ydot = point["x"],point["y"]#point[0],point[1]

	#Ser hvis point er i Quadrant 1,2,3 eller 4.
	if xdot <= lims[0]+(lims[1]-lims[0])/2:
		if ydot <= lims[2] + (lims[3]-lims[2])/2:
			#if "A1" not in rectdict:
			#	rectdict["A1"] = {}
			coords = [lims[0],lims[0]+(lims[1]-lims[0])/2,lims[2],lims[2] + (lims[3]-lims[2])/2]
			if "A1" in rectdict:
				#(lims[0],lims[1]/2,lims[2],lims[3]/2)
				level += 1
				PointCoords(level,point,pointkey,coords,rectdict["A1"])

				#Her, lav en children check faktisk?
				#Det er DENNE rectdict jeg gerne vil spille, det er ikke "A1"
				#Derfor skal jeg nok også bruge lims, IKKE coords.
				#Men det giver da ikke så meget mening at kalde NewSplit her?
				#Kun fordi det er en slags quickfik eller noget?:S

				#NewSplit(level,lims,rectdict)

				#Elelr også kan jeg kalde begge, bare for at være sikker
				#NewSplit(level,coords,rectdict["A1"])

				#for particle in rectdict["Particles"]:

					#rectdict["Particles"].remove(particle)

			else:
				#if "Particles" not in rectdict:
				#	rectdict["Particles"] = []

				#rectdict["Particles"].append(point)

				#Lav dictionary af particles
				if "Partics" not in rectdict:
					rectdict["Partics"] = {}

				ParticleNumber = len(rectdict["Partics"])

				rectdict["Partics"][pointkey] = point
				#rectdict["Partics"].update(point)
				# rectdict["Partics"]["P{}".format(ParticleNumber)] = {"x":point["x"],
				# 													"y":point["y"],
				# 													"m":point["m"],
				# 													"xold":0,
				# 													"yold":0,
				# 													"u":0,
				# 													"v":0,
				# 													"ax":0,
				# 													"ay":0}
				#lims eller coords
				NewSplit(level,lims,rectdict)


				#NewSplit(level,coords,rectdict["A1"])

				#Måske lave en call her? Men nok ikke helt rigtigt, tbh..
				#Lave en call til NewSplit? Hmm? Maybe not...
				# children = len(rectdict["Particles"])
				# if children >= 4:
				# #Doesnt seem to be a maxlevel problem!
				# 	if level < maxlevel+1000: 
				# 		NewSplit(level,coords,rectdict)

		else:
			coords = [lims[0],lims[0]+(lims[1]-lims[0])/2,lims[2]+(lims[3]-lims[2])/2,lims[3]]#(lims[0],lims[1]/2,lims[3]/2,lims[3])
			if "A3" in rectdict:
				level += 1
				#Med particle,
				PointCoords(level,point,pointkey,coords,rectdict["A3"])

				#NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A3"])

			else:
				#if "Particles" not in rectdict:
				#	rectdict["Particles"] = []
				#rectdict["Particles"].append(point)


				#Lav dictionary af particles
				if "Partics" not in rectdict:
					rectdict["Partics"] = {}

				ParticleNumber = len(rectdict["Partics"])
				rectdict["Partics"][pointkey] = point

				NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A3"])



	else:
		if ydot <= lims[2] + (lims[3]-lims[2])/2:
			coords = [lims[0]+(lims[1]-lims[0])/2,lims[1],lims[2],lims[2]+(lims[3]-lims[2])/2]#(lims[1]/2,lims[1],lims[2],lims[3]/2)
			if "A2" in rectdict:
				
				level += 1
				#Med particle,
				PointCoords(level,point,pointkey,coords,rectdict["A2"])

				#NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A2"])

			else:
				#if "Particles" not in rectdict:
				#	rectdict["Particles"] = []
				#rectdict["Particles"].append(point)

				#Lav dictionary af particles
				if "Partics" not in rectdict:
					rectdict["Partics"] = {}

				ParticleNumber = len(rectdict["Partics"])
				rectdict["Partics"][pointkey] = point
				
				NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A2"])

		else:
			coords = [lims[0]+(lims[1]-lims[0])/2,lims[1],lims[2]+(lims[3]-lims[2])/2,lims[3]]#(lims[1]/2,lims[1],lims[3]/2,lims[3])
			if "A4" in rectdict:
				
				level += 1
				#Med particle,
				PointCoords(level,point,pointkey,coords,rectdict["A4"])

				#NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A4"])

			else:
				#if "Particles" not in rectdict:
				#	rectdict["Particles"] = []
				#rectdict["Particles"].append(point)

				#Lav dictionary af particles
				if "Partics" not in rectdict:
					rectdict["Partics"] = {}

				ParticleNumber = len(rectdict["Partics"])
				rectdict["Partics"][pointkey] = point

				NewSplit(level,lims,rectdict)

				#NewSplit(level,coords,rectdict["A4"])

	#Ved ikke helt om det her skal med
	#else:
		#Jeg tror her skal jeg append i rectdict FØR, hmm...
	#	if "Particles" not in rectdict:
	#		rectdict["Particles"] = []
	#	rectdict["Particles"].append(point)


	#Lige nu så fjerner jeg 
	#NewSplit(level,lims,rectdict)

def CenterOfMass(rectdict):

	#EFTER man har sat points ind, skal man traverse tree, og calculate CM etc.
	#Så CenterOfMass
	#skal være independent, kaldes EFTER de andre.

	#Discuss Center Of Mass her, Multipole expansion etc her.

	#Det er basically, jeg kan measure distance d hen til hver CM fra en given particle...
	#Også hvis r/d < 0.5, så kan jeg bruge CM, right?


	if "Partics" in rectdict:
		children = len(rectdict["Partics"])


		#Remove CM, no particles.
		if children == 0:
			rectdict["CM"] = []
		else:
			M = 0
			CMx = 0
			CMy = 0
			for partic in rectdict["Partics"]:
				
				M += rectdict["Partics"][partic]["m"]
				
				CMx += rectdict["Partics"][partic]["x"]*rectdict["Partics"][partic]["m"]
				CMy += rectdict["Partics"][partic]["y"]*rectdict["Partics"][partic]["m"]

			CMx = CMx/M
			CMy = CMy/M


			rectdict["CM"] = [CMx,CMy, M]
			#print(rectdict["CM"])

			plt.scatter(CMx,CMy,c="red")


	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			CenterOfMass(rectdict[area])

def CalcForce(rectdict,level):
	#find a particle, by traversing THE ENTIRE dictionary,
	#Maybe have an identifier to note if "this section has no particles it all"..
	#OR remove sections or something...

	#So i have a particle, THEN i need to traverse entire dictionary again
	
	#I'm gonna need to add numbers ax,ay...

	#Kan det være jeg bør have flere functions, ligesom NewSplit og PointCoords?
	#Der arbejder sammen


	#Måske aller først, bare gør den classical direct sum,
	#og BAGEFTER gør det fancy med CM etc..
	#Så lige nu vil jeg bare gerne loop de loop alle particler


	#Det er jo sådan set kun de deepest levels som vil have particles
	#Level 0 vil ikke have particles når der er 4 eller mere

	#Så, når points>=4, så har vi, at dette IKKE gælder, der er ikke Partics in rectdict..
	if "Partics" in rectdict:
		if len(rectdict["Partics"]) > 0:
			
			

			for partic in rectdict["Partics"]:

				#traverse the whole dict again
				#Skal jeg lave Fx,Fy global måske? Eller carry dem ind i function som argument?

				#Eller måske skal man bruge noget return, og set Fx,Fy = FindParti
				#Så kan være return hjælper

				#Making accel elements, not a good way, if there's many timesteps,
				#then EACH time we will extend the particle list..
				#particle.extend([0,0])

				#Vi laver en carry for testing
				#rectdict["Partics"][partic]["Carry"] = 0

				#Use Particles
				#FindParticleForForce(Dict,rectdict["Partics"][partic])
				
				#Use Multipole Expansion
				FindAreaForForce(Dict, rectdict["Partics"][partic])
				
				#Making velocity...not a good way, if there's many timesteps,
				#then EACH time we will extend the particle list..
				#particle.extend([0,0])


				#Problem here is that I'm updating positions, but they will then be altered for the
				#next particle... so I need xold,x,xnew maybe, so that other particles still can use xold etc
				#Timestep velocity
				#UpdateVelocity(rectdict["Partics"][partic])

				#Timestep position
				UpdatePosition(rectdict["Partics"][partic])

				#print("Printing Carry")
				#print("level {}, carry {}".format(level,rectdict["Partics"][partic]["Carry"]))



			#print("Done with force")
			#print(rectdict["Partics"][partic])
	#Når den har været igennem alle particles her, så går vi et level deeper
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			#print("Going deeper")
			#level += 1
			CalcForce(rectdict[area],level)
			#level -= 1



def FindParticleForForce(rectdict,particle):#,Fx,Fy):
	#global Fx, Fy
	#print(particle)
	if "Partics" in rectdict:
			if len(rectdict["Partics"]) > 0:
				for particle1 in rectdict["Partics"]:

					#Her skal jeg måske passe på den ikke compare "P5" med "Jupiter" selvom de faktisk er ens?
					if particle != rectdict["Partics"][particle1]:

						particle["ax"] = particle["ax"] + (-1)*G*rectdict["Partics"][particle1]["m"]*((particle["x"]-rectdict["Partics"][particle1]["x"])
																/((particle["x"]-rectdict["Partics"][particle1]["x"])**2
																+(particle["y"]-rectdict["Partics"][particle1]["y"])**2)**(1.5))
						
						particle["ay"] = particle["ay"] + (-1)*G*rectdict["Partics"][particle1]["m"]*((particle["y"]-rectdict["Partics"][particle1]["y"])
																/((particle["x"]-rectdict["Partics"][particle1]["x"])**2
																+(particle["y"]-rectdict["Partics"][particle1]["y"])**2)**(1.5))

	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			
			#particle["Carry"] += 1
			FindParticleForForce(rectdict[area],particle)#,Fx,Fy)
	#else:
		#print(Fx,Fy)
		#particle.extend([Fx,Fy])
		
		#pass
		#return Fx,Fy

def FindAreaForForce(rectdict, particle):
	if "CM" in rectdict:
		d = np.sqrt((particle["x"]-rectdict["CM"][0])**2+(particle["y"]-rectdict["CM"][1])**2)
		r = rectdict["r"]

		r=2
		d=1
		if r/d < 0.5:
			particle["ax"] = particle["ax"] + (-1)*G*rectdict["CM"][2]*((particle["x"]-rectdict["CM"][0])
												/((particle["x"]-rectdict["CM"][0])**2
												+(particle["y"]-rectdict["CM"][1])**2)**(1.5))

			particle["ay"] = particle["ay"] + (-1)*G*rectdict["CM"][2]*((particle["y"]-rectdict["CM"][1])
												/((particle["x"]-rectdict["CM"][0])**2
												+(particle["y"]-rectdict["CM"][1])**2)**(1.5))
		else:
			for partic2 in rectdict["Partics"]:
				if particle != rectdict["Partics"][partic2]:
					particle["ax"] = particle["ax"] + (-1)*G*rectdict["Partics"][partic2]["m"]*((particle["x"]-rectdict["Partics"][partic2]["x"])
																	/((particle["x"]-rectdict["Partics"][partic2]["x"])**2
																	+(particle["y"]-rectdict["Partics"][partic2]["y"])**2)**(1.5))

					particle["ay"] = particle["ay"] + (-1)*G*rectdict["Partics"][partic2]["m"]*((particle["y"]-rectdict["Partics"][partic2]["y"])
																	/((particle["x"]-rectdict["Partics"][partic2]["x"])**2
																	+(particle["y"]-rectdict["Partics"][partic2]["y"])**2)**(1.5))

	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			FindAreaForForce(rectdict[area],particle)

def UpdateVelocity(particle):
	particle["u"] = particle["u"]+dt*particle["ax"]
	particle["v"] = particle["v"]+dt*particle["ay"]

def UpdatePosition(particle):
	particle["xnew"] = 2*particle["x"] - particle["xold"] + particle["ax"]*dt**2
	particle["ynew"] = 2*particle["y"] - particle["yold"] + particle["ay"]*dt**2



def CalcForceFirst(rectdict,level):
	if "Partics" in rectdict:
		if len(rectdict["Partics"]) > 0:
			print("level {} we have partics".format(level))
			for partic in rectdict["Partics"]:


				FindParticleForForce(Dict,rectdict["Partics"][partic])

				#UpdateVelocity(rectdict["Partics"][partic]):

				UpdatePositionFirst(rectdict["Partics"][partic])

	#Når den har været igennem alle particles her, så går vi et level deeper
	#Med 4>= particles, så vil ALLE 4 areas være i level 0
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			#print("Going deeper")
			level += 1
			print(area)
			CalcForceFirst(rectdict[area],level)
			level -= 1



def UpdatePositionFirst(particle):
	particle["xnew"] = particle["x"] + particle["u"]*dt + 0.5*particle["ax"]*dt**2
	particle["ynew"] = particle["y"] + particle["v"]*dt + 0.5*particle["ay"]*dt**2		




def UpdateOld(rectdict):
	"""
	After timestep is done, we update the old positions to the current.
	"""
	if "Partics" in rectdict:
		if len(rectdict["Partics"]) > 0:
			for partic in rectdict["Partics"]:

				#traverse the whole dict again
				#Skal jeg lave Fx,Fy global måske? Eller carry dem ind i function som argument?

				#Eller måske skal man bruge noget return, og set Fx,Fy = FindParti
				#Så kan være return hjælper

				#Making accel elements, not a good way, if there's many timesteps,
				#then EACH time we will extend the particle list..
				#particle.extend([0,0])
				
				
				rectdict["Partics"][partic]["xold"] = rectdict["Partics"][partic]["x"]
				rectdict["Partics"][partic]["yold"] = rectdict["Partics"][partic]["y"]
				
				rectdict["Partics"][partic]["x"] = rectdict["Partics"][partic]["xnew"]
				rectdict["Partics"][partic]["y"] = rectdict["Partics"][partic]["ynew"]
				
				rectdict["Partics"][partic]["ax"] = 0
				rectdict["Partics"][partic]["ay"] = 0

				#print("Done with force")
				#print(rectdict["Partics"][partic])
	#Når den har været igennem alle particles her, så går vi et level deeper
	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			#print("Going deeper")
			UpdateOld(rectdict[area])


def ClearDict(rectdict):#,points):
	#global points
	if "Partics" in rectdict:
		if len(rectdict["Partics"]) > 0:
			#for partic in rectdict["Partics"]:
			#	points[partic] = rectdict["Partics"][partic]
			points.update(rectdict["Partics"])
		
		rectdict["Partics"].clear()
	
	if "CM" in rectdict:
		rectdict.pop("CM")

	for area in ["A1","A2","A3","A4"]:
		if area in rectdict:
			ClearDict(rectdict[area])#,points)
			#rectdict[area].clear()
			rectdict.pop(area)

	#rectdict.clear()



G = 6.674*10**(-11)
AU = 1.5*10**11
Lx0 = -6*AU
Ly0 = -6*AU
Lx = 6*AU
Ly = 6*AU
dt = 300000

points = {}#"Sun":{},"Mercury":{},"Venus":{},"Earth":{}, "Jupiter":{}}
#Sun
points["P0"] = {"x": 0, 
					"y": 0, 
					"m": 2*10**30, 
					"u": 0, 
					"v": 0,
					"ax":0,
					"ay":0,
					"xold":0,
					"yold":0,
					"xnew":0,
					"ynew":0}
#Mercury
points["P1"] = {"x": 0, 
					"y": -46001200*10**3, 
					"m": 3.3*10**23, 
					"u": 47.362*10**3, 
					"v": 0,
					"ax":0,
					"ay":0,
					"xold":0,
					"yold":-46001200*10**3,
					"xnew":0,
					"ynew":0}

#Venus
points["P2"] = {"x": -108200000*10**3, 
					"y": 0, 
					"m": 4.87*10**24, 
					"u": 0, 
					"v": -35000,
					"ax":0,
					"ay":0,
					"xold":-108200000*10**3,
					"yold":0,
					"xnew":0,
					"ynew":0}
#earth
points["P3"] = {"x": 1.5*10**11, 
					"y": 0, 
					"m": 6*10**24, 
					"u": 0, 
					"v": 29800,
					"ax":0,
					"ay":0,
					"xold":1.5*10**11,
					"yold":0,
					"xnew":0,
					"ynew":0}

#Mars
points["P4"] = {"x": 0, 
					"y": 1.3814*1.5*10**11, 
					"m": 6.417*10**23, 
					"u": -24.077*10**3, 
					"v": 0,
					"ax":0,
					"ay":0,
					"xold":0,
					"yold":-46001200*10**3,
					"xnew":0,
					"ynew":0}				
#Jupiter
points["P5"] = {"x": 0, 
					"y": 778500000*10**3, 
					"m": 1.898*10**27, 
					"u": -13.07*10**3, 
					"v": 0,
					"ax":0,
					"ay":0,
					"xold":0,
					"yold":778500000*10**3,
					"xnew":0,
					"ynew":0}



maxlevel = 100

for pointkey in points:
	level = 0
	PointCoords(level,points[pointkey],pointkey,[Lx0,Lx,Ly0,Ly],Dict)


PrintDict()
print(Dict)
print(points.keys())

# print("Calculate Center of Mass of cells/areas")
# #CenterOfMass(Dict)

# #HERE DO STEP 1
level = 0
CenterOfMass(Dict)
CalcForceFirst(Dict,level)
UpdateOld(Dict)

print("Print Dict after clear")

points = {}
print(points.keys())
ClearDict(Dict)#,points)
print(Dict)
print(points.keys())

# #Dict.clear()

# PrintDict()

# #Sådan her kan jeg også clear hele dictionary btw!

# print("Print points after ClearDict")
# print(points)
# print("Print Len of points after ClearDict")
# print(len(points))

for pointkey in points:
	level = 0
	PointCoords(level,points[pointkey],pointkey,[Lx0,Lx,Ly0,Ly],Dict)

CenterOfMass(Dict)
#print(Dict)

#print("Calculate Forces")
#level = 0
#CalcForce(Dict,level)


# ###########################################################################
# #Printing and plotting
# print("NOW TRYING PLOT")
# level = 0
# Areas = ["A1","A2","A3","A4"]

#ListPlot = []


#PlotDict(level, Dict)

# #Lets print the dict
# print("PRINTING DICT")
# totalparticles = 0
# level = 0
# levuls = [0,0,0,0,0,0]

# levelarea = []
# PrintNumChildren(Dict,level,levelarea)
# print("Total number of counted particles")
# print(totalparticles)
# print("Total len of points")
# print(len(points))
# print("LET'S LOOK AT LEVELS")
# print(levuls)

# for point in points:
# 	plt.scatter(points[point]["x"],points[point]["y"], color = "blue",s=20)
# plt.axis([Lx0-1,Lx+1,Ly0-1,Ly+1])
# plt.plot([Lx0,Lx],[Ly0,Ly0], color = "black")
# plt.plot([Lx0,Lx],[Ly,Ly], color = "black")
# plt.plot([Lx0,Lx0],[Ly0,Ly], color = "black")
# plt.plot([Lx,Lx],[Ly0,Ly], color = "black")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Quadtree")
#plt.show()


#Anim here
#plt.hold(True)
fig = plt.figure(2)
ax = fig.gca()
ax.set_xlim(-6*AU,6*AU)
ax.set_ylim(-6*AU,6*AU)
ax.set_title('Kappa')
ax.set_ylabel('y [m]')
ax.set_xlabel('x [m]')


def animate(i): #i increment with 1 each step

	global ListPlot
	ax.clear()
	ax.set_xlim(-6*AU,6*AU)
	ax.set_ylim(-6*AU,6*AU)
	ax.set_title('Kappa')
	ax.set_ylabel('y [m]')
	ax.set_xlabel('x [m]')
	
	#PrintDict()
	# totalparticles = 0
	# level = 0
	# levuls = [0,0,0,0,0,0]
	# levelarea = []
	# PrintNumChildren(Dict,level,levelarea)
	# print("Total number of counted particles")
	# print(totalparticles)
	# print("Total len of points")
	# print(len(points))

	#PrintDict()
	
	level = 0
	CalcForce(Dict,level)

	CenterOfMass(Dict)
	
	level = 0
	ListPlot = []
	PlotDict(level, Dict)
	
	for pointkey in points:
		plt.scatter(points[pointkey]["x"],points[pointkey]["y"], color = "blue",s=20)



	UpdateOld(Dict)

	points.clear()
	ClearDict(Dict)
	#print(points)
	
	for pointkey in points:
		level = 0
		PointCoords(level,points[pointkey],pointkey,[Lx0,Lx,Ly0,Ly],Dict)
		

	return None


anim = animation.FuncAnimation(fig, animate, frames=10, interval=100)

plt.show()