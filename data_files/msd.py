import numpy as np

def minImageDist (pos1, pos2, box):
	d = 0
	for i in range(0, 3):
		dx = pos2[i]-pos1[i]
		while (dx < -box[i]/2.0):
			dx += box[i]
		while (dx > box[i]/2.0):
			dx -= box[i]
		d += dx*dx
	return d


file = open("trajectory.xyz", 'r')
pos = []

nFrames = 100
nAtoms = 4000
L = 16.796
for i in range(0, nFrames):
	file.readline()
	file.readline()
	p2 = []
	for j in range(0, nAtoms):
		data = file.readline().strip().split()
		p2.append([float(data[1]), float(data[2]), float(data[3])])
	pos.append(p2)

file.close()
	
msd	= []
nMsd = []
for t1 in range(0, nFrames):
	msd.append(0.0)
	nMsd.append(0)
	
box = [L, L, L]
for t1 in range(0, nFrames):
	print "t1 = "+str(t1)+"\n"
	for t2 in range(t1+1, nFrames):
		dt = t2-t1
		sqDis = 0
		for atom in range(0, nAtoms):
			sqDis += minImageDist(pos[t2][atom], pos[t1][atom], box)
		msd[dt] += sqDis/nAtoms
		nMsd[dt] += 1
		
for t1 in range(0, nFrames):
	if (nMsd[t1] > 0):
		msd[t1] /= nMsd[t1] 
		
file = open('msd.dat', 'w')
for t1 in range(0, nFrames):
	file.write(str(t1)+"\t"+str(msd[t1])+"\n")
file.close()