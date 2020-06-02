import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob

def plot(path, file):
	filename="%s/%s"%(path,file)
	if not (os.path.isfile(filename)):
		print("File doesn't exist") 
		return

	print("Plotting file %s in %s"%(file,path))

	colors = ['ro','ko','bo']

	pars = open('%s/parameters.dat'%(path),'r')
	tmp = pars.readline()
	tmp = pars.readline()
	L = int(tmp.strip().split(':')[1])
	pars.close()

	cwd = os.getcwd()
	os.chdir(path)

	t_array = []

	plt.figure(figsize = (10,10))

	f_in = open(file,'r')
	while True:
		tmp = f_in.readline()
		if tmp and len(tmp.strip().split()) > 1:
			t = float(tmp.strip().split()[0])
			t_array.append(t)
			nclust = int((len(tmp.strip().split())-1)/3)
			for i in range(nclust):
				CoM = float(tmp.strip().split()[i*3+1])
				vel = int(tmp.strip().split()[i*3+2])
				size = int(tmp.strip().split()[i*3+3])
				plt.plot(CoM,t,colors[vel+1], ms = 0.1*size)
		elif not tmp:
			break

	tmax = max(t_array)
	tmin = min(t_array)
	plt.xlim(0,L)
	plt.ylim(tmin,tmax)
	plt.tight_layout()
	plt.savefig('%s.png'%(file.split('.dat')[0]))

	os.chdir(cwd)

def main():

	seed = sys.argv[1]
	file = sys.argv[2]

	dirs = glob.glob('../../output_1DsppBD/%s'%(seed))

	for path in dirs:
		plot(path, file)

if __name__ == '__main__':
	main()