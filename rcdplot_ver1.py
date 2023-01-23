# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib.animation as animation

def leitura(dirpath, fileidx):
	filenumber = str(int(fileidx))
	filename = dirpath+'/SOLUTIONX_T'+filenumber+'.out'
	
	#print(filename)
	fp = open(filename, "r")
	line_time = fp.readline()
	mystr = line_time.split(" ")
	simtime = float(mystr[1])
	line_nx = fp.readline()
	line_ne = fp.readline()
	
	#print(line_time)
	#print(line_ne)
	fp.close()
	
	M=np.loadtxt(filename, skiprows=3)
	
	return M, simtime


def le_simuinfo(dirpath):
	filename = dirpath+'/simulation.info'
	fp = open(filename,"r")
	fp.readline()
	# Read Number of points in x
	line = fp.readline()
	mystr = line.split(":")
	Nx = int(mystr[1])
	# print(type(Nx), Nx)
	
	# Read Domain Length
	line = fp.readline()
	mystr = line.split(":")
	Domain_Length = float(mystr[1])

	# Read Simulation Total Time
	line = fp.readline()
	mystr = line.split(":")
	T_Final = float(mystr[1])

	# Read Delta T
	line = fp.readline()
	mystr = line.split(":")
	Delta_T = float(mystr[1])

	# Read Number of Particles
	line = fp.readline()
	line = fp.readline()
	mystr = line.split(":")
	Ne = int(mystr[1])

	# Read Particle Initial concentrations / Diameter / Density
	line = fp.readline()
	phi0=np.zeros(Ne)
	diam=np.zeros(Ne)
	rho_p=np.zeros(Ne)
	for i in range(Ne):
		line = fp.readline()
		mystr = line.split("  ")
		phi0[i] = float(mystr[0])
		diam[i] = float(mystr[1])
		rho_p[i] = float(mystr[2])
		

	# Read Maximal Total Concentration
	line = fp.readline()
	mystr = line.split(":")
	phi_tot = float(mystr[1])

	# Read Hindered Settling Function Data
	line = fp.readline()
	mystr = line.split(":")
	hs_n = float(mystr[1])
	line = fp.readline()
	mystr = line.split(":")
	hs_lambda = float(mystr[1])

	# Read Fluid Viscosity
	line = fp.readline()
	mystr = line.split(":")
	mu_f = float(mystr[1])

	# Read Fluid Density
	line = fp.readline()
	mystr = line.split(":")
	rho_f = float(mystr[1])

	# Read Number of Output Files
	line=fp.readline()
	line = fp.readline()
	mystr = line.split(":")
	Num_Output = int(mystr[1])

	# Read Cotas 
	line = fp.readline()
	mystr = line.split(":")
	Num_Cotas = int(mystr[1])
	line = fp.readline()
	cotas = np.zeros(Num_Cotas)
	for i in range(Num_Cotas):
		line = fp.readline()
		mystr = line.split(":")
		cotas[i] = float(mystr[1])
	for i in range(7):
		line = fp.readline()
		# print(line)
	
	# Read times
	SimuTimes = np.zeros(Num_Output)
	print("Reading Simulation Times")
	for i in range(Num_Output):
		line = fp.readline()
		# print(i, "Line: ", line)
		SimuTimes[i] = float(line)

	# Return Data...
	return Nx, Domain_Length, T_Final, Ne, phi0, diam, rho_p, \
		phi_tot, hs_n, hs_lambda, mu_f, rho_f, Num_Output, Num_Cotas, cotas, SimuTimes

def le_dadoexperimental():
	filename = 'dadosexperimentais/dados_20porcento.txt'
	fp = open(filename, 'r')
	mystr = fp.readline()

	num_cotas=5
	cota = np.zeros(5)
	cota[0] = 87.25
	cota[1] = 67.25
	cota[2] = 27.25
	cota[3] = 17.25
	cota[4] = 7.25

	fp.close()
	M = np.genfromtxt(filename, skip_header=1)

	# num_tempos = 61
	# M = np.zeros((num_tempos,num_cotas+1))
	# for k in range(num_tempos):
	# 	mystr = fp.readline()
	# 	# remove duplicated spaces, tab, etc.
	# 	mystr = " ".join(mystr.split())  
	# 	values = mystr.split(" ")
	# 	for i in range(len(values)):
	# 		values[i] = values[i].replace(",",".")
	# 		M[k,i] = float(values[i]) 
	# 	#print(values)
	return M, cota, num_cotas

def le_simucotas(dirpath):
	filename = dirpath + "/cotas.data"
	fp = open(filename, 'r')
	#fp.readline()

	# FIRST READ THE POSITIONS IN X OF COTA VALUES
	line = fp.readline()
	mystr = line.split(" ")
	n = len(mystr)-1
	x_cota = np.zeros(n)
	for i in range(n):
		x_cota[i] = float(mystr[i])

	#line = fp.readline()

	# READ COTAS
	#M = np.zeros( (1, n+1) )
	#line = fp.readline()
	
	#data = numpy.genfromtxt(yourFileName,skiprows=n)

	fp.close()
	data = np.genfromtxt(filename, skip_header=2)


	# for k in range(num_tempos):
	# 	line = fp.readline()
	# 	values = line.split(" ")
	# 	for i in range (n+1):
	# 		M[k,i] = float(values[i])
	
	return data, x_cota



def make_video(dirpath):

	plota_experimento = False

	print('Creating Video from Folder : ' + dirpath)
	print('Reading Simulation Info...')
	Nx, Domain_Length, T_Final, Ne, phi0, diam, rho_p, \
		phi_tot, hs_n, hs_lambda, mu_f, rho_f, Num_Output, Num_Cotas, cotas, SimuTimes = le_simuinfo(dirpath)
	print('Simulation Info:  ')
	print('Nx : ', Nx)
	print('Domain Length : ', Domain_Length)
	print('Simulation Final Time : ', T_Final)
	print('Number of Particle      : ', Ne)
	print('Initial Concentrations  : ', phi0)
	print('Diameter of Particles   : ', diam)
	print('Density  of Particles   : ', rho_p)
	print('Maximum Total Concentration : ', phi_tot)
	print('Hindered Settling Func. - Lambda : ', hs_lambda)
	print('Hindered Settling Func. - n      : ', hs_n)
	print('Fluid Viscosity : ', mu_f)
	print('Fluid Density   : ', rho_f)
	print('Number of File Outputs : ', Num_Output)

	if (plota_experimento):
		print('Reading Experimental Data...')
		[data_exp, cota_exp, num_cotas_exp] = le_dadoexperimental()
		data_exp[:,0]=data_exp[:,0]*60  # convert to seconds
		data_exp[:,1:6]=data_exp[:,1:6]/100 # convert to [0,1]
		cota_exp = cota_exp/100 # convert cm to meters

		# Create association among times: simulated and experimental
		nt_exp  = data_exp.shape[0]
		nt_simu = Num_Output 
		idx_Texp = np.zeros(Num_Output, dtype=np.int64)
		for i in range(Num_Output):
			found = False
			k = 0
			while (not found):
				if ((data_exp[k,0] <= SimuTimes[i]) and ( data_exp[k+1,0] >= SimuTimes[i] )):
					idx_Texp[i] = k
					found = True
				k=k+1

	
	fig, ax=plt.subplots()
	ax.axis([-0.1,3.0, 0, 0.6 ])
	ax.set_xlabel('Position (m)')
	ax.set_ylabel('Concentration')
	
	im1, = ax.plot([],[],'g', label='p. 1')
	im2, = ax.plot([],[],'b', label='p. 2')
	im3, = ax.plot([],[],'y', label='p. 3')
	im4, = ax.plot([],[],'navy', label='p. 4')
	im5, = ax.plot([],[],'m', label='p. 5')
	im6, = ax.plot([],[],'indigo', label='p. 6')
	im7, = ax.plot([],[],'orangered', label='p. 7')
	im8, = ax.plot([],[],'cyan', label='p. 8')
	im9, = ax.plot([],[],'crimson', label='p. 9')
	im10, = ax.plot([],[],'turquoise', label='p. 10')


	imtot, = ax.plot([],[],'r', label='Total Concentration')

	if (plota_experimento):
		imexp1, = ax.plot([],[],'o-m', label='exp0')
		imexp2, = ax.plot([],[],'+-k', label='exp1')

	ax.legend()
	ax.grid(True)

	idxexp=0
	value2=0
	def updatefig(i):
		data, simtime = leitura(dirpath, i)
		imtot.set_data(data[:,0], data[:,1]+data[:,2]+data[:,3] + data[:,4]+data[:,5]+data[:,6] + \
			data[:,7]+data[:,8]+data[:,9]+data[:,10])
		
		im1.set_data(data[:,0], data[:,1])
		im2.set_data(data[:,0], data[:,2])
		im3.set_data(data[:,0], data[:,3])
		im4.set_data(data[:,0], data[:,4])
		im5.set_data(data[:,0], data[:,5])
		im6.set_data(data[:,0], data[:,6])
		im7.set_data(data[:,0], data[:,7])
		im8.set_data(data[:,0], data[:,8])
		im9.set_data(data[:,0], data[:,9])
		im10.set_data(data[:,0], data[:,10])
		
		if (plota_experimento):
			mystr = np.array2string(data_exp[idx_Texp[i], 0] ) + '  '+ np.array2string(data_exp[idx_Texp[i]+1, 0]) 
			imexp1.set_data(cota_exp, data_exp[idx_Texp[i], 1:6] )
			imexp2.set_data(cota_exp, data_exp[idx_Texp[i]+1, 1:6] )
			ax.set_title('Simu. Time : '+ "{:.2f}".format(simtime) + ' / Exp. Time: ' + mystr  ) 
		else:
			ax.set_title('Simu. Time : '+ "{:.2f}".format(simtime) )

	ani = animation.FuncAnimation(fig, updatefig, frames= Num_Output,  interval=100, blit = False )
	video_filename = dirpath+'/video.mp4'
	ani.save(video_filename)
	
	
	plt.show()



def plot_cotas(dirpath):
	print('Reading Simulation Data...')
	data, xcota = le_simucotas(dirpath)

	# LEITURA DOS DADOS EXPERIMENTAIS
	print('Reading Experimental Data...')
	[data_exp, cota_exp, num_cotas_exp] = le_dadoexperimental()

	# Plot 1 [COTAS]
	fig, axs = plt.subplots(nrows=2,ncols=2)
	fig.suptitle('Caso 3 particulas')
	axs[0,0].plot(data[:,0]/60, data[:,1], '+-b', label = 'Sim. at x='+str(xcota[0]*100)+' cm')
	axs[0,0].plot(data_exp[:,0], data_exp[:,5]/100, 'o-r', label = 'experimental')
	axs[0,0].set_title('Cota x = '+ str(xcota[0]*100)+' cm')
	axs[0,0].legend()

	axs[0,1].plot(data[:,0]/60, data[:,2], '+-b', label = 'Sim. at x='+str(xcota[1]*100)+' cm')
	axs[0,1].plot(data_exp[:,0], data_exp[:,4]/100, 'o-r', label = 'experimental')
	axs[0,1].set_title('Cota x = '+ str(xcota[1]*100)+' cm')
	axs[0,1].legend()

	axs[1,0].plot(data[:,0]/60, data[:,3], '+-b', label = 'Sim. at x='+str(xcota[2]*100)+' cm')
	axs[1,0].plot(data_exp[:,0], data_exp[:,3]/100, 'o-r', label = 'experimental')
	axs[1,0].set_title('Cota x = '+ str(xcota[2]*100)+' cm')
	axs[1,0].legend()

	axs[1,1].plot(data[:,0]/60, data[:,4], '+-b', label = 'Sim. at x='+str(xcota[3]*100)+' cm')
	axs[1,1].set_title('Cota x = '+ str(xcota[3]*100)+' cm')
	axs[1,1].plot(data_exp[:,0], data_exp[:,2]/100, 'o-r', label = 'experimental')
	axs[1,1].legend()

	# Plot 2 : Plot at give
	plt.figure()
	plt.plot(data[:,0]/60, data[:,5], '+-b', label = 'Sim. at x='+str(xcota[4]*100)+' cm')
	#plt.set_title('Cota x = '+ str(xcota[4]*100)+' cm')
	plt.plot(data_exp[:,0], data_exp[:,3]/100, 'o-r', label = 'experimental')
	plt.legend()
	#plt.plot(x,np.sin(x))

	plt.show()


	
def plot_file(dirpath, fileidx):
	data = leitura(dirpath, fileidx)
	plt.plot(data[:,0], data[:,1],'r')
	plt.plot(data[:,0], data[:,2],'y')
	plt.plot(data[:,0], data[:,3],'g')
	plt.plot(data[:,0], data[:,1]+data[:,2]+data[:,3],'b')
	plt.show()

