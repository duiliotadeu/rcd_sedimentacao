# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import re
import matplotlib.animation as animation

# subrotinas definidas :
#    leitura(dirpath,fileidx) :  leitura do arquivo "SOLUTIONX_T???.out" com informacao da simulacao
#    lesimuinfo(dirpath)  : leitura de informacoes da simulacao
#    le_dadoexperimental() : 
#                           leitura do dado experimental que esta na pasta ,,,,,
#          Subrotina para um experimento especifico.
#    le_simucotas(dirpath) : 
#    make_video(dirpath, plota_experimento=False):
#    plot_cotas(dirpath)
#    plot_file()
#    video_diff()


# leitura: 
# subrotina para leitura com arquivo com o resultado da simulacao 
# que esta no diretorio  dirpath  
# e que tem indice para solucao:   fileidx
# pois cada resultado da simulacao esta indexada com os arquivos numerados
# chamada da funcao:   leitura(dirpath, fileidx)
#   dirpath :  string
#   fileidx : inteiro
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



def make_video(dirpath, plota_experimento=False):
    print('Creating Video from Folder : ' + dirpath)
    # passo 1. Leitura dos dados da Simulacao:
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

    phitot_inicial = 0.0
    for i in range(Ne):
        phitot_inicial = phitot_inicial+phi0[i]
    print('................................')
    print('Phi Tot Inicial: ', phitot_inicial)
    pos_deposicao = Domain_Length*phitot_inicial/phi_tot
    integral_concentration = Domain_Length*phitot_inicial
    print('pos. Deposicao: ', pos_deposicao)
    print('integral Concentracao Inicial: ', integral_concentration)

    if (plota_experimento):
        # Leitura dos dados do experimento
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
    # Set Axis Limits
    ax.axis([0,phi_tot+0.1,  -0.1, Domain_Length ])
    ax.set_ylabel('Position (m)')
    ax.set_xlabel('Concentration')
    
    ims={}
    list_colors = ['g','b', 'orangered', 'y', 'm', 'indigo', 'navy', 'cyan', 'crimson', 'tab:brown']
    for i in range(Ne):
        #mystr = "Part. "+str(i)+' '+str(diam[i]*1000)+'mm'
        mystr = str(diam[i]*1000)+' mm'
        im, = ax.plot([],[], list_colors[i], label=mystr)
        ims[i]=im

    imtot, = ax.plot([],[],'r', label='Concentracao Total')

    if (plota_experimento):
        imexp1, = ax.plot([],[],'o-m', label='exp.t0')
        imexp2, = ax.plot([],[],'+-k', label='exp.t1')

    ax.legend()
    ax.grid(True)

    def updatefig(i):
        data, simtime = leitura(dirpath, i)
        y = data[:,1]
        for k in range(1,Ne):
            y = y + data[:,k+1]
        imtot.set_data(y, data[:,0])
        
        for k in range(Ne):
            ims[k].set_data(data[:,k+1], data[:,0])

        
        if (plota_experimento):
            mystr = np.array2string(data_exp[idx_Texp[i], 0] ) + '  '+ np.array2string(data_exp[idx_Texp[i]+1, 0]) 
            imexp1.set_data(data_exp[idx_Texp[i], cota_exp, 1:6] )
            imexp2.set_data(data_exp[idx_Texp[i]+1, cota_exp,  1:6] )
            ax.set_title('Simulation Time : '+ "{:.2f}".format(simtime) + ' / Exp. Time: ' + mystr  ) 
        else:
            ax.set_title('Simulation Time : '+ "{:.2f}".format(simtime)+' s')

    ani = animation.FuncAnimation(fig, updatefig, frames= Num_Output,  interval=100, blit = False )
    video_filename = dirpath+'/video.mp4'
    ani.save(video_filename, dpi=300)
    plt.show()



def plot_cotas(dirpath):
    print('Reading Simulation Data...')
    data, xcota = le_simucotas(dirpath)

    # LEITURA DOS DADOS EXPERIMENTAIS
    print('Reading Experimental Data...')
    [data_exp, cota_exp, num_cotas_exp] = le_dadoexperimental()

    x = data[:,0]/60

    for i in range(5):  
        plt.figure()
        plt.xlim(0.0, 310)
        plt.ylim(0.0, 0.6)
        plt.plot(x, data[:,i+1], '+-b', label = 'Sim. at x='+ "{:.2f}".format(xcota[i]*100) +' cm')
        plt.plot(data_exp[:,0], data_exp[:,5-i]/100, 'o-r', label = 'experimental')
        plt.title('Cota x = '+ "{:.2f}".format(xcota[i]*100)+' cm')
        plt.xlabel('Time [min.]')
        plt.ylabel('Concentration')
        plt.legend()
        filename = dirpath+'/cota'+"{:.2f}".format(xcota[i]*100) +'.pdf'
        plt.savefig(filename, dpi=300)

    # # Plot 1 [COTAS]
    # fig, axs = plt.subplots(nrows=2,ncols=2)
    
    # fig.suptitle('Caso 3 particulas')
    # axs[0,0].plot(data[:,0]/60, data[:,1], '+-b', label = 'Sim. at x='+str(xcota[0]*100)+' cm')
    # axs[0,0].plot(data_exp[:,0], data_exp[:,5]/100, 'o-r', label = 'experimental')
    # axs[0,0].set_title('Cota x = '+ str(xcota[0]*100)+' cm')
    # axs[0,0].legend()

    # axs[0,1].plot(data[:,0]/60, data[:,2], '+-b', label = 'Sim. at x='+str(xcota[1]*100)+' cm')
    # axs[0,1].plot(data_exp[:,0], data_exp[:,4]/100, 'o-r', label = 'experimental')
    # axs[0,1].set_title('Cota x = '+ str(xcota[1]*100)+' cm')
    # axs[0,1].legend()

    # axs[1,0].plot(data[:,0]/60, data[:,3], '+-b', label = 'Sim. at x='+str(xcota[2]*100)+' cm')
    # axs[1,0].plot(data_exp[:,0], data_exp[:,3]/100, 'o-r', label = 'experimental')
    # axs[1,0].set_title('Cota x = '+ str(xcota[2]*100)+' cm')
    # axs[1,0].legend()

    # axs[1,1].plot(data[:,0]/60, data[:,4], '+-b', label = 'Sim. at x='+str(xcota[3]*100)+' cm')
    # axs[1,1].set_title('Cota x = '+ str(xcota[3]*100)+' cm')
    # axs[1,1].plot(data_exp[:,0], data_exp[:,2]/100, 'o-r', label = 'experimental')
    # axs[1,1].legend()

    # # Plot 2 :
    # plt.figure()
    # plt.plot(data[:,0]/60, data[:,5], '+-b', label = 'Sim. at x='+str(xcota[4]*100)+' cm')
    # #plt.set_title('Cota x = '+ str(xcota[4]*100)+' cm')
    # plt.plot(data_exp[:,0], data_exp[:,3]/100, 'o-r', label = 'experimental')
    # plt.legend()
    #plt.plot(x,np.sin(x))

    plt.show()


    
def plot_file(dirpath, fileidx):
    data,m = leitura(dirpath, fileidx)
    print(type(data))
    plt.plot(data[:,0], data[:,1],'r')
    plt.plot(data[:,0], data[:,2],'y')
    plt.plot(data[:,0], data[:,3],'g')
    plt.plot(data[:,0], data[:,4],'b')
    plt.plot(data[:,0], data[:,5],'r')
    plt.plot(data[:,0], data[:,6],'r')
    plt.plot(data[:,0], data[:,7],'g')
    
    #plt.plot(data[:,0], data[:,1]+data[:,2]+data[:,3],'b')
    #plt.plot(data[:,0], data[:,1]+data[:,2],'b')
    #plt.plot(data[:,0], data[:,1],'b')
    plt.show()


def plot_topo(dirpath):
    # LEITURA DO ARQUIVO
    filename = dirpath+'/topo.data'
    M=np.loadtxt(filename)
    plt.plot(M[:,0], M[:,1],'r')
    plt.show()



def video_diff():
    # Folders of experiments to compare in video 
    dirpath1 = 'glicerina1_3particulas'
    dirpath2 = 'glicerina5_3particulas'
    dirpath3 = 'glicerina6_3particulas'
    num_cases = 3

    # Data Info of 1st Simulation
    Nx, Domain_Length, T_Final, Ne, phi0, diam, rho_p, \
        phi_tot, hs_n, hs_lambda, mu_f, rho_f, Num_Output, Num_Cotas, cotas, SimuTimes = le_simuinfo(dirpath1)



    fig, ax=plt.subplots()
    # Set Axis Limits
    ax.axis([0,phi_tot+0.1,  -0.1, Domain_Length ])
    ax.set_ylabel('Position (m)')
    ax.set_xlabel('Concentration')
    
    ims={}
    list_colors = ['g','b', 'orangered', 'y', 'm', 'indigo', 'navy', 'cyan', 'crimson', 'tab:brown']
    for i in range(num_cases):
        mystr = 'Caso '+str(i+1)
        im, = ax.plot([],[], list_colors[i], label=mystr)
        ims[i]=im
    
    ax.legend()
    ax.grid(True)

    def updatefig(i):

        # simulation 1 
        data, simtime = leitura(dirpath1, i)
        y = data[:,1]
        for k in range(1,Ne):
            y = y + data[:,k+1]
        ims[0].set_data(y, data[:,0])

        # simulation 2
        data, simtime = leitura(dirpath2, i)
        y = data[:,1]
        for k in range(1,Ne):
            y = y + data[:,k+1]
        ims[1].set_data(y, data[:,0])

        # simulation 3
        data, simtime = leitura(dirpath3, i)
        y = data[:,1]
        for k in range(1,Ne):
            y = y + data[:,k+1]
        ims[2].set_data(y, data[:,0])

        ax.set_title('Simulation Time : '+ "{:.2f}".format(simtime)+' s')

    ani = animation.FuncAnimation(fig, updatefig, frames= Num_Output,  interval=100, blit = False )
    video_filename = 'video_comparacao.mp4'
    ani.save(video_filename, dpi=300)
    plt.show()




def compara_topo():
    cores = ['r', 'b', 'm', 'g']
    casos = ['glicerina3_1particula_25p', 'glicerina3_1particula_20p', 'glicerina3_1particula_15p', 'glicerina3_1particula_10p']
    lista_legendas=['25 %', '20 %', '15 %', '10 %']
    
    plt.rcParams.update({'font.size': 26})
    fig, ax = plt.subplots()

    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Altura (m)')
    ax.set_ylim(0, 1)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax.set_xticks([0, 2500, 5000,7500, 10000, 12500, 15000, 17500,20000,22500])

    for i in range(4): 
# LEITURA DO ARQUIVO
        filename = casos[i] + '/topo.data'
        print(filename)
        M=np.loadtxt(filename)
        ax.plot(M[:,0], M[:,1], cores[i], label=lista_legendas[i])
        #plt.show()

    legend = ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Particula do Grupo 2 : 0.079mm\n'+'Variando Concentração Inicial')

    plt.show()



def compara_topo2():
    cores = ['r', 'b', 'm', 'g']
    casos = ['glicerina2_2particula_25p', 'glicerina2_2particula_20p', 'glicerina2_2particula_15p', 'glicerina2_2particula_10p']
    lista_legendas=['25 %', '20 %', '15 %', '10 %']

    plt.rcParams.update({'font.size': 26})
    fig, ax = plt.subplots()

    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Altura (m)')
    ax.set_ylim(0, 1)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax.set_xticks([0, 2500, 5000,7500, 10000, 12500, 15000, 17500,20000,22500])

    for i in range(4): 
# LEITURA DO ARQUIVO
        filename = casos[i] + '/topo.data'
        print(filename)
        M=np.loadtxt(filename)
        ax.plot(M[:,0], M[:,1], cores[i], label=lista_legendas[i])
        #plt.show()

    legend = ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Grupo 2 |  Duas partículas: 0.066mm [50%]  0.092mm [50%]\n'  +'Variando Concentração Inicial')

    plt.show()


def compara_topo3():
    cores = ['r', 'b', 'm', 'g']
    casos = ['glicerina2_3particula_25p', 'glicerina2_3particula_20p', 'glicerina2_3particula_15p', 'glicerina2_3particula_10p']
    lista_legendas=['25 %', '20 %', '15 %', '10 %']

    plt.rcParams.update({'font.size': 26})
    fig, ax = plt.subplots()

    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Altura (m)')
    ax.set_ylim(0, 1)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    ax.set_xticks([0, 2500, 5000,7500, 10000, 12500, 15000, 17500,20000,22500])

    for i in range(4): 
# LEITURA DO ARQUIVO
        filename = casos[i] + '/topo.data'
        print(filename)
        M=np.loadtxt(filename)
        ax.plot(M[:,0], M[:,1], cores[i], label=lista_legendas[i])
        #plt.show()

    legend = ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Grupo 2 |  Três partículas: 0.053mm [25%]  0.079mm [50%] 0.105mm [25%]\n'  +'Variando Concentração Inicial')

    plt.show()

def compara_topo_particulas():
    cores = ['r', 'b', 'm', 'g']
    casos = ['glicerina2_1particula_25p', 'glicerina2_2particula_25p', 'glicerina2_3particula_25p']
    lista_legendas=['1 particula', '2 particulas', '3 particulas']

    plt.rcParams.update({'font.size': 26})
    fig, ax=plt.subplots()
    
    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Altura (m)')

    major_ticks = np.arange(0, 22500, 2500)
    minor_ticks = np.arange(0, 22500, 500)
    
    #ax.set_xticks([0, 2500, 3000, 3500, 4000, 4500, 5000,7500, 10000, 12500, 15000, 17500,20000,22500])
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)

    for i in range(3):
        filename = casos[i]+'/topo.data'
        print(filename)
        M=np.loadtxt(filename)
        ax.plot(M[:,0], M[:,1], cores[i], label=lista_legendas[i])
        
    legend = ax.legend(loc='upper right')
    ax.grid(which='both')
    ax.grid(which='minor', alpha=0.5)
    ax.grid(which='major', alpha=1.0)
    ax.set_title('Grupo2. Comparação com concentração inicial 25%')

    plt.show()

def compara_topo_particulas_zoom():
    cores = ['r', 'b', 'm', 'g']
    casos = ['glicerina2_1particula_25p', 'glicerina2_2particula_25p', 'glicerina2_3particula_25p']
    lista_legendas=['1 particula', '2 particulas', '3 particulas']

    plt.rcParams.update({'font.size': 26})
    fig, ax=plt.subplots()
    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Altura (m)')
    for i in range(3):
        filename = casos[i]+'/topo.data'
        
        print(filename)
        M=np.loadtxt(filename)
        #print(M.shape)
        
        ax.plot(M[100:,0], M[100:,1], cores[i], label=lista_legendas[i])

    legend = ax.legend(loc='upper right')
    ax.grid()
    ax.set_title('Grupo2. Comparação com concentração inicial 25%')

    plt.show()


# CHECK MASS CONSERVATION
