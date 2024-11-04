import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

filename1 = 'SS_time.txt'
#filename2 = 'Times10.txt'
#filename2 = 'CO_Kurucz.dat'
#filename2 = '/groups/astro/kholten/CELS/MARCS/data/molecules/molecules_exomol/AlCl_exomol.dat'
#filename3 = '/groups/astro/kholten/Line_lists_uploaded/12C-31P__MoLLIST.trans' 
#filename = 'NO_HITEMP.dat'
#f1 = open(filename1)
#f2 = open(filename2)
#f3 = open(filename3)

#SpecFile1 = np.loadtxt(filename1, skiprows=19, dtype='float', max_rows=53)
SpecFile1 = np.loadtxt(filename1, dtype='float')
#SpecFile2 = np.loadtxt(filename2, dtype='float')
#SpecFile2 = np.loadtxt(filename2, skiprows=14, dtype='float')
#MARCSOS = np.loadtxt(filename2, skiprows=12, dtype='float')
#LLRune = np.loadtxt(filename2, dtype='float')
#LLKHM = np.loadtxt(filename3, dtype='float')

labels = ['k', 'T', 't_mix', 't_CO', 't_NH3', 't_CO2','t_HCN', 't_PH3']
#plt.plot(SpecFile1[:,0],SpecFile1[:,2]/SpecFile2[:,2])
#plt.show()

plt.plot(SpecFile1[:,1],SpecFile1[:,3])
#for i in range(6):
#	plt.plot(SpecFile1[:,1],SpecFile1[:,i+2],label=labels[i+2])
#	plt.plot(SpecFile2[:,1],SpecFile2[:,i+2],'.',label=labels[i+2])
plt.yscale('log')
#plt.plot(SpecFile2[:,0],SpecFile2[:,1],label='Regular')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,2],label='Incl O3')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,3],label='Excl O3')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,4],label='Excl O3,CH4')
#plt.plot(LLRune[:,3],LLRune[:,2],'+')
#plt.plot(LLKHM[:,3],LLKHM[:,2],'o')

#plt.xlim([10,10.3])
#plt.ylim([0,0.3])
#plt.yscale('log')
#plt.xlabel('Wavelength, microns')
plt.xlabel('Temperature, K')
plt.ylabel('Time for steady state, s')
#plt.legend(loc=1)
plt.show()

