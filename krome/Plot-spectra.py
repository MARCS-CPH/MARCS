import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

filename1 = 'build/fort.66'
#filename2 = 'Times10.txt'
#filename2 = 'CO_Kurucz.dat'
#filename2 = '/groups/astro/kholten/CELS/MARCS/data/molecules/molecules_exomol/AlCl_exomol.dat'
#filename3 = '/groups/astro/kholten/Line_lists_uploaded/12C-31P__MoLLIST.trans' 
#filename = 'NO_HITEMP.dat'
#f1 = open(filename1)
#f2 = open(filename2)
#f3 = open(filename3)

#SpecFile1 = np.loadtxt(filename1, skiprows=19, dtype='float', max_rows=53)
#Header = np.loadtxt(filename1, dtype='str',max_rows=1,comments='None')
#SpecFile1 = np.loadtxt(filename1, dtype='float', skiprows=1)
SpecFile1 = np.loadtxt(filename1, dtype='float')
#SpecFile2 = np.loadtxt(filename2, skiprows=14, dtype='float')
#MARCSOS = np.loadtxt(filename2, skiprows=12, dtype='float')
#LLRune = np.loadtxt(filename2, dtype='float')
#LLKHM = np.loadtxt(filename3, dtype='float')
#print(Header)
#labels = ['k', 'T', 't_mix', 't_CO', 't_NH3', 't_CO2','t_HCN', 't_PH3']
#plt.plot(SpecFile1[:,0],SpecFile1[:,2]/SpecFile2[:,2])
#plt.show()


#labels = ['O2','O','O3']
labels = ['N2', 'O', 'NO', 'N', 'O2']
#colors = ['b', 'o', 'g']
for i in range(len(SpecFile1[1,:])-1):
	plt.plot(SpecFile1[:,0],SpecFile1[:,i+1],label=labels[i])
#plt.plot(SpecFile1[:,0],SpecFile1[:,1],label='O2')
#plt.plot(SpecFile1[:,0],SpecFile1[:,3],label='O3')
#plt.plot(SpecFile1[:,0],SpecFile1[:,2],label='O')

#for i in range(6):
#	plt.plot(SpecFile1[:,1],SpecFile1[:,i+2],label=labels[i+2])
#	plt.plot(SpecFile2[:,1],SpecFile2[:,i+2],'.',label=labels[i+2])
#plt.yscale('log')
#plt.plot(SpecFile2[:,0],SpecFile2[:,1],label='Regular')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,2],label='Incl O3')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,3],label='Excl O3')
#plt.plot(1E4/SpecFile1[:,0],SpecFile1[:,4],label='Excl O3,CH4')
#plt.plot(LLRune[:,3],LLRune[:,2],'+')
#plt.plot(LLKHM[:,3],LLKHM[:,2],'o')

#plt.xlim([10,10.3])
#plt.ylim([0,0.3])
plt.yscale('log')
plt.xscale('log')
#plt.xlabel('Wavelength, microns')
plt.xlabel('Time, yr')
plt.ylabel('Mixing ratio')
#plt.ylim([1E-50,1])
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Number density [molecules cm$^{-3}$]')
plt.title('KROME solver')
#plt.xlim([1E-28,1E7])
plt.show()

print(SpecFile1[-1,1])
print(SpecFile1[-1,2])
print(SpecFile1[-1,3])
