import numpy as np
import matplotlib.pyplot as plt
import glob
import os

files = glob.glob("output*.dat")
salida = open('salida.dat', 'w+')
for f in files:
    temp=os.path.basename(f).split('.')[0].split('_')[1]
    if temp!='0':
       nombre=temp
n_files = len(files)
def energy(data):
    E = data[:,6].sum() + data[:,7].sum() 
    return E

def radius(data):
    E_pot = data[:,6]
    min_pot = np.argmin(E_pot)
    #print "min_pot", min_pot
    x = data[:,0] - data[min_pot, 0]
    y = data[:,1] - data[min_pot, 1]
    z = data[:,2] - data[min_pot, 2]
    r = np.sqrt(x**2 + y**2 +z**2)
    r = np.sort(r)
    return r[1:]



i_snap = 0
data_init = np.loadtxt("output_{}.dat".format(i_snap))
E_init = energy(data_init)
r_init = radius(data_init)


# time ./a.out 450 0.1
i_snap = nombre
data_init = np.loadtxt("output_{}.dat".format(i_snap))
E_final = energy(data_init)
r_final = radius(data_init)
r_final = np.sort(r_final)
log_r_final = np.log10(r_final)

h, c = np.histogram(log_r_final,bins=30)
porcentaje=int(round((E_final-E_init)*100/E_final))
print('El porcentaje de energia conservada es de'+' '+str(porcentaje)+'%')


log_r_center = 0.5 * (c[1:]+c[:-1])
for i in range(len(h)):
    if h[i]!=0:
       salida.write(str(log_r_center[i])+' '+str(np.log10(h[i])-2.0*log_r_center[i])+'\n')
salida.close()
