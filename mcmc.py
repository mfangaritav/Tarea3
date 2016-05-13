from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from math import *
import scipy as sp

data=np.loadtxt("salida.dat")
x=data[:,0]
y=data[:,1]
posel=0
termine=False
for i in range(len(x)):
    if x[i]>1.5 and not termine:
       posel=i
       termine=True
x=data[:posel,0]
y=data[:posel,1]

x1=np.linspace(-1.5,2,100)

def poisson(k, lamb):
    a=(lamb**k/sp.special.gamma(k+1))
    return a * np.exp(-lamb)

def likelihood(x, y, sigmay, a, b,c,d):
    funcion=np.log10((10**a)/(((10**x/10**b)**c)*((1+10**x/10**b)**d)))
    return (np.sum(np.log(1./(np.sqrt(2.*np.pi) * sigmay))) +
            np.sum(-0.5 * (y - funcion)**2 / sigmay**2))

def encontrarpxrc(x,rc):
    posicion=0;
    for i in range(len(x)):
        if x[i]>rc:
           posicion=i-1
    return posicion
def daroptimo(p,a,b,c,d,final,n):
    pini=[a,b,c,d]
    x1=np.linspace(-1.5,2.5,100)
    x2=np.linspace(pini[p],final,n)
    minimo=likelihood(x, y, sigmay, pini[0], pini[1],pini[2],pini[3])
    pminimo=pini[p]
    plt.figure(9)
    for i in range(n):
#        print pini[0], pini[1],pini[2],pini[3]
        pini[p]=pini[p]+(final-pini[p])/n
        if p==1:
           temp=likelihood(x, y, sigmay, pini[0], pini[1],pini[2],pini[3])
           y3=np.log10(10**pini[0]/(((10**x1/10**(pini[1]))**(pini[2]))*((1+10**x1/10**(pini[1]))**(pini[3]))))
           plt.plot(x1,y3)
           plt.scatter(x,y)
        else:
           j=encontrarpxrc(x,pini[1])
           if p==0:              
              temp=likelihood(x[:j], y[:j], sigmay, pini[0], pini[1],pini[2],pini[3])
           if p==2:              
              temp=likelihood(x[:j], y[:j], sigmay, pini[0], pini[1],pini[2],pini[3])
           if p==3:              
              temp=likelihood(x[j:], y[j:], sigmay, pini[0], pini[1],pini[2],pini[3])
        if temp>minimo:
           pminimo=pini[p]
#        print pminimo
    return pminimo
a,b,c,d = 1.0, 0.0,0.0,3.0

astep, bstep,cstep, dstep = 10.0,1.0,0.1,0.5
        
nsteps = 5000
aes=[]
bs=[]
cs=[]
ds=[]
minimo=-100000
pminimo=-1
#plt.figure(5)
for j in range(50):
 a,b,c,d = 1.0, 1.0,1.0,1.0   
 chain = []
 probs = []
 naccept = 0
 sigmay=1

 L_old    = likelihood(x, y, sigmay, a, b,c,d)
 
 prob_old = np.exp(L_old)

 for i in range(nsteps):
     anew = a + np.random.normal() * astep
     bnew = b + np.random.normal() * bstep
     cnew = c + np.random.normal() * cstep
     dnew = d + np.random.normal() * dstep

     L_new    = likelihood(x, y, sigmay, anew, bnew,cnew, dnew)
     prob_new = np.exp(L_new)

     if (prob_new / prob_old > np.random.uniform()):
         a = anew
         b = bnew
         c = cnew
         d = dnew
         L_old = L_new
         prob_old = prob_new
         naccept += 1
     else:
         pass

     chain.append((a,b,c,d))
     probs.append((L_old))
 aa = [a for a,b,c,d in chain]
 aes.append(aa)
 bb = [b for a,b,c,d in chain]
 bs.append(bb)
 cc = [c for a,b,c,d in chain]
 cs.append(cc)
 dd = [d for a,b,c,d in chain]
 ds.append(dd)

 for k in range(15):
     na,binsa=histogram(aa,5+k,normed=True)
     (mua,sigmaa)=norm.fit(aa)

     nb,binsb=histogram(bb,5+k,normed=True)
     (mub,sigmab)=norm.fit(bb)
 
     nc,binsc=histogram(cc,5+k,normed=True)
     (muc,sigmac)=norm.fit(cc)
 
     nd,binsd=histogram(dd,5+k,normed=True)
     (mud,sigmad)=norm.fit(dd)
 y3=np.log10(10**mua/(((10**x1/10**(mub))**(muc))*((1+10**x1/10**(mub))**(mud)))) 
# plt.plot(x1,y3)
# plt.scatter(x,y)
 likeli=likelihood(x, y, sigmay, mua, mub,muc,mud)
 if likeli>minimo:
    minimo=likeli
    pminimo=j
# print minimo,pminimo
plt.figure(1)
na,binsa,p=plt.hist(aes[pminimo],20,normed=True)
(mua,sigmaa)=norm.fit(aes[pminimo])
ya=norm.pdf(binsa,mua,sigmaa)
plt.title('Histograma del logaritmo del parametro Rho_0')
plt.ylabel('Numero de log(Rho_s)\'s')
plt.xlabel('log(Rho_0)')
plt.plot(binsa,ya,'r--',color='black')

plt.figure(2)
nb,binsb,p=plt.hist(bs[pminimo],20,normed=True,color='green')
(mub,sigmab)=norm.fit(bs[pminimo])
yb=norm.pdf(binsb,mub,sigmab)
plt.title('Histograma del logaritmo del parametro rc')
plt.ylabel('Numero de log(rc)\'s')
plt.xlabel('log(rc)')
plt.plot(binsb,yb,'r--',color='black')

plt.figure(3) 
nc,binsc,p=plt.hist(cs[pminimo],20,normed=True,color='gray')
(muc,sigmac)=norm.fit(cs[pminimo])
yc=norm.pdf(binsc,muc,sigmac)
plt.title('Histograma del parametro alpha')
plt.ylabel('Numero de alpha\'s')
plt.xlabel('alpha')
plt.plot(binsc,yc,'r--',color='black')

plt.figure(4)
nd,binsd,p=plt.hist(ds[pminimo],20,normed=True,color='red')
(mud,sigmad)=norm.fit(ds[pminimo])
yd=norm.pdf(binsd,mud,sigmad)
plt.title('Histograma del parametro beta')
plt.ylabel('Numero de beta\'s')
plt.xlabel('beta')
plt.plot(binsd,yd,'r--',color='black')

plt.figure(0)
plt.scatter(x,y)
plt.title('Grafica de los datos y el modelo de MCMC')
plt.ylabel('log(Rho)')
plt.xlabel('log(r)')
y1=np.log10(10**mua/(((10**x1/10**(mub))**(muc))*((1+10**x1/10**(mub))**(mud))))

print 'Rho_0 es '+str(10**mua)+' con una incertidumbre de +- '+str(10**sigmaa)
print 'rc es '+str(10**mub)+' con una incertidumbre de +- '+str(10**sigmab)
print 'alpha es '+str(muc)+' con una incertidumbre de +- '+str(sigmac)
print 'beta es '+str(mud)+' con una incertidumbre de +- '+str(sigmad)
'''
mua1=binsa[np.argmax(na)]+(binsa[np.argmax(na)+1]-binsa[np.argmax(na)])/2
mub1=binsb[np.argmax(nb)]+(binsb[np.argmax(nb)+1]-binsb[np.argmax(nb)])/2
muc1=binsc[np.argmax(nc)]+(binsc[np.argmax(nc)+1]-binsc[np.argmax(nc)])/2
mud1=binsd[np.argmax(nd)]+(binsd[np.argmax(nd)+1]-binsd[np.argmax(nd)])/2

muc1=daroptimo(2,mua1i,mub1i,muc1i,mud1i,binsa[np.argmax(nc)+1],100)
mud1=daroptimo(3,mua1i,mub1i,muc1i,mud1i,binsa[np.argmax(nd)+1],100)
mua1=daroptimo(0,mua1i,mub1i,muc1i,mud1i,binsa[np.argmax(na)+1],100)
mub1=daroptimo(1,mua1i,mub1i,muc1i,mud1i,binsa[np.argmax(nb)+1],100)

mub1=daroptimo(1,mua1i,mub1i,muc1i,mud1i,binsb[np.argmax(nb)+1],100)
mua1=daroptimo(0,mua1i,mub1,muc1i,mud1i,binsa[np.argmax(na)+1],100)
muc1=daroptimo(2,mua1,mub1,muc1i,mud1i,binsc[np.argmax(nc)+1],100)
mud1=daroptimo(3,mua1,mub1,muc1,mud1i,binsd[np.argmax(nd)+1],100)
'''
#y2=np.log10(10**2.1/(((10**x1/10**(-0.28))**(0.0))*((1+10**x1/10**(-0.28))**(4.8/1.7))))

plt.plot(x1,y1)

#plt.plot(x1,y2)
plt.show()
