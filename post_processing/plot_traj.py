import numpy as np 
import matplotlib.pyplot as plt 

a = np.loadtxt('trajectory.txt')
#aa = np.loadtxt('trajectory_old.txt')
#a = np.round(a,5)
b = np.loadtxt('poincare_map.txt').transpose()
#b = np.transpose(b)


ofset = 0 
#print(a[1,1]-a[-1,1])

fig, ax = plt.subplots()
#ax.scatter(a[bb:,0], a[bb:,1], s = 1, c = 'blue')
ax.plot(a[ofset:,0], a[ofset:,1],c = 'blue')
#ax.plot(aa[bb:,0], aa[bb:,1],linestyle='dashed',c = 'green')
ax.scatter(b[:,-1], b[:,0], c='red')
plt.show()

fig, ax = plt.subplots()
ax.plot(a[ofset:,1], a[ofset:,2])
#ax.scatter(b[:,0], b[:,1], c='red')
plt.show()




bb = int(0.1*len(a[:,1]))
bb = len(a[:,1])-1*400


#fig, ax = plt.subplots() #Tworze obiekt fig with axes ax 
#ax.scatter(a[bb:,0], a[bb:,2], s=1) # metoda do axes na rysowanie
#ax.scatter(b[:,2], b[:,1], s=5)

#for i in range(len(a[bb:,0])):
#    print(i)
#    ax.annotate(str(i), (a[(-i),0], a[(-i),2]))
#plt.show()


#fig, ax = plt.subplots() 
#ax.plot(a[bb:,1], a[bb:,2]) 
#ax.scatter(b[:,0], b[:,1])
#plt.show()

#fig, ax = plt.subplots()  
#ax.plot(a[:,0], a[:,2]) 
#ax.plot(a[:,0], a[:,4]) 




#print(b)
#fig, ax = plt.subplots() 
#ax.scatter(b[:,0], b[:,1]) 
#plt.show()


