import numpy as np 
import matplotlib.pyplot as plt 

a = np.loadtxt('results\\one_period.txt')
#a_old = np.loadtxt('one_period_old.txt')
b = np.loadtxt('results\\poincare_map.txt').transpose()

#print(a[1,1]-a[-1,1])



# Indexing: 
    #zacznij od 400 to a[400:]
    #ostatnie 400 to a[:400]

fig, ax = plt.subplots() #Tworze obiekt fig with axes ax 
#ax.plot(a[:400,0], a[:400,1]) # metoda do axes na rysowanie
ax.plot(a[:,0], a[:,1])
plt.show()

fig, ax = plt.subplots()
#ax.plot(a[:400,1], a[:400,2]) 
ax.plot(a[:,1], a[:,2])
#ax.plot(a_old[:,1], a_old[:,2])
ax.scatter(b[10:,0], b[10:,1], c='red')
plt.show()
