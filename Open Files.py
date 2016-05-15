import matplotlib.pyplot as plt
import numpy as np

CalcX = True
Tc = 0.5 * np.log(1 + np.sqrt(2))
avgMags = []
avgMagsSquare = []
avgMagsNorm = []
susceptibility = [[],[],[],[],[]]

fList = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]
lList = [5, 10, 20, 30, 50]

for i in lList:
    avgMags.append(open('avgMag for L = %s.txt' % i, 'r').read().split())
    avgMagsNorm.append(open('avgMag for L = %s.txt' % i, 'r').read().split())

for i in range(len(lList)):
    for j in range(len(fList)):
        avgMagsNorm[i][j] = float(avgMagsNorm[i][j]) / lList[i]**2

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for i in range(len(lList)):
    ax1.plot(fList, avgMagsNorm[i], lw = '2', label = 'L= %s' % lList[i])
    
plt.legend()
plt.xlabel('T/Tc')
plt.ylabel('Magnetization per Lattice Site')
plt.title('Magnetization for LxL Lattices at Different Temperatures')

if CalcX:
    for i in lList:
        avgMagsSquare.append(open('avgMagSquare for L = %s.txt' % i, 'r').read().split())
        
    for i in range(len(lList)):
        for j in range(len(fList)):        
            temp = fList[i] * Tc
            susceptibility[i].append(1 / (temp * fList[j]**2) * 
                        (float(avgMagsSquare[i][j]) - float(avgMags[i][j])**2))
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)       
    for i in range(len(lList)):
        ax2.plot(fList, susceptibility[i], lw = '2', label = 'L = %s' % lList[i])
    
    plt.legend()
    plt.xlabel('T/Tc')
    plt.ylabel('Susceptibility')
    plt.title('Susceptibility for LxL Lattices at Different Temperatures') 
    
plt.plot()
        
    