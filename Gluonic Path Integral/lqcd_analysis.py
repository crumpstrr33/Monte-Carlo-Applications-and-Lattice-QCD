'''
These measurements were calculated using 'gluonic_path_integral.py' for the parameters:
N = 8
N_cor = 1
N_cf = 100
N_ma = 10
eps = 0.24
beta = 1.719 (Used for improved action)
beta = 5.5   (Used for unimproved action)
mu0 = 0.797
c1 = 5 / (3 * mu0**4)
c2 = 1 / (12 * mu0**6)
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

## MC measurements for unimproved action of axa, ax2a and ax3a Wilson
## loops, respecively.
NL1 = [0.497344020312, 0.496450260817, 0.494958096286, 0.494851831912, 0.498917593326, \
       0.493026257220, 0.502275343917, 0.495433804047, 0.499467497290, 0.500157602246, \
       0.497136509812, 0.500525429943, 0.499014443871, 0.499960524011, 0.496154358965, \
       0.492016657021, 0.495288238746, 0.497269714332, 0.493477910734, 0.495886898132, \
       0.494197199973, 0.496312901518, 0.500390272546, 0.498456793718, 0.494058546515]
NL2 = [0.257230865019, 0.259437305415, 0.260486221661, 0.254958551552, 0.260595545840, \
       0.253866862411, 0.268636920195, 0.258482366414, 0.261380015244, 0.265880375940, \
       0.258835229511, 0.264018010242, 0.262779466365, 0.263820759928, 0.256596299210, \
       0.251760961519, 0.259317406373, 0.261605963086, 0.254019367478, 0.258375046051, \
       0.258158044474, 0.257034352038, 0.264048352012, 0.261191478641, 0.257626063514]
NL3 = [0.134890394000, 0.136872268899, 0.137835388765, 0.131402345954, 0.137304957551, \
       0.131398595414, 0.143192635285, 0.135145416013, 0.138336154220, 0.142778870681, \
       0.134728630707, 0.141507103734, 0.139066450142, 0.143051946471, 0.135175955561, \
       0.129554267208, 0.135367546173, 0.140491536126, 0.132013554637, 0.135281543474, \
       0.136710547127, 0.134068247136, 0.138976463669, 0.135811364922, 0.133513709576]

## MC measurements for improved action of axa, ax2a and ax3a Wilson
## loops, respectively.
IL1 = [0.546899888398, 0.543201036559, 0.542129088921, 0.546583939921, 0.547355102538, \
       0.547870737400, 0.544058320769, 0.539918331032, 0.546602695633, 0.549825292447, \
       0.541272211417, 0.543882846736, 0.548054637339, 0.547386516233, 0.546470037113, \
       0.545987722040, 0.549673858748, 0.544914445525, 0.544952740120, 0.546440990028, \
       0.542682159446, 0.545803155112, 0.548100001640, 0.548348736616, 0.545896913836]
IL2 = [0.293063985768, 0.287865575506, 0.283915778653, 0.291320887316, 0.289720390085, \
       0.291516017631, 0.285550222287, 0.284037825658, 0.291080401998, 0.294779954715, \
       0.281906032594, 0.287406340107, 0.293845789113, 0.294918945486, 0.289143705607, \
       0.289962297067, 0.293800229343, 0.286782004509, 0.29126972340, 0.2898528116670, \
       0.288159146506, 0.290645930815, 0.293494653100, 0.294051938582, 0.288672869529]
IL3 = [0.157718199215, 0.154986710066, 0.151262000480, 0.155583812793, 0.154418568623, \
       0.155964440616, 0.150108664766, 0.151352231156, 0.156111355575, 0.161014392580, \
       0.148852751579, 0.154647100522, 0.160696711091, 0.160566719173, 0.155638909685, \
       0.156402956469, 0.158001102554, 0.152843069598, 0.156428188140, 0.154572645439, \
       0.155961238107, 0.155655346574, 0.159939853102, 0.159657783581, 0.153999297913]

## N_cor = 50, so just stretching the range of x so to match with the correct
## number of sweeps
xRange = []
for i in range(25):
    xRange.append(50 * i)
 
## Plots the above data into one plot. The unimproved action is represented by
## squares and the improved action is represented by triangles. With the Wilson
## loops represented by:
##      axa  = blue
##      ax2a = green
##      ax3a = red
plot = True
if plot:
    plt.figure()
    plt.ylim(0, 0.8)
    plt.ylabel('Monte Carlo Measurement')
    plt.xlabel('Sweeps')
    plt.plot(xRange, NL1, 'bs', xRange, NL1, 'k', lw = '1.5')
    plt.plot(xRange, NL2, 'gs', xRange, NL2, 'k', lw = '1.5')
    plt.plot(xRange, NL3, 'rs', xRange, NL3, 'k', lw = '1.5')
    
    plt.plot(xRange, IL1, 'b^', xRange, IL1, 'k', lw = '1.5')
    plt.plot(xRange, IL2, 'g^', xRange, IL2, 'k', lw = '1.5')
    plt.plot(xRange, IL3, 'r^', xRange, IL3, 'k', lw = '1.5')
    
    Ax1A = mpatches.Patch(color = 'blue', label = '$a$ x $a$')
    Ax2A = mpatches.Patch(color = 'green', label = '$a$ x $2a$')
    Ax3A = mpatches.Patch(color = 'red', label = '$a$ x $3a$')
    plt.legend(handles = [Ax1A, Ax2A, Ax3A])

## The mean and standard deviation for the above data sets.
print('Unimproved Action:')
print('    Ax1A: %.5f ± %.5f' % (np.average(NL1), np.std(NL1)))
print('    Ax2A: %.5f ± %.5f' % (np.average(NL2), np.std(NL2)))
print('    Ax3A: %.5f ± %.5f' % (np.average(NL3), np.std(NL3)))
print('\n')

print('Improved Action:')
print('    Ax1A: %.5f ± %.5f' % (np.average(IL1), np.std(IL1)))
print('    Ax2A: %.5f ± %.5f' % (np.average(IL2), np.std(IL2)))
print('    Ax3A: %.5f ± %.5f' % (np.average(IL3), np.std(IL3)))
