'''
These measurements were calculated using 'gluonic_path_integral.py' for the
parameters found in fimp_actes parameters_improved.txt and parameters_unimproved.txt
for the improved and unimproved action, respectively. This program plots the 
data as shown below and prints some basic stats.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


'''
Reads from the fimp_actes created by 'gluonic_path_integral.py' and stores them in
the arrays imp_act and unimp_act for the three Wimp_actson loops for the
improved and unimproved action, respectively. It also reads the total number of 
measurements done for both actions to be used in the plot.
'''
def read_data():
    imp_act = []
    unimp_act = []
    for i in range(1, 4):
        imp_act.append(open('Imp_Ax%sA_Measurements.txt' % i,
                            'r').read().split(','))
        unimp_act.append(open('NoImp_Ax%sA_Measurements.txt' % i,
                              'r').read().split(','))

        a = []
        b = []
        for j in imp_act[i - 1][-1]:
            if j != 'A':
                a.append(j)
            else:
                break
        for j in unimp_act[i - 1][-1]:
            if j != 'A':
                b.append(j)
            else:
                break
        imp_act[i - 1][-1] = ''.join(a)
        unimp_act[i - 1][-1] = ''.join(b)

        for j in range(len(imp_act[i - 1])):
            imp_act[i - 1][j] = float(imp_act[i - 1][j])
        for j in range(len(unimp_act[i - 1])):
            unimp_act[i - 1][j] = float(unimp_act[i - 1][j])

    ## Reads the parameters to find the N_cor value
    imp_param = open('parameters_improved.txt', 'r').read().split('\n')[:-1]
    unimp_param = open('parameters_unimproved.txt', 'r').read().split('\n')[:-1]

    for i in imp_param:
        if i[0:7] == 'N_cor =':
            N_cor_imp = int(i[8:])
        elif i[0:6] == 'N_cf =':
            N_cf_imp = int(i[7:])
    for i in unimp_param:
        if i[0:7] == 'N_cor =':
            N_cor_unimp = int(i[8:])
        elif i[0:6] == 'N_cf =':
            N_cf_unimp = int(i[7:])

    return imp_act, unimp_act, N_cor_imp, N_cor_unimp, N_cf_imp, N_cf_unimp


'''
Calculates the range for each set of data points since each point is N_cor sweeps
from the others and there are a total of N_cf points. Then plots the above data
on plot differentiating the actions by shape (squares for unimproved and
triangles for improved) and the loops by color (blue for AxA, green for Ax2A
and red for Ax3A).
'''
def plot_data(imp_act, unimp_act, N_cor_imp, N_cor_unimp, N_cf_imp, N_cf_unimp):
    x_range_imp = []
    x_range_unimp = []
    
    for i in range(N_cf_imp):
        x_range_imp.append(N_cor_imp * i)
    for i in range(N_cf_unimp):
        x_range_unimp.append(N_cor_unimp * i)

    plt.figure()
    plt.ylim(0, 0.8)
    plt.ylabel('Monte Carlo Measurement')
    plt.xlabel('Sweeps')

    plt.plot(x_range_imp, imp_act[0], 'b^',
             x_range_imp, imp_act[0], 'k', lw = '1.5')
    plt.plot(x_range_imp, imp_act[1], 'g^',
             x_range_imp, imp_act[1], 'k', lw = '1.5')
    plt.plot(x_range_imp, imp_act[2], 'r^',
             x_range_imp, imp_act[2], 'k', lw = '1.5')

    plt.plot(x_range_unimp, unimp_act[0], 'bs',
             x_range_unimp, unimp_act[0], 'k', lw = '1.5')
    plt.plot(x_range_unimp, unimp_act[1], 'gs',
             x_range_unimp, unimp_act[1], 'k', lw = '1.5')
    plt.plot(x_range_unimp, unimp_act[2], 'rs',
             x_range_unimp, unimp_act[2], 'k', lw = '1.5')

    ## Creates the legend
    Ax1A = mpatches.Patch(color = 'blue', label = '$a$ x $a$')
    Ax2A = mpatches.Patch(color = 'green', label = '$a$ x $2a$')
    Ax3A = mpatches.Patch(color = 'red', label = '$a$ x $3a$')
    plt.legend(handles = [Ax1A, Ax2A, Ax3A])


'''
Main function
'''
def main():
    ## Read the data
    imp_act, unimp_act, N_cor_imp, N_cor_unimp, N_cf_imp, N_cf_unimp = read_data()

    ## Plot the data
    plot_data(imp_act, unimp_act, N_cor_imp, N_cor_unimp, N_cf_imp, N_cf_unimp)

    ## The mean and standard deviation for the above data sets.
    print('Unimproved Action:')
    for i in range(3):
        print('    Ax1A: %.5f ± %.5f' % (np.average(unimp_act[i]),
                                         np.std(unimp_act[i])))
    print('\n')
    print('Improved Action:')
    for i in range(3):
        print('    Ax1A: %.5f ± %.5f' % (np.average(imp_act[i]),
                                         np.std(imp_act[i])))

if __name__ == "__main__":
    main()