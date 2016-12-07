'''
These measurements were calculated using 'gluonic_path_integral.py' for the
parameters found in fimp_actes parameters_improved.txt and parameters_unimproved.txt
for the improved and unimproved action, respectively. This program plots the 
data as shown below and prints some basic stats.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


'''
Reads from the files created by 'gluonic_path_integral.py' and stores them
in the arrays imp_act and unimp_act, depending on the action used, for the three
Wilson loops. It also reads the parameters used for the improved and unimproved
action which are stored as dictionaries.
'''
def read_data():
    par_imp = np.load('parameters_improved.npy')[()]
    par_unimp = np.load('parameters_unimproved.npy')[()]

    imp_act = np.array([np.load('imp_Ax1A_measurements.npy'),
                        np.load('imp_Ax2A_measurements.npy'),
                        np.load('imp_Ax3A_measurements.npy')])

    unimp_act = np.array([np.load('unimp_Ax1A_measurements.npy'),
                          np.load('unimp_Ax2A_measurements.npy'),
                          np.load('unimp_Ax3A_measurements.npy')])

    return par_imp, par_unimp, imp_act, unimp_act


'''
Calculates the range for each set of data points since each point is N_cor
sweeps from the others and there are a total of N_cf points. Then plots the
above data on plot differentiating the actions by shape (squares for unimproved
and triangles for improved) and the loops by color (blue for AxA, green for
Ax2A and red for Ax3A).
'''
def plot_data(imp_act, unimp_act, N_cor_imp, N_cor_unimp, N_cf_imp, N_cf_unimp):
    x_range_imp = []
    x_range_unimp = []

    for i in range(N_cf_imp):
        x_range_imp.append(N_cor_imp * i)
    for i in range(N_cf_unimp):
        x_range_unimp.append(N_cor_unimp * i)

    plt.figure()
    plt.ylim(0, 1)
    plt.ylabel('Monte Carlo Measurement')
    plt.xlabel('Sweeps')

    ## Only plots if there is data (data isn't an empty set)
    if imp_act != []:
        plt.plot(x_range_imp, imp_act[0], 'b^',
                 x_range_imp, imp_act[0], 'k', lw = '1.5')
        plt.plot(x_range_imp, imp_act[1], 'g^',
                 x_range_imp, imp_act[1], 'k', lw = '1.5')
        plt.plot(x_range_imp, imp_act[2], 'r^',
                 x_range_imp, imp_act[2], 'k', lw = '1.5')

    if unimp_act != []:
        plt.plot(x_range_unimp, unimp_act[0], 'bs',
                 x_range_unimp, unimp_act[0], 'k', lw = '1.5')
        plt.plot(x_range_unimp, unimp_act[1], 'gs',
                 x_range_unimp, unimp_act[1], 'k', lw = '1.5')
        plt.plot(x_range_unimp, unimp_act[2], 'rs',
                 x_range_unimp, unimp_act[2], 'k', lw = '1.5')

    ## Creates the legend
    legend = []
    legend.append(mpatches.Patch(color = 'blue', label = '$a$ x $a$'))
    legend.append(mpatches.Patch(color = 'green', label = '$a$ x $2a$'))
    legend.append(mpatches.Patch(color = 'red', label = '$a$ x $3a$'))
    if imp_act != []:
        legend.append(mlines.Line2D([], [], marker = '^', color = 'grey',
                               markersize=8, linestyle = '-',
                               label = 'Improved Action'))
    if unimp_act != []:
        legend.append(mlines.Line2D([], [], marker = 's', color = 'grey',
                                 markersize=8, linestyle = '-',
                                 label = 'Unimproved Action'))
    plt.legend(handles = legend)


'''
Main function
'''
def main():
    ## Read the data
    par_imp, par_unimp, imp_act, unimp_act = read_data()

    ## Plot the data
    plot_data(imp_act, unimp_act, par_imp['N_COR'], par_unimp['N_COR'],
              par_imp['N_CF'], par_unimp['N_CF'])

    ## The mean and standard deviation for the above data sets if they exist
    if unimp_act != []:
        print('Unimproved Action:')
        for i in range(3):
            print('    Ax%sA: %.5f ± %.5f' % ((i + 1), np.average(unimp_act[i]),
                            np.std(unimp_act[i]) / np.sqrt(par_unimp['N_CF'])))
        print('\n')
    if imp_act != []:
        print('Improved Action:')
        for i in range(3):
            print('    Ax%sA: %.5f ± %.5f' % ((i + 1), np.average(imp_act[i]),
                            np.std(imp_act[i]) / np.sqrt(par_imp['N_CF'])))


if __name__ == "__main__":
    main()