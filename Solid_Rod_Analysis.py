###################################################################################################################
#
#   Code: Sample_Analysis_SampleFilledVial_V2.py
#
#   Programmer: Payton Linton (linton.93@osu.edu)
#
#   Purpose: This code uses the data for the frequency shift for a solid rod with no holder.
#            Since there is no holder, there are no corrections that need to be made.
#            It should noted that as off now the effective volume has not been calculated and
#            instead I have been hard ocding it in for the different rods I was looking at.
#            This should at some point be fixed. But now is not that time. 
#
#
#
###############################################################################################################



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Anal_func_V3 import cavity_analysis
import warnings
## This is just to supress the divide by 0 warning
warnings.filterwarnings("ignore")


Plastics = ['Rexolite']
modes = [3, 5, 7, 9]
colors = ['b', 'orange', 'g', 'r', 'm', 'k', 'c']
epsr_RMSE = np.zeros((len(Plastics), len(modes)))
epsi_RMSE = np.zeros((len(Plastics), len(modes)))
loss_RMSE = np.zeros((len(Plastics), len(modes)))
eps_r_solid = np.zeros((len(Plastics), len(modes)))
eps_i_solid = np.zeros((len(Plastics), len(modes)))
loss_solid = np.zeros((len(Plastics), len(modes)))
b = np.zeros((len(Plastics), len(modes)))
m = np.zeros((len(Plastics), len(modes)))
mm = np.zeros((len(Plastics), len(modes)))
ranges = [[10,19], [9, 18], [7,16], [5,14]]
fig, axs = plt.subplots(1,2)

real_range = [[11,19], [11,19], [11,19], [11,19]] ## Teflon
imag_range = [[6,12], [4,16], [9,19], [9,16]] ## Teflon

real_range = [[11,19], [11,19], [11,19], [11,19]] ## Rexolite
imag_range = [[7,15], [4,16], [13,19], [9,16]] ## Rexolite
for ii in range(len(Plastics)):
    #eps_r_solid = []
    #bruh = [[9,18], [9,16], [9,14], [9,12]]
    for i in range(len(modes)):
        #print(modes[i])
        df = pd.read_csv('./Example Data/Solid Rexolite/' + Plastics[ii] + 'Rod1_PostBake_Peak' + str(modes[i]) + '.csv')

        xdata1, real_data1, imag_data1, loss_tan1 = cavity_analysis(df, df['Aeff'][0], 
                                                                        plot = False, 
                                                                        real_range = real_range[i], imag_range = imag_range[i], Solid = True, Corr = True)
        eps_r_solid[ii,i] = real_data1[2]
        eps_i_solid[ii,i] = imag_data1[2]

        loss_solid[ii,i] = loss_tan1[1]
        b[ii, i] = real_data1[5]
        m[ii, i] = real_data1[6]
        mm[ii, i] = 2*(real_data1[2] - 1)
        epsi_RMSE[ii,i] = imag_data1[3]
        epsr_RMSE[ii, i] = real_data1[4]
        loss_RMSE[ii, i] = loss_tan1[2]

    print(Plastics[ii], np.round(eps_r_solid[ii, :],3))
    print(np.round(epsr_RMSE[ii, :], 4))
    print(np.round(eps_i_solid[ii, :], 5))
    print(np.round(epsi_RMSE[ii, :], 6))
    print(np.round(loss_solid[ii, :], 7))
    print(np.round(loss_RMSE[ii, :], 7))

    axs[0].errorbar(modes, np.array(eps_r_solid[ii, :]), yerr = epsr_RMSE[ii, :], marker = 'o', linestyle = '-', label = 'Solid ' + Plastics[ii] + ' Rod1')
    axs[1].plot(modes, loss_solid[ii, :], marker = 'o', linestyle = '-', label = 'Solid ' + Plastics[ii] + ' Rod1')


axs[0].set_xlabel('Mode', fontsize = 12)
axs[0].set_ylabel('Measured Permittivity', fontsize = 12)
axs[0].axhline(1, linestyle = '', color = 'k')
axs[0].axhline(2.54, linestyle = '--', color = 'g', label = 'Canonical Rexolite')
axs[0].grid()
axs[1].set_xlabel('Mode', fontsize = 12)
axs[1].set_ylabel('Measured tan($\delta$)', fontsize = 12)
#axs[1].plot([3, 5, 7, 9], loss_Teflon, marker = 'o', linestyle = '--', color = 'r', label = 'vial7 Teflon Rod1')
axs[0].legend(bbox_to_anchor=(1.0, 0.4), prop = {'size': 10}, ncol = 2)
axs[1].grid()
#axs[0].set_xscale('log')
axs[1].set_yscale('log')
axs[1].axhline(4.5e-4, linestyle = '--', color = 'g')

plt.suptitle('Measured Permittivty and Loss Tangent for Solid Plastic Samples', fontsize = 15, fontweight = 'bold')
plt.show()


