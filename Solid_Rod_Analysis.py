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
from Sample_Analysis_w_correction_func_V3 import cavity_analysis
import warnings
## This is just to supress the divide by 0 warning
warnings.filterwarnings("ignore")


Plastics = ['Teflon', 'KelF', 'Rexolite', 'Acrylic', 'PVDF', 'Acetal', 'PET']
#Plastics = ['Rexolite']
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
        df = pd.read_csv('./Measurements/' + Plastics[ii] + '/' + Plastics[ii] + 'Rod1_PostBake_Peak' + str(modes[i]) + '.csv')
        #ranges2 = [[12,19], [12, 19], [8, 17], [7, 12]] ##[[45%, 85%], [45%, 82%], [23%, 72%], [18%, 45%]] where these are percent of sample insertion not vial insertion
        xdata1, real_data1, imag_data1, loss_tan1 = cavity_analysis(df, df['Aeff'][0], 
                                                                        plot = False, 
                                                                        real_range = real_range[i], imag_range = imag_range[i], Solid = True, Corr = True)
        eps_r_solid[ii,i] = real_data1[2]
        eps_i_solid[ii,i] = imag_data1[2]
        #if (modes[i] == 3):
            #print(real_data1[2])

        loss_solid[ii,i] = loss_tan1[1]
        b[ii, i] = real_data1[5]
        m[ii, i] = real_data1[6]
        mm[ii, i] = 2*(real_data1[2] - 1)
        #print(m[ii,i], mm[ii, i])
        epsi_RMSE[ii,i] = imag_data1[3]
        epsr_RMSE[ii, i] = real_data1[4]
        loss_RMSE[ii, i] = loss_tan1[2]
    print(Plastics[ii], np.round(eps_r_solid[ii, :],3))
    print(np.round(epsr_RMSE[ii, :], 4))
    print(np.round(eps_i_solid[ii, :], 5))
    print(np.round(epsi_RMSE[ii, :], 6))
    print(np.round(loss_solid[ii, :], 7))
    print(np.round(loss_RMSE[ii, :], 7))
    #print(np.round(real_data1[7], 3))
    axs[0].errorbar(modes, np.array(eps_r_solid[ii, :]), yerr = epsr_RMSE[ii, :], marker = 'o', linestyle = '-', label = 'Solid ' + Plastics[ii] + ' Rod1')
    axs[1].plot(modes, loss_solid[ii, :], marker = 'o', linestyle = '-', label = 'Solid ' + Plastics[ii] + ' Rod1')
    #axs[2].plot(modes, (np.array(eps_r_solid) - 1)**(-1), marker = 'o', linestyle = '--')



#axs[0].plot([3, 5, 7, 9], np.divide(np.array(eps_r),np.array(eps_r_Teflon)), marker = 'o', linestyle = '--', color = 'b')
#axs[0].legend()
axs[0].set_xlabel('Mode', fontsize = 12)
axs[0].set_ylabel('Measured Permittivity', fontsize = 12)
axs[0].axhline(1, linestyle = '', color = 'k')
axs[0].axhline(2.54, linestyle = '--', color = 'g', label = 'Canonical Rexolite')
axs[0].axhline(2.055, linestyle = '--', color = 'b', label = 'Canonical Teflon')
axs[0].axhline(2.375, linestyle = '--', color = 'orange', label = 'Canonical KelF')
axs[0].axhline(2.634, linestyle = '--', color = 'r', label = 'Canonical Acrylic')
axs[0].grid()
axs[1].set_xlabel('Mode', fontsize = 12)
axs[1].set_ylabel('Measured tan($\delta$)', fontsize = 12)
#axs[1].plot([3, 5, 7, 9], loss_Teflon, marker = 'o', linestyle = '--', color = 'r', label = 'vial7 Teflon Rod1')
axs[0].legend(bbox_to_anchor=(1.0, 0.4), prop = {'size': 10}, ncol = 2)
axs[1].grid()
#axs[0].set_xscale('log')
axs[1].set_yscale('log')
axs[1].axhline(4.5e-4, linestyle = '--', color = 'g')
axs[1].axhline(2e-4, linestyle = '--', color = 'b')
axs[1].axhline(3e-3, linestyle = '--', color = 'orange')
axs[1].axhline(7.2e-3, linestyle = '--', color = 'r')
plt.suptitle('Measured Permittivty and Loss Tangent for Solid Plastic Samples', fontsize = 15, fontweight = 'bold')

canon = np.array([2.055, 2.375, 2.54, 2.634])
plt.show()

for i in range(len(Plastics[0:4])):
    plt.plot(modes, np.abs((eps_r_solid[i, :] - canon[i]))/canon[i]*100, marker = 'o', linestyle = '--', label = Plastics[i])
plt.xlabel('Mode', fontsize = 12)
plt.ylabel('Percent Error [%]', fontsize = 12)
plt.title('Percent error from canonical value using Eigenmode Approach', fontsize = 12, fontweight = 'bold')
plt.grid()
plt.show()

'''plt.plot(modes, m.transpose())
plt.plot(modes, mm.transpose(), linestyle = '--')
plt.show()'''
#print(epsr_RMSE)
'''for ii in range(len(Plastics)):
    plt.plot(modes, epsr_RMSE[ii, :], marker ='o', label = Plastics[ii] + ': RMSE', color = colors[ii])
    plt.plot(modes, del_eps_r[ii, :],marker = 'D', linestyle = '--', label = Plastics[ii] + ': Propogated Error', color = colors[ii])
plt.xlabel('Mode')
plt.ylabel('RMSE')
plt.legend()
plt.show()'''


