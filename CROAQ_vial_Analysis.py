import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


peaks = [1, 3, 5]
epr = []
eprstd = []
epi = []
epistd = []
elstd = []

epr_err = []
epi_err = []
el_err = []
fl = []

## PET
#r = [[52, 60], [52, 60], [52, 60], [52, 60]]
#ri = [[50, 58], [49, 58], [50, 57], [50, 57]]

## PET BESSEL
#r = [[50, 60], [50, 60], [50, 60], [50, 63]]
#ri = [[50, 58], [49, 58], [50, 57], [50, 57]]

## REXOLITE
#r = [[52, 60], [52, 61], [52, 61], [52, 61]]
#ri = [[55, 72], [67, 72], [66, 72], [65, 72]]

## Acrylic
#r = [[55, 66], [55, 68], [55, 68], [55, 68]]
#ri = [[50, 62], [55, 70], [55, 70], [55, 70]]

## Acrylic BESSEL
#r = [[55, 60], [58, 62], [52, 62], [48, 58]]
#ri = [[55, 60], [55, 60], [55, 60], [55, 60]]

## Acetal
#r = [[52, 60], [52, 63], [55, 63], [55, 63]]
#ri = [[48, 55], [48, 55], [48, 55], [48, 55]]

## PVDF
r = [[50, 62], [48, 62], [50, 60], [50, 62]]
ri = [[50, 55], [50, 55], [50, 55], [50, 55]]

## KelF
#r = [[53, 63], [52, 62], [53, 63], [52, 62]]
#ri = [[55, 63], [55, 67], [55, 65], [50, 60]]

## Teflon
#r = [[60, 70], [60, 68], [60, 68], [60, 68]]
#ri = [[30, 40], [40, 60], [55, 68], [55, 68]]

## LUNAR SAMPLES HAVE SLIGHTLY SHORTER BASE SO I EXPECT TO NEED TO MOVE DOWN A BIT

## SAMPLE 1
#r = [[50, 59], [50, 60], [48, 59], [50, 60]]
#ri = [[50, 57], [50, 57], [50, 57], [48, 54]]

## SAMPLE 2
#r = [[49, 59], [50, 60], [48, 58], [50, 60]]
#ri = [[50, 57], [49, 58], [50, 57], [49, 56]]

## SAMPLE 3
#r = [[52, 61], [50, 59], [48, 60], [50, 59]]
#ri = [[50, 57], [49, 59], [50, 57], [49, 56]]

## SAMPLE 4
#r = [[48, 58], [48, 58], [47, 56], [47, 58]]
#ri = [[50, 57], [49, 59], [50, 57], [49, 56]]

## SAMPLE 5
#r = [[50, 58], [51, 60], [48, 57], [50, 58]]
#ri = [[49, 60], [49, 59], [50, 57], [50, 58]]

## SAMPLE 5 t2
#r = [[50, 58], [51, 60], [48, 57], [50, 58]]
#ri = [[45, 63], [49, 63], [49, 63], [49, 63]]

## SAMPLE 7
#r = [[51, 61], [51, 63], [50, 59], [50, 59]]
#ri = [[52, 62], [52, 62], [52, 62], [53, 62]]

## SAMPLE 8
#r = [[53, 63], [52, 63], [47, 57], [48, 59]]
#ri = [[53, 63], [53, 63], [53, 63], [53, 63]]

## SAMPLE 9
#r = [[47, 57], [47, 57], [47, 57], [48, 59]]
#ri = [[47, 57], [47, 57], [47, 57], [47, 57]]

## SAMPLE 10
#r = [[48, 58], [48, 58], [47, 57], [47, 56]]
#ri = [[51, 59], [51, 59], [51, 59], [51, 59]]

## FOR LUNAR SAMPLE VIAL 6
#r = [[81, 92], [79, 89], [79, 92], [83, 95]]
#ri = [[80, 94], [80, 95], [80, 90], [80, 90]]

## FOR LUNAR SAMPLE 4 IN UHF
#r = [[35, 55], [35, 55], [35, 55]]
#ri = [[35, 65], [35, 65], [35, 65]]

## FOR LUNAR SAMPLE 5 IN UHF
#r = [[55, 80], [55, 80], [43, 60]]
#ri = [[40, 75], [40, 75], [40, 75]]

## FOR LUNAR SAMPLE 5 IN UHF BESSEL
#r = [[50, 70], [50, 70], [50, 70]]
#ri = [[50, 70], [50, 70], [50, 70]]

## FOR LUNAR SAMPLE 7 IN UHF
#r = [[60, 90], [50, 90], [45, 75]]
#ri = [[40, 85], [55, 110], [50, 90]]


## FOR LUNAR SAMPLE 2 IN 3-5 GHz
#r = [[0,0], [38, 50], [38, 52], [40, 51], [40, 50]]
#ri = [[0,0], [35, 45], [35, 45], [35, 45], [35, 45]]

## FOR LUNAR SAMPLE 1 IN WR159
#r = [[0,0], [40, 65], [40, 65], [35, 47], [40, 50]]
#ri = [[0,0], [35, 50], [25, 33], [27, 37], [35, 45]]

## FOR LUNAR SAMPLE 2 IN WR159
#r = [[0,0], [40, 65], [40, 65], [30, 45], [40, 50]]
#ri = [[0,0], [35, 50], [25, 33], [25, 35], [35, 45]]

## FOR Rexolite IN WR159
#r = [[0,0], [40, 65], [40, 65], [30, 45], [40, 50]]
#ri = [[0,0], [35, 50], [15, 30], [15, 30], [35, 45]]

Plot = True
vialepsr = 2.54
vialepsi = 4.0e-4
for iii in peaks:
    #df = pd.read_csv('../CROAQ/Data/Rexolite_test2_Peak' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak' + str(iii) + '_t2/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_P' + str(iii) + '_082725/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Rexolite_091025_P' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/RexoliteRod_N5230C_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample5_t2_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample2_uwave_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample5_UHF_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample2_WR159_p' + str(iii) + '/trial1.csv')

    #df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/Rexolite_RexVial_WR159_p' + str(iii) + '/trial1.csv')
    df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/PVDF1_RexVial_p' + str(iii) + '/trial1.csv')
    f = np.array(df['Frequency (GHz)'])
    if (iii == 1):
        Q = np.array(df['Q spline'])
    else:
        Q = np.array(df['Q spline'])
    h = np.array(df['height (mm)'])
    Aeff = df['Aeff_E'][0]
    Aeff_s = df['Aeff_Es'][0]
    Vc = df['WG_length (mm)'][0]*df['WG_d (mm)'][0]*df['WG_a (mm)'][0]

    f0 = f[0]
    #print(f0)
    fl.append(f0)
    Q0 = Q[0]
    yr = (f0 - f)/f0
    yi = 1/Q - 1/Q0
    x = Aeff*(h - 0.25*25.4)/(Vc)

    p, V = np.polyfit(x[r[iii//2][0]:r[iii//2][1]], yr[r[iii//2][0]:r[iii//2][1]], 1, cov = True)

    m = p[0]
    b = p[1]
    #b=0
    #print(b)

    #dyr = np.gradient(yr, h)
    #plt.plot(h, dyr, marker = 'o', linestyle = '--')
    #plt.show()

    if (Plot):
        plt.plot(x, yr, marker = 'o', linestyle = '')
        plt.plot(x, m*x + b)
        plt.xlabel('Effective Inserted Volume Fraction, x', fontsize = 14)
        plt.ylabel('Resonant Frequency Shift, yr', fontsize = 14)
        plt.axvline(x[r[iii//2][0]])
        plt.axvline(x[r[iii//2][1]])
        #print("b: ", b)
        plt.show()

    pi, Vi = np.polyfit(x[ri[iii//2][0]:ri[iii//2][1]], yi[ri[iii//2][0]:ri[iii//2][1]], 1, cov = True)

    mi = pi[0]
    bi = pi[1]

    print(bi)


    resr =np.array(yr[r[iii//2][0]:r[iii//2][1]]) - m*np.array(x[r[iii//2][0]:r[iii//2][1]]) - b
    resr_rms = np.sqrt(np.mean(resr**2))

    resi = np.array(yi[ri[iii//2][0]:ri[iii//2][1]]) - mi*np.array(x[ri[iii//2][0]:ri[iii//2][1]]) - bi
    resi_rms = np.sqrt(np.mean(resi**2))
    #print(resi_rms, np.sqrt(Vi[1][1]), pi[1])
    num_resr = np.sqrt(resr_rms**2 + V[1][1])
    num_resi = np.sqrt(resi_rms**2 + Vi[1][1])

    if (Plot):
        plt.plot(x, yi, marker = 'o', linestyle = '')
        plt.plot(x, mi*x + bi)
        plt.axvline(x[ri[iii//2][0]])
        plt.axvline(x[ri[iii//2][1]])
        plt.show()

    epsr = (yr - b)/(2*x)+1
    epsi = (yi - bi)/(4*x)

    

    #plt.errorbar(x[ri[iii//2][0]:ri[iii//2][1]], epsi[ri[iii//2][0]:ri[iii//2][1]], yerr = num_resi/(4*x[ri[iii//2][0]:ri[iii//2][1]]))
    #plt.show()

    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = -0.008*(epsr0 - 1) + 1.014
    elif (iii == 5):
        a = 0*(epsr0 - 1) + 1.014
    elif (iii == 7):
        a = 0.025*(epsr0 - 1) + 0.995'''

    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.001*(epsr0 - 1) + 1.0
    elif (iii == 5):
        a = -0.003*(epsr0 - 1) + 1.01
    elif (iii == 7):
        a = 0.002*(epsr0 - 1) + 1.036'''
    
    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.00*(epsr0 - 1) + 1.00
    elif (iii == 5):
        a = -0.000*(epsr0 - 1) + 1.004
    elif (iii == 7):
        a = 0.00*(epsr0 - 1) + 1.05'''

    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.001*(epsr0 - 1) + 0.996
    elif (iii == 5):
        a = -0.007*(epsr0 - 1) + 1.037
    elif (iii == 7):
        a = 0.008*(epsr0 - 1) + 1.065'''

    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.001*(epsr0 - 1) + 0.995
    elif (iii == 5):
        a = -0.005*(epsr0 - 1) + 1.033
    elif (iii == 7):
        a = 0.021*(epsr0 - 1) + 1.034'''

    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.002*(epsr0 - 1) + 0.994
    elif (iii == 5):
        a = -0.01*(epsr0 - 1) + 1.04
    elif (iii == 7):
        a = 0.019*(epsr0 - 1) + 1.036'''


    ## USE THIS FOR IN VIAL MEASUREMENTS IN 3-5 GHz CAVITY
    '''if (iii == 3):
        a = 1
    elif (iii == 5):
        a = 0.0027*epsr0 + 1.00058
    elif (iii == 7):
        a = 0.0122*epsr0 + 0.9991
    elif (iii == 9):
        a = 0.0273*epsr0 + 0.9901'''




    ## USE THIS FOR THE 1-3 GHz IN VIAL MEASUREMENTS
    if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.003*(epsr0 - 1) + 0.991
    elif (iii == 5):
        a = 0.002*(epsr0 - 1) + 1.012
    elif (iii == 7):
        a = 0.026*(epsr0 - 1) + 1.019


    ## THIS IS TESTING A BESSEL FUNCTION APPROACH TO a for 1-3GHz
    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.0023212721997093715 * epsr0 +  0.9975125566697594
    elif (iii == 5):
        a = 0.004439300980754807 * epsr0 +  0.9952314785943057
    elif (iii == 7):
        a = 0.007701748596717931 * epsr0 +  0.9917318693592986'''


    ## THIS IS TESTING A BESSEL FUNCTION APPROACH TO a for 3-5GHz
    '''if (iii == 3):
        a = 1
        a = 0.008938441495112217 * epsr0 +  0.9935026850317076
    elif (iii == 5):
        a = 0.013236204709964203 * epsr0 +  0.9902083769280944
    elif (iii == 7):
        a = 0.020237909277574656 * epsr0 +  0.9845763339115715
    elif (iii == 9):
        a = 0.030704207475220558 * epsr0 +  0.975549266779338'''


    ## THIS IS TESTING A BESSEL FUNCTION APPROACH TO a for 0.5-1GHz
    '''if (iii == 1):
        a = 1
    elif (iii == 3):
        a = 0.00034451219091939845 * epsr0 +  0.9996328333320491
    elif (iii == 5):
        a = 0.0007703319342376751 * epsr0 +  0.9991772598777916

    if (iii == 1):
        a = 1
            #a = 0.010476595627288565 * epsr0 +  0.9992015473314041
    elif (iii == 3):
        a = 0.01376857195549889 * epsr0 +  0.986188595768463
    elif (iii == 5):
        a = 0.027177560550619484 * epsr0 +  0.9722818782707711
    elif (iii == 7):
        a = 0.062372236578592394 * epsr0 +  0.9342106244235927'''

    #a = 1
    epsr = (yr - b)/(2*x*a)+1
    epsi = (yi - bi)/(4*x*a)

    epsr = (epsr - vialepsr)*Aeff/Aeff_s + vialepsr
    epsi = (epsi - vialepsi*vialepsr)*Aeff/Aeff_s + vialepsi*vialepsr

    #print(np.round(np.mean(num_resr/(4*a*x[r[iii//2][0]:r[iii//2][1]])), 7), np.round(np.mean(num_resi/(4*a*x[ri[iii//2][0]:ri[iii//2][1]])), 7))
    epr_err.append(np.round(np.mean(num_resr/(2*a*x[r[iii//2][0]:r[iii//2][1]])), 5))
    epi_err.append(np.round(np.mean(num_resi/(4*a*x[ri[iii//2][0]:ri[iii//2][1]])), 7))

    '''fig, axs = plt.subplots(1, 1)
    axs.plot(epsr, marker = 'o', linestyle = '')
    axs.axhline(np.mean(epsr[r[iii//2][0]: r[iii//2][1]]), color = 'k', linestyle = '--')
    #axs[0].set_xlim(0.00005, np.max(x))
    axs.set_xlim(r[iii//2][0], r[iii//2][1])
    axs.set_ylim(2.5,3.5)
    axs.axvline(r[iii//2][0], linestyle = '--')
    axs.axvline(r[iii//2][1], linestyle = '--')
    axs.set_xlabel("Step", fontsize = 15)
    axs.set_ylabel("Real Permittivity", fontsize = 15)
    axs.tick_params(axis='both', which='major', labelsize=12)
    axs.grid()
    
    #axs[1].plot(x, epsi, marker = 'o', linestyle = '')
    #axs[1].set_xlim(0.00005, np.max(x))
    #axs[1].axvline(x[ri[iii//2][0]], linestyle = '--')
    #axs[1].axvline(x[ri[iii//2][1]], linestyle = '--')
    #axs[1].set_ylim(1e-4, 0.1)
    #axs[1].set_yscale('log')
    #axs[1].grid(which = 'both')
    plt.show()'''

    #plt.plot(x[20:], epsr[20:], marker = 'o', linestyle = '')
    #plt.show()
    #plt.plot(x[20:], epsi[20:], marker = 'o', linestyle = '')
    #plt.show()
    if (iii == 1):
        #a = 1
        epsr0 = np.mean(epsr[r[iii//2][0]:r[iii//2][1]])
    epr.append(np.mean(epsr[r[iii//2][0]:r[iii//2][1]]))
    eprstd.append(np.std(epsr[r[iii//2][0]:r[iii//2][1]]))
    #print('epsr std: ', np.std(epsr[r[iii//2][0]:r[iii//2][1]]))
    epistd.append(np.std(epsi[ri[iii//2][0]:ri[iii//2][1]]))
    epi.append(np.mean(epsi[ri[iii//2][0]:ri[iii//2][1]])/np.mean(epsr[r[iii//2][0]:r[iii//2][1]]))

    el_err.append(np.mean(epsi[ri[iii//2][0]:ri[iii//2][1]])/np.mean(epsr[r[iii//2][0]:r[iii//2][1]])*np.sqrt((np.mean(num_resr/(4*a*x[r[iii//2][0]:r[iii//2][1]]))/np.mean(epsr[r[iii//2][0]:r[iii//2][1]]))**2 + (np.mean(num_resi/(4*a*x[ri[iii//2][0]:ri[iii//2][1]]))/np.mean(epsi[ri[iii//2][0]:ri[iii//2][1]]))**2))
    elstd.append(np.mean(epsi[ri[iii//2][0]:ri[iii//2][1]])/np.mean(epsr[r[iii//2][0]:r[iii//2][1]])*np.sqrt((np.std(epsr[r[iii//2][0]:r[iii//2][1]])/np.mean(epsr[r[iii//2][0]:r[iii//2][1]]))**2 + (np.std(epsi[ri[iii//2][0]:ri[iii//2][1]])/np.mean(epsi[ri[iii//2][0]:ri[iii//2][1]]))**2))
    #print('loss std: ', np.std(epsi[ri[iii//2][0]:ri[iii//2][1]])/np.mean(epsr[r[iii//2][0]:r[iii//2][1]]))
    #print('Real Permittivity', np.mean(epsr[r[0]:r[1]]))
    #print('loss tangent', np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]]))
    #print("")

    #plt.plot(x, epsr, marker = 'o', linestyle = '')
    #plt.xlabel('Frequency (GHz)', fontsize = 15)
    #plt.ylabel('Real Permittivity', fontsize = 15)
    #plt.title('$TE_{101}$, Rexolite Rod', fontsize = 15)
    #plt.axhline(2.54)
    #plt.axvline(x[r[0]])
    #plt.axvline(x[r[1]])
    #plt.grid()
    #plt.ylim(1, 10)

    #plt.show()

    #plt.plot(x, np.array(epsi)/np.array(epsr), marker = 'o', linestyle = '')
    #plt.xlabel('Frequency (GHz)', fontsize = 15)
    #plt.ylabel('loss tangent', fontsize = 15)
    #plt.title('Rexolite Rod', fontsize = 15)
    #plt.axhline(4.5e-4)
    #plt.grid()
    #plt.axvline(x[ri[0]])
    #plt.axvline(x[ri[1]])
    #plt.yscale('log')
    #plt.ylim(1e-5, 1e-3)

    #plt.show()

fig, axs = plt.subplots(2, 1)
axs[0].plot(fl, epr, marker = 'o', linestyle = '--')
axs[0].set_xlabel('Frequency (GHz)', fontsize = 15)
axs[0].set_ylabel('Real permitivity', fontsize = 15)
axs[0].set_title('Lunar Sample', fontsize = 15)
axs[0].grid()
axs[0].axhline(2.54)
axs[0].set_ylim(1,5)

#print(epi)
axs[1].plot(fl, epi, marker = 'o', linestyle = '--')
axs[1].set_xlabel('Frequency (GHz)', fontsize = 15)
axs[1].set_ylabel('Loss Tangent', fontsize = 15)
#axs[1].set_title('Lunar Sample', fontsize = 15)
axs[1].grid()
axs[1].axhline(4.0e-4)
axs[1].set_ylim(1e-4,1)
axs[1].set_yscale('log')

plt.tight_layout()
plt.show()
print('Real: ', np.round(np.array(epr), 5), np.round(np.array(eprstd), 5))
print(np.round(np.array(epr_err), 5))
print('Imag: ', np.round(np.array(epi)*np.array(epr), 5), np.round(np.array(epistd), 7))
print(np.round(np.array(epi_err), 7))
print('Loss: ', np.round(np.array(epi), 7), np.round(np.array(elstd), 8))
print(np.round(np.array(el_err), 8))
#plt.plot(fl, np.array(epi)/np.array(epr), marker = 'o', linestyle = '--')
#plt.xlabel('Frequency (GHz)', fontsize = 15)
#plt.ylabel('Imag permitivity', fontsize = 15)
#plt.title('Rexolite Rod', fontsize = 15)
#plt.grid()
#plt.yscale('log')
#plt.axhline(4.5e-4)
#plt.ylim(1e-5,1)
#plt.show()
