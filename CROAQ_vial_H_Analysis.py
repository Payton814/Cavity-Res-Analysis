import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


peaks = [2,4, 6]
r = [35, 50]
ri = [35, 50]

e = [[3.98002, 3.97026, 3.99693, 3.97222],
     [3.69846, 3.67172, 3.70142, 3.68187],
     [3.26016, 3.28068, 3.29367, 3.2756],
     [3.27369, 3.28166, 3.2883, 3.26116],
     [3.61046, 3.46686, 3.61109, 3.5784 ],
     [2.93017, 2.94165, 2.94362, 2.92943 ],
     [3.27531, 3.2738, 3.28498, 3.2739],
     [3.75994, 3.7404, 3.75054, 3.74953 ],
     [3.78218, 3.77593, 3.79879, 3.79155 ],
     [3.62115, 3.6301, 3.67086, 3.63661]]

el = [[0.01153, 0.01132, 0.01151, 0.01201],
      [0.0104 , 0.00977, 0.01008, 0.01041],
      [0.00331, 0.00348, 0.0039 , 0.00437],
      [0.00331, 0.00358, 0.00399, 0.00446],
      [0.00181, 0.00192, 0.0022 , 0.00267],
      [0.00161, 0.00173, 0.00188, 0.00195],
      [0.00242, 0.00269, 0.00324, 0.0038],
      [0.00141, 0.00163, 0.00194, 0.00239 ],
      [0.00149, 0.00159, 0.00171, 0.00201],
      [0.0016,  0.00175, 0.00196, 0.00241 ]]


#e = [3.96, 3.98, 4.00, 4.04]
#el = [0.0113, 0.0113, 0.0115, 0.0122]

#e = [3.67, 3.67, 3.69, 3.74]
#el = [0.0098, 0.0097, 0.0098, 0.0010]

#e = [3.62, 3.65, 3.68, 3.73]
#el = [0.00143, 0.00167, 0.00196, 0.00235]

#e = [3.24, 3.27, 3.34, 3.39]
#el = [0.00317, 0.00342, 0.00382, 0.00436]

#e = [3.275, 3.272, 3.345, 3.527]
#el = [0.00349, 0.00359, 0.00394, 0.00448]


## UHF SAMPLE
#e = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[3.59,3.6,0,0], [3.425, 3.427, 0, 0], [0,0,0,0], [3.2, 3.23, 0, 0]]
#el = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0.003,0.0032,0,0], [0,0,0,0], [0,0,0,0], [0.002324,  0.0020137, 0, 0]]

## 3-5 GHz Cavity
e = [[0,0,0,0], [3.76452, 3.77799, 3.78271, 3.79827]]
el = [[0,0,0,0],[0.0106476, 0.0107805, 0.0113982, 0.0117952]]

e = [[3.62,3.62,3.62,0]]
el = [[0.0101,0.152,0.152,0]]

Efield_Corr = True
vialepsr = 2.54
viallosstan = 4.0e-4
fig, axs = plt.subplots(2, 1)
#fig2, axs2 = plt.subplots(1, 2)
samplenum = [1]
for ii in samplenum:
    epr = []
    epi = []
    fl = []
    eprstd = []
    epistd = []
    elstd = []
    epr_err = []
    epi_err = []
    el_err = []
    for iii in peaks:
        if (iii == 2):
            r = [32, 50]
            ri = [32, 50]
            ## FOR UHF
            #r = [25, 85]
            #ri = [30, 65]

            ## FOR UHF SAMPLE 7
            #r = [30, 80]
            #ri = [40, 110]

            ## FOR SAMPLE 1 WR195
            r = [30, 65]
            ri = [30, 65]
        else:
            r = [25, 50]
            ri = [25, 50]

            ## FOR UHF
            #r = [30, 85]
            #ri = [30, 75]

            ## FOR UHF SAMPLE 7
            #r = [40, 80]
            #ri = [40, 110]

            ## FOR SAMPLE 1 WR195
            r = [30, 65]
            ri = [20, 40]

        if (ii == 6):
            r = [50, 95]
            ri = [65, 95]
        #df = pd.read_csv('../CROAQ/Data/Rexolite_test2_Peak' + str(iii) + '/trial1.csv')
        #df = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak' + str(iii) + '_t2/trial1.csv')
        #df = pd.read_csv('../CROAQ/Data/Rexolite_P' + str(iii) + '_082725/trial1.csv')
        #df = pd.read_csv('./Example Data/CROAQ_Data/Rexolite_091025_P' + str(iii) + '/trial1.csv')
        #df = pd.read_csv('./Example Data/CROAQ_Data/RexoliteRod_N5230C_p' + str(iii) + '/trial1.csv')
        #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample' + str(ii) + '_p' + str(iii) + '/trial1.csv')
        #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample' + str(ii) + '_uwave_p' + str(iii) + '/trial1.csv')
        df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample' + str(ii) + '_WR159_p' + str(iii) + '/trial1.csv')
        #df = pd.read_csv('./Example Data/CROAQ_Data/LunarSample' + str(ii) + '_UHF_p' + str(iii) + '/trial1.csv')
        f = np.array(df['Frequency (GHz)'])
        Q = np.array(df['Q spline'])
        h = np.array(df['height (mm)'])
        Aeff = df['Aeff_E'][0]
        Aeff_s = df['Aeff_Es'][0]
        Aeff_H = df['Aeff_H'][0]
        Aeff_Hs = df['Aeff_Hs'][0]
        Vc = df['WG_length (mm)'][0]*df['WG_d (mm)'][0]*df['WG_a (mm)'][0]
        print(f[0])
        f0 = f[0]
        fl.append(f0)
        Q0 = Q[0]
        yr = (f0 - f)/f0
        yi = 1/Q - 1/Q0
        x = Aeff_Hs*(h - 0.25*25.4)/(Vc)

        p, V = np.polyfit(x[r[0]:r[1]], yr[r[0]:r[1]], 1, cov = True)

        m = p[0]
        b = p[1]
        #print(b)

        dyr = np.gradient(yr, h)
        #plt.plot(h, dyr, marker = 'o', linestyle = '--')
        #plt.show()
        fig2, axs2 = plt.subplots(1, 2)
        axs2[0].plot(x, yr, marker = 'o', linestyle = '')
        axs2[0].plot(x, m*x + b)
        axs2[0].axvline(x[r[0]])
        axs2[0].axvline(x[r[1]])
        #print("b: ", b)

        pi, Vi = np.polyfit(x[ri[0]:ri[1]], yi[ri[0]:ri[1]], 1, cov = True)

        mi = pi[0]
        bi = pi[1]

        axs2[1].plot(x, yi, marker = 'o', linestyle = '')
        axs2[1].plot(x, mi*x + bi)
        axs2[1].axvline(x[ri[0]])
        axs2[1].axvline(x[ri[1]])
        plt.show()

        resr =np.array(yr[r[0]:r[1]]) - m*np.array(x[r[0]:r[1]]) - b
        resr_rms = np.sqrt(np.mean(resr**2))

        resi = np.array(yi[ri[0]:ri[1]]) - mi*np.array(x[ri[0]:ri[1]]) - bi
        resi_rms = np.sqrt(np.mean(resi**2))
        #print(resi_rms, np.sqrt(Vi[1][1]), pi[1])
        num_resr = np.sqrt(resr_rms**2 + V[1][1])
        num_resi = np.sqrt(resi_rms**2 + Vi[1][1])


        ## USE THIS FOR IN VIAL MEASUREMENTS IN 3-5 GHz CAVITY
        if (iii == 2):
            a = 1
        elif (iii == 4):
            a = 0.0027*e[ii - 1][iii//2-1] + 1.00058
        elif (iii == 6):
            a = 0.0122*e[ii - 1][iii//2-1] + 0.9991
        elif (iii == 8):
            a = 0.0273*e[ii - 1][iii//2-1] + 0.9901

        #a = 1
        epsr = (yr - b)/(2*x*a)+1
        epsi = (yi - bi)/(4*x*a)

        if (Efield_Corr == False):
            e[ii - 1][iii//2-1] = 0
            el[ii - 1][iii//2-1] = 0
        print(e[ii - 1][iii//2-1])
        epsr = epsr - (vialepsr - 1)*Aeff/Aeff_Hs - (e[ii - 1][iii//2-1] - vialepsr)*Aeff_s/Aeff_Hs
        epsi = epsi - viallosstan*vialepsr*Aeff/Aeff_Hs - (e[ii - 1][iii//2-1]*el[ii - 1][iii//2-1] - viallosstan*vialepsr)*Aeff_s/Aeff_Hs

        if (iii == 1):
            a = 1
            epsr0 = np.mean(epsr[r[0]:r[1]])
        epr.append(np.mean(epsr[r[0]:r[1]]))
        epi.append(np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]]))
        eprstd.append(np.std(epsr[r[0]:r[1]]))
        epistd.append(np.std(epsi[ri[0]:ri[1]]))
        elstd.append(np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]])*np.sqrt((np.std(epsr[r[0]:r[1]])/np.mean(epsr[r[0]:r[1]]))**2 + (np.std(epsi[ri[0]:ri[1]])/np.mean(epsi[ri[0]:ri[1]]))**2))

        epr_err.append(np.round(np.mean(num_resr/(2*a*x[r[0]:r[1]])), 5))
        epi_err.append(np.round(np.mean(num_resi/(4*a*x[ri[0]:ri[1]])), 7))
        el_err.append(np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]])*np.sqrt((np.mean(num_resr/(4*a*x[r[0]:r[1]]))/np.mean(epsr[r[0]:r[1]]))**2 + (np.mean(num_resi/(4*a*x[ri[0]:ri[1]]))/np.mean(epsi[ri[0]:ri[1]]))**2))
        #print('Real Permeability', np.mean(epsr[r[0]:r[1]]))
        #print('Imaginary Permeability', np.mean(epsi[ri[0]:ri[1]]))
        #print('loss tangent', np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]]))

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

    #fig, axs = plt.subplots(2, 1)
    axs[0].plot(fl, np.round(epr, 3), marker = 'o', linestyle = '--', label = ii)
    axs[0].set_xlabel('Frequency (GHz)', fontsize = 15)
    axs[0].set_ylabel('Real permitivity', fontsize = 15)
    axs[0].set_title('Lunar Sample', fontsize = 15)

    axs[0].set_ylim(0.5,1.5)

    axs[1].plot(fl, epi, marker = 'o', linestyle = '--')
    axs[1].set_xlabel('Frequency (GHz)', fontsize = 15)
    axs[1].set_ylabel('Loss Tangent', fontsize = 15)
    #axs[1].set_title('Lunar Sample', fontsize = 15)
    #axs[1].axhline(4.5e-4)
    axs[1].set_ylim(1e-4,1)
    axs[1].set_yscale('log')

plt.tight_layout()
axs[0].grid()
axs[1].grid()
axs[0].legend(ncols = 2, loc = 'upper right')
plt.show()
print('real: ', np.round(np.array(epr), 3), np.round(np.array(eprstd), 6), np.round(np.array(epr_err), 6))
print('loss tan: ', np.round(np.array(epi), 5), np.round(np.array(elstd), 8), np.round(np.array(el_err), 8))
print('imag: ', np.round(np.array(epi)*np.array(epr), 5), np.round(np.array(epistd), 8), np.round(np.array(epi_err), 8))
#plt.plot(fl, np.array(epi)/np.array(epr), marker = 'o', linestyle = '--')
#plt.xlabel('Frequency (GHz)', fontsize = 15)
#plt.ylabel('Imag permitivity', fontsize = 15)
#plt.title('Rexolite Rod', fontsize = 15)
#plt.grid()
#plt.yscale('log')
#plt.axhline(4.5e-4)
#plt.ylim(1e-5,1)
#plt.show()
