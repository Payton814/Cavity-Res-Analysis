import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


peaks = [1, 3, 5, 7]
epr = []
epi = []
fl = []
r = [[45, 51], [45, 51], [45, 51], [45, 51]]
ri = [[43, 52], [43, 52], [43, 52], [43, 52]]

Cavity = "Model"
for iii in peaks:
    k = iii//2 + 0
    #df = pd.read_csv('../CROAQ/Data/Rexolite_test2_Peak' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak' + str(iii) + '_t2/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_P' + str(iii) + '_082725/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Rexolite_091025_P' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/RexoliteRod_N5230C_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/AcrylicSolid_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/AcrylicSolid_uwave_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/SolidRexolite_WR159_p' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/AcrylicSolid_p' + str(iii) + '/trial1.csv')
    df = pd.read_csv('./Example Data/CROAQ_Data/Plastics/Rexolite_101325_p' + str(iii) + '/trial1.csv')
    f = np.array(df['Frequency (GHz)'])
    Q = np.array(df['Q spline'])
    h = np.array(df['height (mm)'])
    Aeff = df['Aeff_E'][0]
    Vc = df['WG_length (mm)'][0]*df['WG_d (mm)'][0]*df['WG_a (mm)'][0]
    f0 = f[0]
    fl.append(f0)
    Q0 = Q[0]
    yr = (f0 - f)/f0
    yi = 1/Q - 1/Q0
    ## 1-3 GHz cavity lid is 0.25 in thick
    x = Aeff*(h - 0.25*25.4)/(Vc)

    ## Microwave (3-5 GHz) cavity wall is 1.675mm thick
    #x = Aeff*(h - 1)/(Vc)

    #print('percent inserted: ', (h - 0.25*25.4)[r[0]]/29.08*100, (h - 0.25*25.4)[r[1]]/29.08*100)
    p, V = np.polyfit(x[r[k][0]:r[k][1]], yr[r[k][0]:r[k][1]], 1, cov = True)

    m = p[0]
    b = p[1]
    #print(b)

    dyr = np.gradient(yr, h)
    #plt.plot(h, dyr, marker = 'o', linestyle = '--')
    #plt.show()

    #plt.plot(x, yr, marker = 'o', linestyle = '')
    #plt.plot(x, m*x + b)
    #plt.axvline(x[r[k][0]])
    #plt.axvline(x[r[k][1]])
    print("b: ", b)
   # plt.show()

    pi, Vi = np.polyfit(x[ri[k][0]:ri[k][1]], yi[ri[k][0]:ri[k][1]], 1, cov = True)

    mi = pi[0]
    bi = pi[1]

    #plt.plot(x, yi, marker = 'o', linestyle = '')
    #plt.plot(x, mi*x + bi)
    #plt.axvline(x[ri[k][0]])
    #plt.axvline(x[ri[k][1]])
    #plt.show()

    epsr = (yr - b)/(2*x)+1
    epsi = (yi - bi)/(4*x)

    if (Cavity == 'WR229'):
        if (iii == 3):
            #a = 0.992475
            a = 1
            #a = -0.00721*3.82 + 1.0114
        elif (iii == 5):
            #a = 1.0085302942
            #a = -.001807988*3.82 + 1.013277
            a = -0.001*epsr0 + 1.02
        elif (iii == 7):
            #a = 1.03939121857
            #a = 0.002710683*3.82 + 1.0322745
            a = 0.012*epsr0 + 1
        elif (iii == 9):
            #a = 1.06628986
            #a = 0.01370348*3.82 + 1.030312
            #a = 1.19
            a = 0.039*epsr0 + 0.957
    elif (Cavity == "UHF1to3"):
        if (iii == 1):
            a = 1
        elif (iii == 3):
            a = -0.008*(epsr0 - 1) + 1.014
        elif (iii == 5):
            a = 0*(epsr0 - 1) + 1.014
        elif (iii == 7):
            a = 0.025*(epsr0 - 1) + 0.995
    elif (Cavity == "Model"):
        if (iii == 1):
            a = 1
            epsr = (yr - b)/(2*x)+1
            epsi = (yi - bi)/(4*x)
        elif (iii == 3):
            #a = 0.010356 * epsr0 + 0.989667
            a = 0.011292 * epsr0 + 0.988734
            alpha = 0.011292
            beta = 1
        elif (iii == 5):
            #a = 0.021773 * epsr0 + 0.978257
            a = 0.023582 * epsr0 + 0.976451
            alpha = 0.023582
            beta = 1
        elif (iii == 7):
            #a = 0.049707 * epsr0 + 0.950261
            a = 0.056063 * epsr0 + 0.943929
            #a = 0.042286 * epsr0 + 0.957639
            #a = 0.03 * epsr0 + 0.97
            alpha = 0.056063
            beta = 1
    else:
        a = 1

    #a = 1

    #a_s = np.cos(np.pi*9*iii/(2*6.5*25.4))/np.cos(np.pi*9*iii*np.sqrt(2.54)/(2*6.5*25.4))
    #a_h = 1 - np.pi*9**2*f0/(4*300*20.193)
    #a = a_s*a_h

    epsr = (yr - b)/(2*x*a)+1
    epsi = (yi - bi)/(4*x*a)

    if (iii != 1):
        eps = (yr - b)/(2*alpha*x)
        rt = eps**2 - 4*beta/alpha
        print(eps)
        print(4*beta/alpha)

        epsr = (1 + eps + np.sqrt(rt))/2
        rt = (eps/2)**2 - 4*beta/alpha
        epsi = (eps/2 + np.sqrt(rt))/2

    if (iii == 1):
        a = 1
        epsr0 = np.mean(epsr[r[k][0]:r[k][1]])
    if (Cavity == 'WR229' and iii == 3):
        epsr0 = np.mean(epsr[r[k][0]:r[k][1]])
    epr.append(np.mean(epsr[r[k][0]:r[k][1]]))
    epi.append(np.mean(epsi[ri[k][0]:ri[k][1]])/np.mean(epsr[r[k][0]:r[k][1]]))
    #print('Real Permittivity', np.mean(epsr[r[k][0]:r[k][1]]))
    #print('loss tangent', np.mean(epsi[ri[k][0]:ri[k][1]])/np.mean(epsr[r[k][0]:r[k][1]]))

    #plt.plot(x, epsr, marker = 'o', linestyle = '')
    #plt.xlabel('Frequency (GHz)', fontsize = 15)
    #plt.ylabel('Real Permittivity', fontsize = 15)
    #plt.title('$TE_{101}$, Rexolite Rod', fontsize = 15)
    #plt.axhline(2.54)
    #plt.axvline(x[r[k][0]])
    #plt.axvline(x[r[k][1]])
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
    k = k + 1

fig, axs = plt.subplots(2, 1)
axs[0].plot(fl, epr, marker = 'o', linestyle = '--')
axs[0].set_xlabel('Frequency (GHz)', fontsize = 15)
axs[0].set_ylabel('Real permitivity', fontsize = 15)
axs[0].set_title('Lunar Sample', fontsize = 15)
axs[0].grid()
axs[0].axhline(2.53)
axs[0].set_ylim(1,5)

#print(epi)
axs[1].plot(fl, epi, marker = 'o', linestyle = '--')
axs[1].set_xlabel('Frequency (GHz)', fontsize = 15)
axs[1].set_ylabel('Loss Tangent', fontsize = 15)
#axs[1].set_title('Lunar Sample', fontsize = 15)
axs[1].grid()
axs[1].axhline(4.5e-4)
axs[1].set_ylim(1e-4,1)
axs[1].set_yscale('log')

plt.tight_layout()
plt.show()

print(np.array(epr))
print(np.array(epi))
print(np.array(epi)*np.array(epr))
