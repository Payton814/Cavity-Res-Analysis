import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


peaks = [1, 3, 5, 7]
epr = []
epi = []
fl = []
r = [50, 60]
ri = [50, 60]
for iii in peaks:
    #df = pd.read_csv('../CROAQ/Data/Rexolite_test2_Peak' + str(iii) + '/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak' + str(iii) + '_t2/trial1.csv')
    #df = pd.read_csv('../CROAQ/Data/Rexolite_P' + str(iii) + '_082725/trial1.csv')
    #df = pd.read_csv('./Example Data/CROAQ_Data/Rexolite_091025_P' + str(iii) + '/trial1.csv')
    df = pd.read_csv('./Example Data/CROAQ_Data/RexoliteRod_N5230C_p' + str(iii) + '/trial1.csv')
    f = np.array(df['Frequency (GHz)'])
    Q = np.array(df['Q raw'])
    h = np.array(df['height (mm)'])
    Aeff = df['Aeff_E'][0]
    Vc = df['WG_length (mm)'][0]*df['WG_d (mm)'][0]*df['WG_a (mm)'][0]

    f0 = f[0]
    fl.append(f0)
    Q0 = Q[0]
    yr = (f0 - f)/f0
    yi = 1/Q - 1/Q0
    x = Aeff*(h - 0.25*25.4)/(Vc)

    p, V = np.polyfit(x[r[0]:r[1]], yr[r[0]:r[1]], 1, cov = True)

    m = p[0]
    b = p[1]
    print(b)

    dyr = np.gradient(yr, h)
    #plt.plot(h, dyr, marker = 'o', linestyle = '--')
    #plt.show()

    #plt.plot(x, yr, marker = 'o', linestyle = '')
    #plt.plot(x, m*x + b)
    #plt.axvline(x[r[0]])
    #plt.axvline(x[r[1]])
    print("b: ", b)
    #plt.show()

    pi, Vi = np.polyfit(x[ri[0]:ri[1]], yi[ri[0]:ri[1]], 1, cov = True)

    mi = pi[0]
    bi = pi[1]

    #plt.plot(x, yi, marker = 'o', linestyle = '')
    #plt.plot(x, mi*x + bi)
    #plt.axvline(x[ri[0]])
    #plt.axvline(x[ri[1]])
    #plt.show()

    epsr = (yr - b)/(2*x)+1
    epsi = (yi - bi)/(4*x)

    if (iii == 1):
        a = 1
    elif (iii == 3):
        a = -0.008*(epsr0 - 1) + 1.014
    elif (iii == 5):
        a = 0*(epsr0 - 1) + 1.014
    elif (iii == 7):
        a = 0.025*(epsr0 - 1) + 0.995

    #a = 1
    epsr = (yr - b)/(2*x*a)+1
    epsi = (yi - bi)/(4*x*a)
    if (iii == 1):
        a = 1
        epsr0 = np.mean(epsr[r[0]:r[1]])
    epr.append(np.mean(epsr[r[0]:r[1]]))
    epi.append(np.mean(epsi[ri[0]:ri[1]]))
    print('Real Permittivity', np.mean(epsr[r[0]:r[1]]))
    print('loss tangent', np.mean(epsi[ri[0]:ri[1]])/np.mean(epsr[r[0]:r[1]]))

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

plt.plot(fl, epr, marker = 'o', linestyle = '--')
plt.xlabel('Frequency (GHz)', fontsize = 15)
plt.ylabel('Real permitivity', fontsize = 15)
plt.title('Rexolite Rod', fontsize = 15)
plt.grid()
plt.axhline(2.53)
plt.ylim(1,3)
plt.show()

#plt.plot(fl, np.array(epi)/np.array(epr), marker = 'o', linestyle = '--')
#plt.xlabel('Frequency (GHz)', fontsize = 15)
#plt.ylabel('Imag permitivity', fontsize = 15)
#plt.title('Rexolite Rod', fontsize = 15)
#plt.grid()
#plt.yscale('log')
#plt.axhline(4.5e-4)
#plt.ylim(1e-5,1)
#plt.show()
