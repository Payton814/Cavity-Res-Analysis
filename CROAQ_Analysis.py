import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


peaks = [3]
epr = []
epi = []
fl = []
for iii in peaks:
    df = pd.read_csv('../CROAQ/Data/Rexolite_test2_Peak' + str(iii) + '/trial1.csv')
    df = pd.read_csv('../CROAQ/Data/Rexolite_Cav2_Peak' + str(iii) + '/trial1.csv')
    f = np.array(df['Frequency (GHz)'])
    #Q = np.array(df['Q spline'])
    h = np.array(df['height (mm)'])
    Aeff = df['Aeff_E'][0]
    Vc = df['WG_length (mm)'][0]*df['WG_d (mm)'][0]*df['WG_a (mm)'][0]

    f0 = f[0]
    fl.append(f0)
    #Q0 = Q[0]
    yr = (f0 - f)/f0
    #yi = 1/Q - 1/Q0
    x = Aeff*(h - 0.25*25.4)/(Vc)

    p, V = np.polyfit(x[45:55], yr[45:55], 1, cov = True)

    m = p[0]
    b = p[1]
    print(b)

    plt.plot(x, yr, marker = 'o', linestyle = '')
    plt.plot(x, m*x + b)
    print("b: ", b)
    plt.show()

    #pi, Vi = np.polyfit(x[45:55], yi[45:55], 1, cov = True)

    #mi = pi[0]
    #bi = pi[1]

    #plt.plot(x, yi, marker = 'o', linestyle = '')
    #plt.plot(x, mi*x + bi)
    #plt.show()

    epsr = (yr - b)/(2*x)+1
    #epsi = (yi - bi)/(4*x)

    epr.append(np.mean(epsr[45:55]))
    #epi.append(np.mean(epsi[45:55]))
    print('Real Permittivity', np.mean(epsr[45:55]))
    #print('loss tangent', np.mean(epsi[45:55])/np.mean(epsr[45:55]))

    plt.plot(x, epsr, marker = 'o', linestyle = '')
    plt.xlabel('Frequency (GHz)', fontsize = 15)
    plt.ylabel('Real Permittivity', fontsize = 15)
    plt.title('$TE_{101}$, Rexolite Rod', fontsize = 15)
    plt.axhline(2.54)
    plt.grid()
    plt.ylim(1, 3)

    plt.show()

    #plt.plot(x, np.array(epsi)/np.array(epsr), marker = 'o', linestyle = '')
    #plt.xlabel('Frequency (GHz)', fontsize = 15)
    #plt.ylabel('loss tangent', fontsize = 15)
    #plt.title('Rexolite Rod', fontsize = 15)
    #plt.axhline(4.5e-4)
    #plt.grid()
    #plt.yscale('log')

    #plt.show()

plt.plot(fl, epr, marker = 'o', linestyle = '--')
plt.xlabel('Frequency (GHz)', fontsize = 15)
plt.ylabel('Real permitivity', fontsize = 15)
plt.title('Rexolite Rod', fontsize = 15)
plt.grid()
plt.axhline(2.53)
plt.ylim(1,3)
plt.show()

'''plt.plot(fl, np.array(epi)/np.array(epr), marker = 'o', linestyle = '--')
plt.xlabel('Frequency (GHz)', fontsize = 15)
plt.ylabel('Imag permitivity', fontsize = 15)
plt.title('Rexolite Rod', fontsize = 15)
plt.grid()
plt.yscale('log')
plt.axhline(4.5e-4)
plt.ylim(1e-5,1)
plt.show()'''
