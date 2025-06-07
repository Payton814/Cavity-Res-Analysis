import numpy as np

def getMui(DataFrame, imWin, hoffset = 0, holder = False, DEBUG = False):
    ## Get the resonant frequency
    try:
        Q0 = DataFrame['Q spline'][0]
        Q = np.array(DataFrame['Q spline'])
    except:
        Q0 = DataFrame['Q'][0]
        Q = np.array(DataFrame['Q'])
    Aeff_H = DataFrame['Aeff_Hs'][0]
    Aeff_E = DataFrame['Aeff_E'][0]
    epsis = DataFrame['epsis'][0]
    if(holder):
        Aeff_s = DataFrame['Aeff_Es'][0]
        epsih = DataFrame['epsih'][0]
    h = np.array(DataFrame['height (mm)']) - hoffset
    Vc = DataFrame['WG_length (mm)'][0]*DataFrame['WG_d (mm)'][0]*DataFrame['WG_a (mm)'][0]

    ##Calculate resonant shift and volume ratio
    yi = 1/Q - 1/Q0
    x = h*Aeff_H/Vc

    ## Fit data to determine nonuniform field effects
    p, V = np.polyfit(x[imWin[0]:imWin[1]], yi[imWin[0]:imWin[1]], 1, cov = True)
    m = p[0]
    b = p[1]
    print(b)

    if(holder):
        mui = (yi - b)/(4*x) - epsih*Aeff_E/Aeff_H - epsis*Aeff_s/Aeff_H
        #mui = (yi - b)/(4*x)
    else:
        mui = (yi - b)/(4*x) - epsis*Aeff_E/Aeff_H

    if (DEBUG):
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize = (10,5))
        axs[0].plot(x, yi, linestyle = '', marker = 'o')
        axs[0].plot(x, m*x + b)

        axs[1].plot(x, mui, linestyle = '', marker = 'o')
        axs[1].set_ylim(1e-7, 1)
        axs[1].set_yscale('log')
        plt.show()
    
    return np.mean(mui[imWin[0]:imWin[1]])