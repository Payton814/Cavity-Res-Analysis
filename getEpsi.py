import numpy as np

def getEpsi(DataFrame, imWin, hoffset = 0, holder = False, DEBUG = False):
    ## Get the resonant frequency
    try:
        Q0 = DataFrame['Q spline'][0]
        Q = np.array(DataFrame['Q spline'])
    except:
        Q0 = DataFrame['Q'][0]
        Q = np.array(DataFrame['Q'])
    Aeff_E = DataFrame['Aeff_E'][0]
    if(holder):
        Aeff_s = DataFrame['Aeff_Es'][0]
        epsih = DataFrame['epsih'][0]
    h = np.array(DataFrame['height (mm)']) - hoffset
    Vc = DataFrame['WG_length (mm)'][0]*DataFrame['WG_d (mm)'][0]*DataFrame['WG_a (mm)'][0]

    ##Calculate resonant shift and volume ratio
    yi = 1/Q - 1/Q0
    x = h*Aeff_E/Vc

    ## Fit data to determine nonuniform field effects
    p, V = np.polyfit(x[imWin[0]:imWin[1]], yi[imWin[0]:imWin[1]], 1, cov = True)
    m = p[0]
    b = p[1]

    if(holder):
        epsi = ((yi - b)/(4*x) - epsih)*Aeff_E/Aeff_s + epsih
    epsi = (yi - b)/(4*x)

    if (DEBUG):
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize = (10,5))
        axs[0].plot(x, yi, linestyle = '', marker = 'o')
        axs[0].plot(x, m*x + b)

        axs[1].plot(x, epsi, linestyle = '', marker = 'o')
        axs[1].set_ylim(1e-5, 1)
        axs[1].set_yscale('log')
        plt.show()
    
    return np.mean(epsi[imWin[0]:imWin[1]])