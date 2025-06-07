import numpy as np

def getEpsr(DataFrame, reWin, hoffset = 0, holder = False, DEBUG = False):
    ## Get the resonant frequency
    f0 = DataFrame.iloc[0,2]
    f = np.array(DataFrame.iloc[:,2])
    Aeff_E = DataFrame['Aeff_E'][0]
    if(holder):
        Aeff_s = DataFrame['Aeff_Es'][0]
        epsrh = DataFrame['epsrh'][0]
    h = np.array(DataFrame['height (mm)']) - hoffset
    Vc = DataFrame['WG_length (mm)'][0]*DataFrame['WG_d (mm)'][0]*DataFrame['WG_a (mm)'][0]
    #print(Vc)
    #Vc = 123.647244*25.4*25.4*25.4
    ##Calculate resonant shift and volume ratio
    yr = (f0 - f)/f0
    x = h*Aeff_E/Vc

    ## Fit data to determine nonuniform field effects
    p, V = np.polyfit(x[reWin[0]:reWin[1]], yr[reWin[0]:reWin[1]], 1, cov = True)
    m = p[0]
    b = p[1]

    if (holder):
        epsr = ((yr - b)/(2*x)+1-epsrh)*Aeff_E/Aeff_s + epsrh
    else:
        epsr = (yr - b)/(2*x) + 1
        print(np.std(epsr[reWin[0]:reWin[1]]))

    if (DEBUG):
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize = (10,5))
        axs[0].plot(x, yr, linestyle = '', marker = 'o')
        axs[0].axvline(x[reWin[0]], linestyle = '--')
        axs[0].axvline(x[reWin[1]], linestyle = '--')
        axs[0].plot(x, m*x + b)

        axs[1].plot(x, epsr, linestyle = '', marker = 'o')
        axs[1].axvline(x[reWin[0]], linestyle = '--')
        axs[1].axvline(x[reWin[1]], linestyle = '--')
        axs[1].set_ylim(1, 5)
        #axs[1].set_yscale('log')
        plt.show()
    
    return np.mean(epsr[reWin[0]:reWin[1]])



