import numpy as np

def getMur(DataFrame, reWin, hoffset = 0, holder = False, DEBUG = False):
    ## Get the resonant frequency
    f0 = DataFrame.iloc[0,2]
    f = np.array(DataFrame.iloc[:,2])
    Aeff_H = DataFrame['Aeff_Hs'][0]
    Aeff_E = DataFrame['Aeff_E'][0]
    epsrs = DataFrame['epsrs'][0]
    if(holder):
        Aeff_s = DataFrame['Aeff_Es'][0]
        epsrh = DataFrame['epsrh'][0]
    h = np.array(DataFrame['height (mm)']) - hoffset
    Vc = DataFrame['WG_length (mm)'][0]*DataFrame['WG_d (mm)'][0]*DataFrame['WG_a (mm)'][0]

    ##Calculate resonant shift and volume ratio
    yr = (f0 - f)/f0
    x = h*Aeff_H/Vc

    ## Fit data to determine nonuniform field effects
    p, V = np.polyfit(x[reWin[0]:reWin[1]], yr[reWin[0]:reWin[1]], 1, cov = True)
    m = p[0]
    b = p[1]

    if(holder):
        mur = (yr - b)/(2*x) + 1 - (epsrh - 1)*Aeff_E/Aeff_H - (epsrs - epsrh)*Aeff_s/Aeff_H
        #mur = (yr - b)/(2*x) + 1
    else:
        mur = (yr - b)/(2*x) + 1 - (2.54 - 1)*Aeff_E/Aeff_H

    if (DEBUG):
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize = (10,5))
        axs[0].plot(x, yr, linestyle = '', marker = 'o')
        axs[0].plot(x, m*x + b)

        axs[1].plot(x, mur, linestyle = '', marker = 'o')
        axs[1].set_ylim(0,2)
        #axs[1].set_yscale('log')
        plt.show()
    
    return np.mean(mur[reWin[0]:reWin[1]])