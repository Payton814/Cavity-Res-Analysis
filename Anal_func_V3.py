###################################################################################################################
#
#   Code: Sample_Analysis_func.py
#
#   Programmer: Payton Linton (linton.93@osu.edu)
#
#   Purpose: This function uses the data for the frequency shift. It assumes a solid rod
#            therefore it does not apply any frequency corrections.
#            The effective volume correction should be calculated seperately and is an input.
#            If no volume correction is given then the real volume will be used instead.
#
#   Inputs: dataframe: this is the csv file that has all the data. This function makes references specific columns in the data
#                      so make sure the csv is formatted correctly
#
#           Sample_Aeff: This is the effective area for the sample which comes from integrating the power over the sample cross section
#
#           plot: If this is made true then this function will plot the frequency and Q shifts, as well as the real permittivity and
#                 loss tangent.
#
#           holder_params: If the sample is inside of a holder then this is the effective area of the holder
#
#           real_range: The real permittivity relies on fitting a slope to the frequency shift, this chooses
#                       the range that a linear fit will be found for the frequency shift, as well as the permittivity
#                       values that are averaged
#
#
#
###############################################################################################################


def cavity_analysis(dataframe, Sample_Aeff = None, plot = False, holder_params = None,
                     real_range = [5, 18], imag_range = [5, 18], getDensity = False, Solid = False,
                     Powder = False, Corr = False, epsr3 = None):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    ## The following if statement is for samples that were placed inside of a holder. For the small microwave cavity
    ## the holder was always rexolite. So currently the rexolite measurement is what is hardcoded. This would need to 
    ## change for a holder of a different material.
    if (Solid == False):
        df = pd.read_csv('./Measurements/Rexolite/RexoliteRod1_PostBake_Peak' + str(int(dataframe['Mode'][0])) + '.csv')
        xdataSolid, real_dataSolid, imag_dataSolid, loss_tanSolid = cavity_analysis(df, df['Aeff'][0], 
                                                                        plot = False, 
                                                                        real_range = [10,19], imag_range = [10,17], Solid = True)

    ## If the sample is placed inside of a holder then the effective area of the holder must be given
    if (holder_params != None):
        Aeff = holder_params[0]
    elif (Sample_Aeff != None):
        Aeff = Sample_Aeff
    else:
        Aeff = None
   
    f3 = dataframe['Frequency (GHz)']
    sample_name = dataframe['Sample'][0]
    sample_height = np.array(dataframe['height (mm)'])
    Q = np.array(dataframe['Q'])


    ID = dataframe['ID (mm)'][0]
    if(getDensity == True):
        SampleV = np.pi*(ID)**2*dataframe['sample height (mm)'][0]/4
        Density = dataframe['Sample Weight (mg)'][0]/SampleV
    OD = dataframe['OD base (mm)'][0]
    if (Solid == False and Powder == False):
        Sample_OD = dataframe['Sample OD (mm)'][0]
    V_cavity = dataframe['WG_length (mm)'][0]*dataframe['WG_d (mm)'][0]*dataframe['WG_a (mm)'][0]

    f0 = f3[0] ## The unperturbed resonance is assumed to be the first value in the data since the sample had not yet been inserted
    y_r = []
    y_i = []
    x = []
    ph = [] ## this is the insertion percent
    WG_wall_thickness = dataframe['WG_wall_thickness (mm)'][0]
    for i in range(len(f3)):
        ## Need to make sure we are using the correct sample inserted height
        ## 0 is defined as the top of the cavity, but the sample doesnt enter the cavity
        ## until it passes through the wall thickness.
        ## The base of the vial also enters first, adding another offset to when the sample
        ## enters
        if (sample_height[i] - WG_wall_thickness > 0):
            if (Aeff != None):
                Vs = (Aeff)*(sample_height[i] - WG_wall_thickness)
            else:
                Vs = np.pi*OD**2*(sample_height[i] - WG_wall_thickness)/4
        ## If the sample has not entered the cavity set the inserted valume to 0
        ## This will cause a divide by 0 error for loss tangent, but can be ignored
        else:
            Vs = 0
        ph.append((sample_height[i] - WG_wall_thickness)/29.08)
        y_r.append((f0 - f3[i])/f0)
        y_i.append((f3[i]/f0)*(1/Q[i] - 1/Q[0]))
        x.append(Vs/V_cavity)

    #print((sample_height[real_range[0]+1] - WG_wall_thickness)/sample_height[-1])
    #print((sample_height[real_range[1]] - WG_wall_thickness)/sample_height[-1])
    p, V = np.polyfit(x[real_range[0]:real_range[1]], y_r[real_range[0]:real_range[1]], 1, cov = True)
    pi, Vi = np.polyfit(x[imag_range[0]:imag_range[1]], y_i[imag_range[0]:imag_range[1]], 1, cov = True)

    m = p[0]
    b = p[1]
    #b=0
    
    mi = pi[0]
    bi = pi[1]
    #bi=0

    eps_r = []
    Del_eps_r_list = []
    eps_i = []

    if (holder_params == None):
        try:
            ## Accounts for the increased E-field strength in samples
            a0 = dataframe['EperE0'][0] - 1
            a = 1 + (a0)/2
            #a = dataframe['EperE0'][0]
        except:
            #print("No magnitude correction factor found")
            a = 1
            if (epsr3 != None):
                if (dataframe['Mode'][0] == 3):
                        ## For the lowest mode the trend is so close to constant 1
                        ## I assume the lowest mode is the 'true' permittivity value
                        a = 1
                elif (dataframe['Mode'][0] == 5):
                        a = -0.001*epsr3 + 1.02
                elif (dataframe['Mode'][0] == 7):
                        a = 0.012*epsr3 + 1
                elif (dataframe['Mode'][0] == 9):
                        a = 0.039*epsr3 + 0.957
        if (Corr == False):
            a = 1 ## If this is uncommented, then the correction will not be applied
    else:
        ## If holder params were given then this correction is handled later.
        a = 1
    for ii in range(len(x)):
        ## Determine the permittivity value based on Orloff's technique
        eps_rm = ((y_r[ii] - b)/(2*x[ii]*a) + 1)
        eps_r.append(eps_rm)
        eps_i.append((y_i[ii] - bi)/(4*x[ii]*a))

    Del_eps_r = np.sqrt(np.dot(Del_eps_r_list[real_range[0]: real_range[1]], Del_eps_r_list[real_range[0]: real_range[1]]))/(real_range[1] - real_range[0])

    if (real_range == [17,24]):
        print(np.mean(eps_r[real_range[0]:real_range[1]]))

    loss_tan = []
    for i in range(len(eps_r)):
        loss_tan.append(eps_i[i]/eps_r[i])


    ## The next block applies the correction to how you deal with the holder for the sample
    if (holder_params != None):
        for i in range(len(eps_r)):
            ##eps_r[i] = (eps_r[i] - 2.54)*Aeff/Sample_Aeff + 2.54 + 2*1.54*1.35491/Sample_Aeff
            if (epsr3 != None):
                if (dataframe['Mode'][0] == 3):
                    #a = 0.992475
                    a = 1
                    if (dataframe['Sample Name'][0] == 'K5'):
                        a = 0.98807
                    #a = -0.00721*3.82 + 1.0114
                elif (dataframe['Mode'][0] == 5):
                    
                    #a = 1.0085302942
                    #a = -.001807988*3.82 + 1.013277
                    #a = -0.001*epsr3 + 1.02
                    a = 0.0027*epsr3 + 1.00058
                    if (dataframe['Sample Name'][0] == 'K5'):
                        a = 1.01455
                elif (dataframe['Mode'][0] == 7):
                    #a = 1.03939121857
                    #a = 0.002710683*3.82 + 1.0322745
                    #a = 0.012*epsr3 + 1
                    a = 0.0122*epsr3 + 0.9991
                    if (dataframe['Sample Name'][0] == 'K5'):
                        a = 1.0474
                elif (dataframe['Mode'][0] == 9):
                    #a = 1.06628986
                    #a = 0.01370348*3.82 + 1.030312
                    #a = 1.19
                    #a = 0.039*epsr3 + 0.957
                    a = 0.0273*epsr3 + 0.9901
                    if (dataframe['Sample Name'][0] == 'K5'):
                        a = 1.104277
            else:
                try:
                    a = dataframe['a'][0]
                except:
                    a = 1
            if (Corr == False):
                a = 1
            #print('Using a = ' + str(a) + ' for holder step')
            eps_r[i] = ((eps_r[i] - holder_params[1])*Aeff/Sample_Aeff + (holder_params[1]-1))/a + 1
            #eps_r[i] = (eps_r[i] - holder_params[1])*Aeff/Sample_Aeff + holder_params[1]
            eps_i[i] = ((eps_i[i] - holder_params[2])*Aeff/Sample_Aeff + holder_params[2])/a
    if (holder_params == None and Sample_Aeff == None and Solid == False):
        for i in range(len(eps_r)):
            ##eps_r[i] = (eps_r[i] - 2.54)*Aeff/Sample_Aeff + 2.54 + 2*1.54*1.35491/Sample_Aeff
            eps_r[i] = (eps_r[i] - holder_params[1])*OD**2/Sample_OD**2 + holder_params[1]
            eps_i[i] = (eps_i[i] - holder_params[2])*OD**2/Sample_OD**2 + holder_params[2]

    #### Calculate the RMSE error
    epsr_RMSE = np.sqrt(np.square(np.array(eps_r[real_range[0]:real_range[1]]) - np.mean(eps_r[real_range[0]:real_range[1]])).mean()/(real_range[1] - real_range[0]))
    epsi_RMSE = np.sqrt(np.square(np.array(eps_i[imag_range[0]:imag_range[1]]) - np.mean(eps_i[imag_range[0]:imag_range[1]])).mean()/[imag_range[1] - imag_range[0]])

    loss_tan = []
    for i in range(len(eps_r)):
        loss_tan.append(eps_i[i]/eps_r[i])

    loss_tan_RMSE = np.sqrt(np.square(np.array(loss_tan[imag_range[0]:imag_range[1]]) - np.mean(loss_tan[imag_range[0]:imag_range[1]])).mean()/(imag_range[1] - imag_range[0]))



    ######### TESTING A CORRECTION IDEA ##############


    ## This block plots a 2x2 grid of plots showing the frequency shift as a function of depth,
    ## Q shift as a function of depth, and apparent permittivity and loss tangent as a function of depth. 
    if (plot == True):
        fig, axs = plt.subplots(2, 2)
        #axs[0,0].plot(step3, f3, marker = 'o')
        #axs[0,0].set_title('Frequency shift vs. step number')
        #axs[0,0].set_xlabel('Step number')
        #axs[0,0].set_ylabel('$(f_s - f_0)/f_0$')
        axs[0,0].plot(x, y_r, marker = 'o', label = 'Data')
        axs[0,0].plot(x, m*np.array(x) + b, label = 'Best Fit Line')
        axs[0,0].legend()
        axs[0,0].set_title('frequency shift vs fractional volume')
        axs[0,0].set_xlabel('fractional volume')
        axs[0,0].set_ylabel('$(f_s - f_0)/f_0$')
        axs[0,0].grid()
        axs[0,0].axvline(x[real_range[0]], color = 'k', linestyle = '--')
        axs[0,0].axvline(x[real_range[1]], color = 'k', linestyle = '--')

        axs[0,1].plot(x, y_i, marker = 'o', label = 'Data')
        axs[0,1].plot(x, mi*np.array(x) + bi, label = 'Best Fit Line')
        axs[0,1].legend()
        axs[0,1].set_title('Q shift vs fractional volume')
        axs[0,1].set_xlabel('fractional volume')
        axs[0,1].set_ylabel('$1/Q_s - 1/Q_0$')
        axs[0,1].grid()
        axs[0,1].axvline(x[imag_range[0]], color = 'k', linestyle = '--')
        axs[0,1].axvline(x[imag_range[1]], color = 'k', linestyle = '--')

        axs[1,0].axhline(np.mean(eps_r[real_range[0]:real_range[1]]), color = 'r', label = '$\epsilon_r$ = ' + str(np.round(np.mean(eps_r[real_range[0]:real_range[1]]), 2)))
        axs[1,0].axvline(sample_height[real_range[0]], color = 'k', linestyle = '--')
        axs[1,0].axvline(sample_height[real_range[1]], color = 'k', linestyle = '--')
        axs[1,0].plot(sample_height, eps_r, marker = 'o')
        axs[1,0].set_title('Real Permittivity vs sample depth')
        axs[1,0].set_ylabel('$\epsilon_r$')
        axs[1,0].set_xlabel('sample depth (mm)')
        axs[1,0].grid()
        axs[1,0].legend()

        axs[1,1].axhline(np.mean(eps_i[imag_range[0]:imag_range[1]]), color = 'r', label = '$tan(\delta)$ = ' + str(np.round(np.mean(eps_i[imag_range[0]:imag_range[1]]), 5)))
        axs[1,1].axvline(sample_height[imag_range[0]], color = 'k', linestyle = '--')
        axs[1,1].axvline(sample_height[imag_range[1]], color = 'k', linestyle = '--')
        axs[1,1].plot(sample_height, eps_i, marker = 'o')
        axs[1,1].set_title('Imaginary Permittivity vs sample depth')
        axs[1,1].set_ylabel('$\epsilon_i$, imaginary permittivity')
        axs[1,1].set_xlabel('sample depth (mm)')
        axs[1,1].grid()
        axs[1,1].legend()
        fig.suptitle( str(sample_name) +' using TE10' + str(int(dataframe['Mode'][0])), fontweight = 'bold')
        plt.tight_layout()
        plt.show()

    if(getDensity == True):
        return (x, sample_height), (y_r, eps_r, np.mean(eps_r[real_range[0]:real_range[1]]), Del_eps_r, epsr_RMSE, b, m, (ph[real_range[0]+1], ph[real_range[1]])), (y_i, eps_i, np.mean(eps_i[imag_range[0]:imag_range[1]]), epsi_RMSE, bi, mi), (loss_tan, np.mean(loss_tan[imag_range[0]:imag_range[1]]), loss_tan_RMSE), Density
    else:
        return (x, sample_height), (y_r, eps_r, np.mean(eps_r[real_range[0]:real_range[1]]), Del_eps_r, epsr_RMSE, b, m, (ph[real_range[0]+1], ph[real_range[1]])), (y_i, eps_i, np.mean(eps_i[imag_range[0]:imag_range[1]]), epsi_RMSE, bi, mi), (loss_tan, np.mean(loss_tan[imag_range[0]:imag_range[1]]), loss_tan_RMSE)