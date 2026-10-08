from EigenUtils import find_linear_region
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sample_definitions import SAMPLES

from cavity_sample_field_neighbor_normalized import (
    Cavity,
    Sample,
    solve_loaded_cavity,
    loaded_field,
    loaded_field_magnitude,
    field_inside_sample,
    center_field,
    field_enhancement_over_empty,
    empty_te10n_field,
    neighboring_antinode_fields,
)

f = [1.04, 1.39, 1.91, 2.5]
LS1epsr = [3.98, 3.97, 3.97, 3.97]
Rexolite = np.array([2.53, 2.54, 2.54, 2.54])
Teflon = np.array([2.03, 2.04, 2.04, 2.05])
KelF = np.array([2.35, 2.35, 2.35, 2.35])
Acrylic = np.array([2.58, 2.58, 2.58, 2.59])
Acetal = np.array([2.93, 2.93, 2.94, 2.92])
PET = np.array([3.08, 3.08, 3.10, 3.09])
PVDF = np.array([2.83, 2.80, 2.77, 2.76])


'''SAMPLE = "Teflon1"
SAMPLE_RANGE = [1.9, 2.1]
sample = Sample(radius_mm=6.79/2,
        eps_r=2.04,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''

'''SAMPLE = "KelF1"
SAMPLE_RANGE = [2.2, 2.4]
sample = Sample(radius_mm=6.85/2,
        eps_r=2.35,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''

'''SAMPLE = "Rexolite"
SAMPLE_RANGE = [2.4, 2.6]
sample = Sample(radius_mm=6.87/2,
        eps_r=2.54,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''
#SAMPLE_iRANGE = [1e-4, 2e-2]

'''SAMPLE = "Acrylic"
SAMPLE_RANGE = [2.5, 2.8]
sample = Sample(radius_mm=6.79/2,
        eps_r=2.56,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''
#SAMPLE_iRANGE = [1e-3, 2e-1]

'''SAMPLE = "Acetal1"
SAMPLE_RANGE = [2.8, 3.1]
sample = Sample(radius_mm=6.85/2,
        eps_r=2.9,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''
#SAMPLE_iRANGE = [1e-3, 2e-1]

'''SAMPLE = "PET"
SAMPLE_RANGE = [2.9, 3.2]
#LunarSample1 Rs= 7.06/2, Rv = 9.09/2
sample = Sample(radius_mm=6.8/2,
        eps_r=3.08,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''


'''SAMPLE = "PVDF1"
SAMPLE_RANGE = [2.7, 2.9]
sample = Sample(radius_mm=6.87/2,
        eps_r=2.83,
        vial_outer_radius_mm=9.08/2,
        eps_vial=2.54,)'''
#SAMPLE_iRANGE = [1e-3, 2e-1]

#SAMPLE = "LunarSample1"
#SAMPLE_RANGE = [3.8, 4.1]
#SAMPLE_iRANGE = [1e-3, 2e-1]

samples_to_run = [
    "Teflon1",
    "Rexolite",
    "KelF1",
    "Acrylic",
    "Acetal1",
    "PVDF1",
    "PET"
]

#samples_to_run = [
#    "PVDF1"
#]



for SAMPLE in samples_to_run:

    config = SAMPLES[SAMPLE]

    # --------------------------------------------------------------
    # Permittivity range associated with this sample
    # --------------------------------------------------------------

    SAMPLE_RANGE = config["sample_range"]


    # --------------------------------------------------------------
    # Build the Sample object
    # --------------------------------------------------------------

    sample = Sample(
        radius_mm=config["radius_mm"],
        eps_r=config["eps_r"],
        vial_outer_radius_mm=config["vial_outer_radius_mm"],
        eps_vial=config["eps_vial"],
    )

    epr = []
    print("")
    print(SAMPLE)
    for n in [1, 3, 5, 7]:
        k = n//2

        print("Mode: ", n)

        if (n==1):
            r = [40,57]
            ri = [58, 68]
            iR2 = 0.85
            #R2 = 0.999535 ## 99.9535% is 3.5 sigma

            ## For LunarSample1 R2 seems to be good for 0.999975
            R2 = 0.999923 ## 99.9535% is 3.5 sigma
            startterm = 45
        elif (n==3):
            r = [50,64]
            ri = [51, 63]
            iR2 = 0.85
            #R2 = 0.9995 ## 99.9937% is 4.0 sigma
            ## For LunarSample1 R2 seems to be good for 0.999937
            R2 = 0.99977 ## 99.9535% is 3.5 sigma
            startterm = 45
        elif (n==5):
            r = [53,63]
            ri = [50, 55]
            iR2 = 0.9
            #R2 = 0.9995    ## 99.9937% is 4.0 sigma
            ## For LunarSample1 R2 seems to be good for 0.999975
            R2 = 0.99985 ## 99.9535% is 3.5 sigma
            startterm = 48
        else:
            r = [55,66]
            ri = [55, 68]
            iR2 = 0.9
            #R2 = 0.99995   ## 99.9937% is 4.0 sigma
            ## For LunarSample1 R2 seems to be good for 0.99998
            R2 = 0.9997 ## 99.9535% is 3.5 sigma
            startterm = 48
        SAMPLE_FILE = "./WR229/" + SAMPLE + "/" + SAMPLE + "Rod1_PostBake_Peak" + str(n) + ".csv"
        SAMPLE_FILE = "./Example Data/CROAQ_Data/" + SAMPLE + "Solid_uwave_p" + str(n) + "/trial1.csv"
        SAMPLE_FILE = "./Example Data/CROAQ_Data/Plastics/" + SAMPLE + "Solid_p" + str(n) + "/trial1.csv"
        SAMPLE_FILE = "./Example Data/CROAQ_Data/Plastics/" + SAMPLE + "_RexVial_p" + str(n) + "/trial1.csv"
        #SAMPLE_FILE = "./Example Data/CROAQ_Data/" + SAMPLE  + "_p" + str(n) + "/trial1.csv"
        df = pd.read_csv(SAMPLE_FILE)
        a = df['WG_d (mm)'].to_numpy()[0]
        ell = df['WG_length (mm)'].to_numpy()[0]
        b = df['WG_a (mm)'].to_numpy()[0]
        Aeff = df['Aeff_E'].to_numpy()[0]
        #Rs = df['ID (mm)'].to_numpy()[0]/2
        #Rv = df['OD base (mm)'].to_numpy()[0]/2


        #LunarSample1 Rs= 7.06/2, Rv = 9.09/2
        #Rs = 6.78/2
        #Rv = 9.08/2
        #Rs = 4.5

        cavity = Cavity(
            a_mm=a,
            ell_mm=ell,
        )

        Vc = a*ell*b

        yr = (df['Frequency (GHz)'].to_numpy()[0] - df['Frequency (GHz)'].to_numpy())/df['Frequency (GHz)'].to_numpy()[0]
        yi = 1/df['Q spline'].to_numpy() - 1/df['Q spline'].to_numpy()[0]
        #h = df['height (mm)'].to_numpy() - df['WG_wall_thickness (mm)'].to_numpy()[0]
        h = df['height (mm)'].to_numpy() - 0.25*25.4



        p, V = np.polyfit(h[r[0]:r[1]], yr[r[0]:r[1]], 1, cov=True)
        #pi, Vi = np.polyfit(h[ri[0]:ri[1]], yi[ri[0]:ri[1]], 1, cov=True)

        try:
            real_result = find_linear_region(h, yr, 5, R2, startterm)
        except:
            print("Given R2 threshold could not be met...")
            print("Reducing threshold to 0.999535")
            R2 = 0.999535
            real_result = find_linear_region(h, yr, 5, R2, startterm)

        try:
            imag_result = find_linear_region(h, yi, 7, iR2, startterm)
        except:
            imag_result = find_linear_region(h, yi, 7, 0.99, startterm)



        print(real_result["indices"], imag_result["indices"], R2)

        #plt.plot(h, yr, marker = 'o')
        #plt.axvline(h[r[1]])
        #plt.axvline(h[r[0]])
        #plt.axvline(h[real_result["start_index"]], linestyle = '--', color = 'k')
        #plt.axvline(h[real_result["end_index"]],  linestyle = '--', color = 'k')
        #plt.axvline(h[30], color = 'k', linestyle = '--')
        #plt.show()

        #plt.plot(h, yi, marker = 'o')
        #plt.axvline(h[imag_result["start_index"]])
        #plt.axvline(h[imag_result["end_index"]])
        #plt.axvline(h[30], color = 'k', linestyle = '--')
        #plt.show()


        #mr = p[0]
        br = p[1]
        br = real_result["intercept"]
        bi = imag_result["intercept"]

        #print(bi)

        #mi = pi[0]
        #bi = pi[1]
        #print(b)

        M1 = Vc*(yr - br)/(2*h)
        #M1 = np.mean(M1[r[0]:r[1]])
        M1 = np.mean(M1[real_result["start_index"]:real_result["end_index"]])

        Q1 = Vc*(yi - bi)/(4*h)
        #Q1 = np.mean(Q1[ri[0]:ri[1]])
        Q1 = np.mean(Q1[imag_result["start_index"]:imag_result["end_index"]])

        old_eps = Vc*(yr - br)/(2*h*Aeff) + 1
        #print("Old epsr: ", np.mean(old_eps[r[0]:r[1]]))

        #epsr = 2.54
        Z = []
        Z2 = []

        Av = []
        As = []
        AIDv = []

        epsrs = np.linspace(SAMPLE_RANGE[0], SAMPLE_RANGE[1], 31)
        for epsr in epsrs:

            sample.eps_r = epsr
            #sample = Sample(
            #radius_mm=Rs,
            #eps_r=epsr,
            #vial_outer_radius_mm=Rv,
            #eps_vial=2.54,)

            solution = solve_loaded_cavity(
            cavity=cavity,
            sample=sample,
            n=n,
            m_max=21,
            p_max=41,)

            N = int(1e5)
            x = 2*sample.radius_mm*np.random.random(N) - sample.radius_mm*np.random.random(N)
            z = 2*sample.radius_mm*np.random.random(N) - sample.radius_mm*np.random.random(N)
            mask = (x**2 + z**2 <= sample.radius_mm**2)

            x = x[mask]
            z = z[mask]

            E = loaded_field(
            solution,
            x,
            z,)
            '''E = solve_loaded_cavity(
                a,
                ell,
                n,
                Rs,
                eps_sample=epsr,
                vial_outer_radius_mm=None,
                eps_vial=1.0,
                m_max=21,
                p_max=41,
                normalization="target_coefficient",
            )["E"]'''

            E_empty = np.cos(np.pi*x/a)*np.cos(n*np.pi*z/ell)

            ## Monty
            I = np.mean(E*E_empty)*np.pi*sample.radius_mm**2
            As.append(I)

            x = 2*sample.vial_outer_radius_mm*np.random.random(N) - sample.vial_outer_radius_mm*np.random.random(N)
            z = 2*sample.vial_outer_radius_mm*np.random.random(N) - sample.vial_outer_radius_mm*np.random.random(N)
            mask = (x**2 + z**2 <= sample.vial_outer_radius_mm**2)

            x = x[mask]
            z = z[mask]

            E = loaded_field(
            solution,
            x,
            z,)

            E_empty = np.cos(np.pi*x/a)*np.cos(n*np.pi*z/ell)

            ## Monty
            Iv = np.mean(E*E_empty)*np.pi*sample.vial_outer_radius_mm**2
            Av.append(Iv)

            x = 7.03*np.random.random(N) - 7.03/2*np.random.random(N)
            z = 7.03*np.random.random(N) - 7.03/2*np.random.random(N)
            mask = (x**2 + z**2 <= (7.03/2)**2)

            x = x[mask]
            z = z[mask]

            E = loaded_field(
            solution,
            x,
            z,)

            E_empty = np.cos(np.pi*x/a)*np.cos(n*np.pi*z/ell)

            ## Monty
            IIDv = np.mean(E*E_empty)*np.pi*(7.03/2)**2

            AIDv.append(IIDv - I)





            ## The last term ostensibly accounts for the fact that the sample does not fill the entire volume
            #M2 = (epsr - 2.54)*I + (2.54 - 1)*Iv
            M2 = (epsr - sample.eps_vial)*I + (sample.eps_vial - 1)*Iv - (sample.eps_vial - 1)*(IIDv - I)

            Z.append((M1 - M2)**2)
            #Z2.append((m - M2*2/Vc))
            #print(epsr)

            #plt.ylim(-1, 1)

        print(np.round(Av[np.argmin(Z)], 6), np.round(As[np.argmin(Z)], 5), np.round(AIDv[np.argmin(Z)], 5))

        epsis = np.linspace(1e-5, 1, 10001)
        for epsi in epsis:
            Q2 = (sample.eps_vial*3.5e-4*Av[np.argmin(Z)] + (epsi - sample.eps_vial*3.5e-4)*As[np.argmin(Z)]) - sample.eps_vial*3.5e-4*(AIDv[np.argmin(Z)])
            Z2.append((Q1 - Q2)**2)

        epr.append(epsrs[np.argmin(Z)])
        #plt.plot(epsis, Z2, marker = 'o')
        #plt.show()
        #plt.plot(epsrs, Z, marker = 'o')
        print("New epsr: ", np.round(epsrs[np.argmin(Z)], 3))
        print("New loss tangent: ", epsis[np.argmin(Z2)]/epsrs[np.argmin(Z)])
        ##print(epsrs[np.argmin(Z2)])
        #plt.show()
    if (SAMPLE == "Teflon1"):
        plt.plot(f, 100*(np.array(epr) - Teflon)/Teflon, marker = 'o', linestyle = '--', label = 'Teflon1')
    elif(SAMPLE == "Rexolite"):
        plt.plot(f, 100*(np.array(epr) - Rexolite)/Rexolite, marker = 'o', linestyle = '--', label = 'Rexolite')
    elif(SAMPLE == "KelF1"):
        plt.plot(f, 100*(np.array(epr) - KelF)/KelF, marker = 'o', linestyle = '--', label = 'KelF1')
    elif(SAMPLE == "Acrylic"):
        plt.plot(f, 100*(np.array(epr) - Acrylic)/Acrylic, marker = 'o', linestyle = '--', label = 'Acrylic')
    elif(SAMPLE == "Acetal1"):
        plt.plot(f, 100*(np.array(epr) - Acetal)/Acetal, marker = 'o', linestyle = '--', label = 'Acetal1')
    elif(SAMPLE == "PVDF1"):
        plt.plot(f, 100*(np.array(epr) - PVDF)/PVDF, marker = 'o', linestyle = '--', label = 'PVDF1')
    elif(SAMPLE == "PET"):
        plt.plot(f, 100*(np.array(epr) - PET)/PET, marker = 'o', linestyle = '--', label = 'PET')






#plt.plot(f, LS1epsr, marker = 'o')
#plt.plot(f, epr, marker = 'o')
plt.grid()
plt.xlabel("Frequency (GHz)", fontsize= 15)
plt.ylabel("Percent Difference from Previous Measurements (%)", fontsize = 15)

#plt.ylim(1, 5)
plt.legend()
plt.show()




