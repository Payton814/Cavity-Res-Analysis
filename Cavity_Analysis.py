import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from getEpsr import getEpsr
from getEpsi import getEpsi
from getMur import getMur
from getMui import getMui

#df = pd.read_csv('../CROAQ/Data/Rexolite_Rod2_Cav2_Peak7_t1/trial1.csv')
df = pd.read_csv('../CROAQ/Data/Rexolite_Cav2_Peak7_t3/trial1.csv')
df1 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak1_t2/trial1.csv')
df3 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak3_t2/trial1.csv')
df5 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak5_t2/trial1.csv')
df7 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak7_t2/trial1.csv')
df2 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak2_t2/trial1.csv')
df4 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak4_t2/trial1.csv')
df6 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak6_t2/trial1.csv')
df8 = pd.read_csv('../CROAQ/Data/Rexolite_Rod3_Cav2_Peak8_t2/trial1.csv')
df = [df1, df2, df3, df4, df5, df6, df7, df8]

df1 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak1_t1/trial1.csv')
df2 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak2_t1/trial1.csv')
df3 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak3_t1/trial1.csv')
df4 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak4_t1/trial1.csv')
df5 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak5_t1/trial1.csv')
df6 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak6_t1/trial1.csv')
df7 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak7_t1/trial1.csv')
df8 = pd.read_csv('../CROAQ/Data/LSP-2_Cav2_Peak8_t1/trial1.csv')

df1 = pd.read_csv('../CROAQ/data/FillTest_Rex_Cav1_Peak7_t2/trial1.csv')
df1 = pd.read_csv('./Lunar_Data/3LR_14163,774/start_out/Peak1_start_out.csv')
df3 = pd.read_csv('./Lunar_Data/3LR_14163,774/start_out/Peak3_start_out.csv')
df5 = pd.read_csv('./Lunar_Data/3LR_14163,774/start_out/Peak5_start_out.csv')
df2 = pd.read_csv('./Lunar_Data/3LR_14163,774/start_out/Peak2_start_out.csv')
df4 = pd.read_csv('./Lunar_Data/3LR_14163,774/start_out/Peak4_start_out.csv')
df = [df1, df2, df3, df4, df5, df6, df7, df8]
df = [df1, df2, df3, df4, df5]
epsrl = []
losstanl = []
murl = []
losstanmul = []
fev = []
fodd = []
for ii in range(len(df)):
    if (ii%2==0):
        fodd.append(df[ii].iloc[0, 2])
        epsr = getEpsr(df[ii], [15,32], 13, holder = True, DEBUG=True)
        epsi = getEpsi(df[ii], [18,32], 13, holder = True, DEBUG=False)
        #print(len(epsr))
        print('epsr:', epsr)
        print('epsi:', epsi)
        print('loss tangent:', epsi/epsr)
        epsrl.append(epsr)
        losstanl.append(epsi/epsr)
    else:
        fev.append(df[ii].iloc[0, 2])
        mur = getMur(df[ii], [20,35], 13, holder=True, DEBUG=False)
        mui = getMui(df[ii], [20,35], 13, holder=True, DEBUG=False)
        print('mur:', mur)
        print('mui:', mui)
        print('loss tangent:', mui/mur)
        murl.append(mur)
        losstanmul.append(mui/mur)

fig, axs = plt.subplots(2, 2)
try:
    axs[0,0].plot(fodd, epsrl, marker = 'o', label = 'Real Permittivity')
    axs[0,0].set_ylim(0, 5)
    #axs[0,0].axhline(2.54)
    axs[0,0].grid()
    axs[0,0].legend()
    axs[0,1].plot(fodd, losstanl, marker = 'o', label = 'Electric Loss tangent')
    axs[0,1].grid()
    axs[0,1].legend()
    axs[0,1].set_yscale('log')
    axs[0,1].set_ylim(1e-5, 1)
except:
    print("No odd modes")

try:
    axs[1,0].plot(fev, murl, marker = 'o', label = 'Real Permeability')
    axs[1,0].set_ylim(0, 2)
    axs[1,0].grid()
    axs[1,0].legend()
    axs[1,1].plot(fev, losstanmul, marker = 'o', label = 'Magnetic Loss tangent')
    axs[1,1].grid()
    axs[1,1].legend()
    axs[1,1].set_yscale('log')
    axs[1,1].set_ylim(1e-6, 1e-1)
except:
    print("No even modes")

plt.suptitle('LSP-2 Sample')
plt.show()


#mur = getMur(df2, [55,65], 0.25*25.4, DEBUG=True)
#mui = getMui(df2, [55,65], 0.25*25.4, DEBUG=True)

#print('mur:', mur)
#print('mui:', mui)
#print('loss tangent:', mui/mur)

plt.errorbar([0.645265625, 0.85355625, 1.171276875, 1.51549375], [2.55, 2.57, 2.55, 2.41], yerr = [0.45, 0.18, 0.079, 0.2], marker = 'o', linestyle = '--')
plt.grid()
plt.xlabel('Frequency (GHz)')
plt.ylabel('Real Permittivity')
plt.axhline(2.54)
plt.title('Solid Rexolite Rod in Filled Cavity ($\epsilon_{b,r}$ ~ 2.635)')
plt.show()


