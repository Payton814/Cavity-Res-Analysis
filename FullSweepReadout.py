import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

Folder = '../CROAQ/data/Cav1Filled_FullSweep'
files = glob.glob(Folder + '/*.csv')

for i in range(len(files)):
    df = pd.read_csv(files[i])
    plt.plot(df.iloc[:, 0], df.iloc[:, 1], color = 'r')

plt.grid()
plt.xlabel('Frequency (GHz)')
plt.show()