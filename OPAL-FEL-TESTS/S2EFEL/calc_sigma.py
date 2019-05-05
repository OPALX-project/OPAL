import numpy as np
import pandas as pd

data = pd.read_csv('initial-profile.tsv', sep = '\t', names = ['x', 'px', 'y', 'py', 'z', 'pz'])
data = data.dropna()

sigx = np.std(data['x'])
sigy = np.std(data['y'])
sigz = np.std(data['z'])


print(sigx, sigy, sigz)

