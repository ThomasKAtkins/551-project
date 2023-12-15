# code by thomas

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dfs = []

for locus in 'ABC':
    dfs.append(pd.read_csv(f'data/NMDPfrequencies/nmdpfrequencies{locus}.csv', index_col=0))

df = pd.concat(dfs, axis=1)

# remove 'g' from allele names
df.index = df.index.str.replace('g', '')

# average all columns with _freq
mean_freq = df.filter(regex='_freq').mean(axis=1).to_frame()
mean_freq.columns = ['mean_freq']
mean_freq = mean_freq.sort_values(by='mean_freq', ascending=False)
mean_freq = mean_freq[mean_freq['mean_freq'] > 0.01]
print(mean_freq)
plt.bar(mean_freq.index, mean_freq['mean_freq'])
plt.xticks(rotation=90)
plt.ylabel('Frequency')
plt.xlabel('Allele')
plt.show()