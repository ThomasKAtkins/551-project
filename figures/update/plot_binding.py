# code by thomas

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('data/NetMHCpanOutput/netmhcpan_a0201.csv')

fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6, 8))
axs[0].plot(df['Pos'], df['Score_EL'], color='firebrick')
axs[1].plot(df['Pos'], df['Score_EL'], color='navy')
axs[1].set_yscale('log')
plt.xlabel('Position in SARS-CoV-2 Spike')
axs[0].set_ylabel('NetMHCpan Eluted Ligand Score')
axs[1].set_ylabel('NetMHCpan Eluted Ligand Score')
plt.tight_layout()
plt.savefig('figures/update/netmhcpan_a0201.pdf')