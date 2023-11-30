import pandas as pd
import numpy as np

allele_list = pd.read_csv('data/NetMHCpan/alleles.txt', header=None)[0].tolist()
num_peptides = pd.read_csv('data/NetMHCpan/input_peptides.txt', header=None).shape[0]
print(num_peptides)

score_EL_matrix = np.zeros((num_peptides, len(allele_list)))
rank_EL_matrix = np.zeros((num_peptides, len(allele_list)))

i = 0
for allele in allele_list:
    df = pd.read_csv(f'data/NetMHCpan/output/{allele}_out.txt', delim_whitespace=True, skiprows=50, skipfooter=5, engine='python', header=None,
                        names=['Pos', 'MHC', 'Peptide', 'Core', 'Of', 'Gp', 'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score_EL', '%Rank_EL', 'BindLevel', 'bind_level'])
    score_EL_matrix[:, i] = df['Score_EL'].values
    rank_EL_matrix[:, i] = df['%Rank_EL'].values
    i += 1

np.savetxt('data/NetMHCpan/score_EL_matrix.csv', score_EL_matrix, delimiter=',')
np.savetxt('data/NetMHCpan/rank_EL_matrix.csv', rank_EL_matrix, delimiter=',')
    