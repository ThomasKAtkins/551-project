import pandas as pd

file = 'data/sars-cov-2-sequences/aligned_spike_transl_filtered.fa'

with open(file, 'r') as f:
    lines = f.readlines()

sequences = []
new_seq = True

for line in lines:
    if line[0]=='>':
        new_seq = True
        continue
    if new_seq:
        sequences.append(line[:-1]) # remove newline character
        new_seq = False
    else:
        sequences[-1] += line[:-1]

peptides = []
pos = []
for sequence in sequences:
    for i in range(len(sequence)-8):
        if 'X' in sequence[i:i+9]:
            continue
        peptides.append(sequence[i:i+9])
        pos.append(i)

peptide_df = pd.DataFrame({'peptide': peptides, 'pos': pos})
peptide_df['key'] = peptide_df['peptide'] + '_' + peptide_df['pos'].astype(str)
peptide_df = peptide_df.groupby('key').count()
peptide_df['count'] = peptide_df['peptide']
peptide_df['freq'] = peptide_df['count']/len(sequences)
peptide_df['peptide'] = peptide_df.index.str.split('_').str[0]
peptide_df['pos'] = peptide_df.index.str.split('_').str[1].astype(int)
peptide_df = peptide_df.sort_values(['pos', 'freq'], ascending=[True, False])

with open("data/sars-cov-2-sequences/ancestral.fa") as f:
    ancestral = f.readlines()[0]

peptide_df['ancestral'] = False

for pos in peptide_df['pos'].unique():
    ancestral_peptide = ancestral[pos:pos+9]
    peptide_df.loc[(peptide_df['pos']==pos) & (peptide_df['peptide']==ancestral_peptide), 'ancestral'] = True

peptide_df.to_csv('data/NetMHCpan/peptides.csv', index=False)
peptide_df['peptide'].to_csv('data/NetMHCpan/input_peptides.txt', index=False, header=False)
