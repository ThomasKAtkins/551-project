import pandas as pd

dfs = []

for c in 'ABC':
    df = pd.read_csv(f'data/NMDPFrequencies/nmdpfrequencies{c}.csv')
    # remove trailing 'g' from allele names
    df['allele'] = [a[:-1] if a[-1] == 'g' else a for a in df['allele']]
    categories = ['AFA', 'API', 'NAM', 'HIS', 'CAU']
    df = df[['allele'] + [c + "_freq" for c in categories]]
    df['max_freq'] = df[[c + "_freq" for c in categories]].max(axis=1)
    # select only alleles with max_freq > 0.01
    df = df[df['max_freq'] > 0.01]
    # sort by allele
    df = df.sort_values(by='allele')
    dfs.append(df)

df = pd.concat(dfs, axis=0)
df.to_csv('data/NetMHCpan/alleles.csv', index=False)
df['allele'].to_csv('data/NetMHCpan/alleles.txt', index=False, header=False)