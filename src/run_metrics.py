# ML

import numpy as np
from collections import Counter
import pandas as pd
from scipy.stats import entropy
import matplotlib.pyplot as plt

def fasta_to_arr(file_path):
    seqs = []
    current_seq = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seqs.append(current_seq)
                current_seq = ""
            else:
                current_seq += line

        if current_seq:
            seqs.append(current_seq)

    return seqs

def get_freq_dist(seqs):
    my_AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    counts = {i: Counter({key: 0 for key in my_AAs}) for i in range(max(map(len, seqs)))}
    
    # Construct counts
    for i in range(max(map(len, seqs))):
        for s in seqs:
            if i < len(s) and s[i] != "X":
                counts[i][s[i]] += 1
    
    vec_of_freqs = {i: [] for i in range(max(map(len, seqs)))}
    vec_of_freqs_bayes = {i: [] for i in range(max(map(len, seqs)))}
    alpha = 0.5
    for i in range(len(counts)):
        # Extract counts
        sorted_items = [t[1] for t in counts[i].items()]

        # Turn counts into frequencies
        total_count = sum(sorted_items)
        vec_of_freqs[i] = np.array([0 for _ in sorted_items]) if total_count == 0 else np.array([count / total_count for count in sorted_items])
        vec_of_freqs_bayes[i] = vec_of_freqs[i].copy()

        # Bayesian smoothing
        vec_of_freqs_bayes[i] = vec_of_freqs_bayes[i] + alpha
        vec_of_freqs_bayes[i] = vec_of_freqs_bayes[i] / np.sum(vec_of_freqs_bayes[i])
    
    return vec_of_freqs, vec_of_freqs_bayes

def KL(a, b):
    return entropy(a, b)

def Shannon(a):
    return entropy(a)

# Swiss-Prot
bg_1 = [0.0777, 0.0157, 0.053, 0.0656, 0.0405,  0.0691, 0.0227, 0.0591, 0.0595, 0.096, 0.0238, 0.0427, 0.0469, 0.0393, 0.0526, 0.0694, 0.055, 0.0667, 0.0118, 0.0311]
# Complete Human Genebank
bg_2 = [0.0704, 0.0231, 0.0484, 0.0692, 0.0378, 0.0675, 0.0256, 0.0450, 0.0565, 0.0984, 0.0237, 0.0368, 0.0610, 0.0465, 0.0552, 0.0799, 0.0534, 0.0613, 0.0121, 0.0282]

arr = fasta_to_arr("data\\sars-cov-2-sequences\\filtered_aligned_spike_transl.fa")
freq_dist, freq_dist_bayes = get_freq_dist(arr)

KL_bg1 = []
KL_bg2 = []
Shannon_entr = []
for i in range(len(freq_dist)):
    KL_bg1.append(KL(bg_1, freq_dist_bayes[i]))
    KL_bg2.append(KL(bg_2, freq_dist_bayes[i]))
    Shannon_entr.append(Shannon(freq_dist[i]))

out_df_KL = pd.DataFrame({"Swiss-Prot" : KL_bg1, "Complete Human Genebank" : KL_bg2})
out_df_Shannon = pd.DataFrame({"Shannon" : Shannon_entr})
out_df_KL.to_csv("data\\KL\\KL_divs.csv", index=False)
out_df_Shannon.to_csv("data\\entropy\\manual\\Shannon.csv", index=False)

abolish_binding = [455, 456, 459, 474, 475, 486, 490, 493, 499]
enhance_binding = [439, 452, 470, 484, 498, 501]

# Plot; limit to RBD
plt.bar(list(range(len(KL_bg2)))[306:527:1], KL_bg2[306:527:1], color="blue")
plt.bar(abolish_binding, [KL_bg1[x] for x in abolish_binding], color="red", label="abolish binding")
plt.bar(enhance_binding, [KL_bg1[x] for x in enhance_binding], color="green", label="enhance binding")
plt.legend()
plt.xlabel("Sequence position in alignment", fontsize=18)
plt.ylabel("KL Divergence", fontsize=18)
plt.title("KL divergence by position (background: Complete Human Genebank)", fontsize=20)
plt.savefig("data\\KL\\KL_SwissProt.png")
plt.show()

plt.bar(list(range(len(KL_bg2))), KL_bg2)
plt.xlabel("Sequence position in alignment")
plt.ylabel("KL Divergence")
plt.title("KL divergence by position (background: Complete Human Genebank)")
plt.savefig("data\\KL\\KL_CHG.png")
plt.show()

plt.bar(list(range(len(Shannon_entr)))[306:527:1], Shannon_entr[306:527:1], color="blue")
plt.bar(abolish_binding, [Shannon_entr[x] for x in abolish_binding], color="red")
plt.bar(enhance_binding, [Shannon_entr[x] for x in enhance_binding], color="green")
plt.scatter(abolish_binding, [Shannon_entr[x] for x in abolish_binding], color="red", label="abolish binding")
plt.scatter(enhance_binding, [Shannon_entr[x] for x in enhance_binding], color="green", label="enhance binding")
plt.legend()
plt.xlabel("Sequence position in alignment", fontsize=18)
plt.ylabel("Shannon Entropy", fontsize=18)
plt.ylim([-0.01, 0.15])
plt.title("Shannon entropy by position", fontsize=20)
plt.savefig("data\\entropy\\manual\\Shannon.png")
plt.show()
