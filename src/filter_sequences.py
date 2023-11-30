

file = 'data/sars-cov-2-sequences/aligned_spike_transl.fa'

with open(file, 'r') as f:
    lines = f.readlines()

headers = []
sequences = []
new_seq = True

for line in lines:
    if line[0]=='>':
        headers.append(line)
        new_seq = True
        continue
    if new_seq:
        sequences.append(line)
        new_seq = False
    else:
        sequences[-1] += line

good_headers = []
good_sequences = []
for header, sequence in zip(headers, sequences):
    # if more than 5 Xs in a row, skip
    if 'XXXXX' in sequence:
        continue
    good_headers.append(header)
    good_sequences.append(sequence)

out_file = 'data/sars-cov-2-sequences/aligned_spike_transl_filtered.fa'
with open(out_file, 'w') as f:
    for header, sequence in zip(good_headers, good_sequences):
        f.write(header)
        # write sequence in 80 character chunks
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80])
    
    