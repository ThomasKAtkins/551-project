# code by Masha Larina

def ltfivedels(seq):
    consecutive_count = 0
    for char in seq:
        if char == 'X':
            consecutive_count += 1
            if consecutive_count > 5:
                return False
        else:
            consecutive_count = 0
    return True

def fasta_to_arr(file_path):
    names = []
    seqs = []
    current_name = ""
    current_seq = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Save the previous sequence, if any
                if current_seq:
                    if(ltfivedels(current_seq)):
                        names.append(current_name[1:])
                        seqs.append(current_seq)
                # Start a new sequence
                current_name = line
                current_seq = ""
            else:
                # Concatenate lines of the current sequence
                current_seq += line

        # Save the last sequence
        if current_seq:
            names.append(current_name[1:])
            seqs.append(current_seq)

    return names, seqs

names, seqs = fasta_to_arr("data\\sars-cov-2-sequences\\aligned_spike_transl.fa")

ofile = open("data\\sars-cov-2-sequences\\filtered_aligned_spike_transl.fa", "w")

for i in range(len(seqs)):
    ofile.write(">" + names[i] + "\n" +seqs[i] + "\n")

ofile.close()
