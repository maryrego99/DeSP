import numpy as np
from Encode.Helper_Functions import preprocess
from Encode.DNAFountain import DNAFountain, Glass
from Model.Model import DNA_Channel_Model
from Model.config import DEFAULT_PASSER
import os

# -------------------- Settings -------------------- #
file_name = "files/lena.jpg"   # Make sure this file exists
chunk_size = 20
alpha = 0.5                    # Redundancy level
crc_length = 4                 # CRC32 = 4 bytes
n_trials = 1                   # Can bump up if needed

# -------------------- Preprocess -------------------- #
data, pad = preprocess(file_name, chunk_size)

# -------------------- Setup Error Model -------------------- #
args = DEFAULT_PASSER
args.syn_sub_prob = 0.008 / 3     # 10% substitution, split across 3 bases
# args.syn_ins_prob = 0.01
# args.syn_del_prob = 0.01
# args.syn_yield = 0.98
# args.syn_number = 30
# args.decay_loss_rate = 0.3
# args.seq_copies = 1000
# args.sam_ratio = 0.01
# args.sam_to_number = 25

channel = DNA_Channel_Model(None, args)

# -------------------- Encode -------------------- #
fountain = DNAFountain(data, alpha=alpha, rs=0)  # rs=0 means no RS
fountain.encode()
fountain.save("crc/out_crc.dna")

# -------------------- Corrupt with Channel Model -------------------- #
# Read DNA strands
with open("crc/out_crc.dna", "r") as f:
    dna_strands = [line.strip() for line in f.readlines() if line.strip()]

# Pass through the channel to simulate errors
corrupted_outputs = channel(dna_strands, print_state=True)

# Extract actual DNA strings from nested structure
corrupted_strands = []
for dna in corrupted_outputs:
    for _, _, corrupted_seq in dna['re']:
        corrupted_strands.append(corrupted_seq)

# Save corrupted strands
with open("crc/out_crc_corrupted.dna", "w") as f:
    for strand in corrupted_strands:
        f.write(strand + "\n")

# -------------------- Decode -------------------- #
glass = Glass("crc/out_crc_corrupted.dna", chunk_num=len(data), rs=0)
ret, solve_num, total_reads, chunks_done, crc_fails, _, _ = glass.decode()

# If decoding was successful, save output
if glass.isDone():
    output_path = "crc/recovered_output.jpg"
    glass.save(output_path, pad=pad)
    print(f"Decoded file saved to: {output_path}")
else:
    print("Decoding incomplete. File not saved.")

# -------------------- Report -------------------- #
print("\n===== Simulation Report (CRC Only) =====")
print(f"Total Droplets Sequenced      : {total_reads}")
print(f"CRC-Failed Reads (Discarded)  : {crc_fails}")
print(f"Valid Reads Passed CRC        : {total_reads - crc_fails}")
print(f"Chunks Recovered              : {chunks_done} / {len(data)}")
print(f"Decoded Successfully?         : {'Yes' if glass.isDone() else 'No'}")