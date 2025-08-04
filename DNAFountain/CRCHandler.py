import time
import zlib
import binascii
import numpy as np
from itertools import combinations
from Encode.Helper_Functions import dna_to_int_array, int_array_to_dna

def flip_bit(bitstring, index):
        flipped = list(bitstring)
        flipped[index] = '1' if flipped[index] == '0' else '0'
        return ''.join(flipped)

def bytes_to_bitstring(b):
    return ''.join(f'{byte:08b}' for byte in b)

def bitstring_to_bytes(s):
    return int(s, 2).to_bytes(len(s) // 8, byteorder='big')

def crc32_check(data_bytes):
        return binascii.crc32(data_bytes).to_bytes(4, byteorder='big')

def is_valid_crc(received_bytes):
    data, received_crc = received_bytes[:-4], received_bytes[-4:]
    return crc32_check(data) == received_crc

def grand_crc_repair(dna_string, max_flips=2):
    byte_data = bytes(dna_to_int_array(dna_string))
    bitstring = bytes_to_bitstring(byte_data)
    length = len(bitstring)

    for num_flips in range(1, max_flips + 1):
        for indices in combinations(range(length), num_flips):
            guess = bitstring
            for idx in indices:
                guess = flip_bit(guess, idx)
            guess_bytes = bitstring_to_bytes(guess)
            if is_valid_crc(guess_bytes):
                repaired = list(guess_bytes)
                # print("Correct guess and repaired. Converting to DNA")
                return int_array_to_dna(repaired)
        # print("Grand did not repair")
    return None

def heuristic_grand_crc_repair(dna_string, max_flips=2):
    byte_data = bytes(dna_to_int_array(dna_string))
    bitstring = bytes_to_bitstring(byte_data)
    length = len(bitstring)

    # reliability high in center, low at edges
    positions = np.arange(length)
    pos_reliability = 1.0 - ((positions - length / 2) / (length / 2)) ** 2 * 0.8
    pos_reliability = np.clip(pos_reliability, 0.1, 1.0)
    bit_error_probs = 1.0 - pos_reliability

    for num_flips in range(1, max_flips + 1):
        sorted_indices = np.argsort(-bit_error_probs)  # descending probability of error

        top_indices = sorted_indices[:20]
        for indices in combinations(top_indices, num_flips):
            guess = bitstring
            for idx in indices:
                guess = flip_bit(guess, idx)
            guess_bytes = bitstring_to_bytes(guess)
            if is_valid_crc(guess_bytes):
                repaired = list(guess_bytes)
                # print("Correct guess and repaired. Converting to DNA")
                return int_array_to_dna(repaired)
            # print("Grand did not repair")

    return None

BASES = ['A', 'C', 'G', 'T']
BASE_TO_INDEX = {}
for i, base in enumerate(BASES):
    BASE_TO_INDEX[base] = i

BIT_TO_BASE = {'00':'A', '01':'C', '10':'G', '11':'T'}
BASE_TO_BIT = {}
for k, v in BIT_TO_BASE.items():
    BASE_TO_BIT[v] = k

def basewise_grand_crc_repair(dna_string, TM_matrix, max_flips=1, top_k_bases=10):
    byte_data = bytes(dna_to_int_array(dna_string))
    bitstring = bytes_to_bitstring(byte_data)
    length = len(bitstring)

    positions = np.arange(length)
    pos_reliability = 1.0 - ((positions - length / 2) / (length / 2)) ** 2 * 0.8
    pos_reliability = np.clip(pos_reliability, 0.1, 1.0)
    bit_error_probs = 1.0 - pos_reliability

    #2 bits at a time
    base_scores = []
    for i in range(0, length, 2):
        bits = bitstring[i:i+2]
        if bits not in BIT_TO_BASE:
            print(f"Inalid bits at {i}")
            continue
        base = BIT_TO_BASE[bits]
        base_idx = BASE_TO_INDEX[base]
        bit_score = (bit_error_probs[i] + bit_error_probs[i+1])/2 #position based error probability
        base_scores.append((i, base, base_idx, bit_score))

    # print("Base scores:")
    # print(*base_scores)

    base_scores.sort(key=lambda x:-x[3])
    top_bases = base_scores[:top_k_bases] #take out limit 
    print(f'Top bases: {top_bases}')

    for num_flips in range(1, max_flips + 1):
        for base_indices in combinations(base_scores, num_flips):
            flipped = bitstring
            for (i, base, base_idx, bit_score) in base_indices:
                #based on TM matrix
                possible_substitutions = []
                
                for target_index in range(4):
                    target_base = BASES[target_index]
                    if target_base == base:
                        continue

                    substitution_probability = TM_matrix[base_idx][target_index]
                    possible_substitutions.append((target_base, substitution_probability))

                possible_substitutions.sort(key=lambda x: -x[1])

                for alt_base, prob in possible_substitutions:
                    alt_bits = BASE_TO_BIT[alt_base]
                    candidate = flipped[:i] + alt_bits + flipped[i+2:]
                    guess_bytes = bitstring_to_bytes(candidate)
                    if is_valid_crc(guess_bytes):
                        print("valid crc [basewise grand]")
                        return int_array_to_dna(list(guess_bytes))

    return None