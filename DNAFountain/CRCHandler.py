from collections import Counter
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

def bit_string_length(dna_string):
    byte_data = bytes(dna_to_int_array(dna_string))
    bitstring = bytes_to_bitstring(byte_data)
    return bitstring, len(bitstring)

def pos_error_prob(length):
    # reliability high in center, low at edges
    positions = np.arange(length)
    pos_reliability = 1.0 - ((positions - length / 2) / (length / 2)) ** 2 * 0.8
    pos_reliability = np.clip(pos_reliability, 0.1, 1.0)
    bit_error_prob = 1.0 - pos_reliability
    return bit_error_prob


#Modularizing the code (extracting common logic), made decoding runtime slower
#This was affecting analysis results
#Hence separate functions for each variant

# (1)
def bitwise_heuristic_edges(dna_string, max_flips=2, strategy="brute", top_k=20, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    search_space = list(range(10)) + list(range(length-10, length)) # target edges directly without sorting

    for num_flips in range(1, max_flips + 1):
        for indices in combinations(search_space, num_flips): # sorted_indices or top_indices
            guess = bitstring
            for i in indices:
                guess = flip_bit(guess, i)
            guess_bytes = bitstring_to_bytes(guess)
            if is_valid_crc(guess_bytes):
                repaired = list(guess_bytes)
                if position_error_counter is not None:
                    position_error_counter.update(indices)
                # print("Correct guess and repaired. Converting to DNA")
                return int_array_to_dna(repaired)
            # print("Grand did not repair")

    return None


# (2)
def bitwise_heuristic_topk(dna_string, max_flips=2, top_k=20, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    bit_error_prob = pos_error_prob(length)
    sorted_indices = np.argsort(-bit_error_prob)
    top_indices = sorted_indices[:top_k]

    for num_flips in range(1, max_flips + 1):
        sorted_indices = np.argsort(-bit_error_prob)  # descending probability of error

        top_indices = sorted_indices[:20] # use sorted_indices to remove top k

        for indices in combinations(top_indices, num_flips): # sorted_indices or top_indices
            guess = bitstring
            for i in indices:
                guess = flip_bit(guess, i)
            guess_bytes = bitstring_to_bytes(guess)
            if is_valid_crc(guess_bytes):
                repaired = list(guess_bytes)
                if position_error_counter is not None:
                    position_error_counter.update(indices)
                # print("Correct guess and repaired. Converting to DNA")
                return int_array_to_dna(repaired)
            # print("Grand did not repair")

    return None

# (3)
def bitwise_heuristic_nok(dna_string, max_flips=2, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    bit_error_prob = pos_error_prob(length)
    sorted_indices = np.argsort(-bit_error_prob)

    for num_flips in range(1, max_flips + 1):
        sorted_indices = np.argsort(-bit_error_prob)  # descending probability of error

        for indices in combinations(sorted_indices, num_flips): # sorted_indices or top_indices
            guess = bitstring
            for i in indices:
                guess = flip_bit(guess, i)
            guess_bytes = bitstring_to_bytes(guess)
            if is_valid_crc(guess_bytes):
                repaired = list(guess_bytes)
                if position_error_counter is not None:
                    position_error_counter.update(indices)
                # print("Correct guess and repaired. Converting to DNA")
                return int_array_to_dna(repaired)
            # print("Grand did not repair")

# (4)
def bitwise_brute(dna_string, max_flips=2, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    for num_flips in range(1, max_flips + 1):
        for indices in combinations(range(length), num_flips):
            guess = bitstring
            for i in indices:
                guess = flip_bit(guess, i)
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

# (5)
def basewise_brute(dna_string, max_flips=1, top_k_bases=10):
    bitstring, length = bit_string_length(dna_string)

    #2 bits at a time
    base_positions = []
    for i in range(0, length, 2):
        bits = bitstring[i:i+2]
        if bits not in BIT_TO_BASE:
            print(f"Inalid bits at {i}")
            continue
        base = BIT_TO_BASE[bits]
        base_index = BASE_TO_INDEX[base]
        base_positions.append((i, base, base_index)) 

    for num_flips in range(1, max_flips + 1):
        for base_indices in combinations(base_positions, num_flips):
            def recursive_substitute(index, current_bits):
                if index == len(base_indices):
                    guess_bytes = bitstring_to_bytes(current_bits)
                    if is_valid_crc(guess_bytes):
                        print("valid crc [basewise brute]")
                        return int_array_to_dna(list(guess_bytes))
                    return None

                i, base, base_index = base_indices[index]
                for target_idx in range(4):
                    target_base = BASES[target_idx]
                    if target_base == base:
                        continue
                    alt_bits = BASE_TO_BIT[target_base]
                    new_bits = current_bits[:i] + alt_bits + current_bits[i+2:]
                    result = recursive_substitute(index + 1, new_bits)
                    if result:
                        return result
                return None

            result = recursive_substitute(0, bitstring)
            if result:
                return result

    return None


# (6)
def basewise_heuristic_topk(dna_string, TM_matrix, max_flips=1, top_k=10, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    bit_error_prob = pos_error_prob(length)

    #2 bits at a time
    base_scores = []
    for i in range(0, length, 2):
        bits = bitstring[i:i+2]
        if bits not in BIT_TO_BASE:
            print(f"Inalid bits at {i}")
            continue
        base = BIT_TO_BASE[bits]
        base_index = BASE_TO_INDEX[base]
        bit_score = (bit_error_prob[i] + bit_error_prob[i+1])/2 #position based error probability
        base_scores.append((i, base, base_index, bit_score))

    base_scores.sort(key=lambda x:-x[3])
    top_bases = base_scores[:top_k] 
    # top_bases = base_scores[:10] + base_scores[-10:]
    # print(f'Top bases: {top_bases}')

    for num_flips in range(1, max_flips + 1):
        for base_indices in combinations(top_bases, num_flips): #use base_scores to take out limit , else top_k_bases
            def recursive_substitute(index, current_bits):
                if index == len(base_indices):
                    guess_bytes = bitstring_to_bytes(current_bits)
                    if is_valid_crc(guess_bytes):
                        print("valid crc [basewise grand - topk]")
                        if position_error_counter is not None:
                            flipped_base_positions = [i // 2 for (i, _, _, _) in base_indices]
                            position_error_counter.update(flipped_base_positions)
                        return int_array_to_dna(list(guess_bytes))
                    return None

                i, base, base_index, _ = base_indices[index]
                possible_substitutions = []

                for target_index in range(4):
                    target_base = BASES[target_index]
                    if target_base == base:
                        continue
                    substitution_prob = TM_matrix[base_index][target_index]
                    possible_substitutions.append((target_base, substitution_prob))

                # Sort substitutions by likelihood (descending)
                possible_substitutions.sort(key=lambda x: -x[1])

                for alt_base, _ in possible_substitutions:
                    alt_bits = BASE_TO_BIT[alt_base]
                    new_bits = current_bits[:i] + alt_bits + current_bits[i+2:]
                    result = recursive_substitute(index + 1, new_bits)
                    if result:
                        return result

                return None

            result = recursive_substitute(0, bitstring)
            if result:
                return result

    return None


# (7)
def basewise_heuristic_nok(dna_string, TM_matrix, max_flips=1, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    bit_error_prob = pos_error_prob(length)

    #2 bits at a time
    base_scores = []
    for i in range(0, length, 2):
        bits = bitstring[i:i+2]
        if bits not in BIT_TO_BASE:
            print(f"Inalid bits at {i}")
            continue
        base = BIT_TO_BASE[bits]
        base_index = BASE_TO_INDEX[base]
        bit_score = (bit_error_prob[i] + bit_error_prob[i+1])/2 #position based error probability
        base_scores.append((i, base, base_index, bit_score))

    base_scores.sort(key=lambda x:-x[3])
    # top_bases = base_scores[:10] + base_scores[-10:]
    # print(f'Top bases: {top_bases}')

    for num_flips in range(1, max_flips + 1):
        for base_indices in combinations(base_scores, num_flips): #use base_scores to take out limit , else top_k_bases
            def recursive_substitute(index, current_bits):
                if index == len(base_indices):
                    guess_bytes = bitstring_to_bytes(current_bits)
                    if is_valid_crc(guess_bytes):
                        print("valid crc [basewise grand - nok]")
                        if position_error_counter is not None:
                            flipped_base_positions = [i // 2 for (i, _, _, _) in base_indices]
                            position_error_counter.update(flipped_base_positions)
                        return int_array_to_dna(list(guess_bytes))
                    return None

                i, base, base_index, _ = base_indices[index]
                possible_substitutions = []

                for target_index in range(4):
                    target_base = BASES[target_index]
                    if target_base == base:
                        continue
                    substitution_prob = TM_matrix[base_index][target_index]
                    possible_substitutions.append((target_base, substitution_prob))

                # Sort substitutions by likelihood (descending)
                possible_substitutions.sort(key=lambda x: -x[1])

                for alt_base, _ in possible_substitutions:
                    alt_bits = BASE_TO_BIT[alt_base]
                    new_bits = current_bits[:i] + alt_bits + current_bits[i+2:]
                    result = recursive_substitute(index + 1, new_bits)
                    if result:
                        return result

                return None

            result = recursive_substitute(0, bitstring)
            if result:
                return result

    return None


# (8)
def basewise_heuristic_edges(dna_string, TM_matrix, max_flips=2, position_error_counter=None):
    bitstring, length = bit_string_length(dna_string)

    # bit_error_prob = pos_error_prob(length)

    #2 bits at a time
    base_scores = []
    for i in range(0, length, 2):
        bits = bitstring[i:i+2]
        if bits not in BIT_TO_BASE:
            print(f"Inalid bits at {i}")
            continue
        base = BIT_TO_BASE[bits]
        base_index = BASE_TO_INDEX[base]
        # bit_score = (bit_error_prob[i] + bit_error_prob[i+1])/2 #position based error probability
        base_scores.append((i, base, base_index))

    top_bases = base_scores[:10] + base_scores[-10:]
    # print(f'Top bases: {top_bases}')

    for num_flips in range(1, max_flips + 1):
        for base_indices in combinations(top_bases, num_flips): #use base_scores to take out limit, else top_k_bases
            def recursive_substitute(index, current_bits):
                if index == len(base_indices):
                    guess_bytes = bitstring_to_bytes(current_bits)
                    if is_valid_crc(guess_bytes):
                        # print("valid crc [basewise grand - recursive]")
                        if position_error_counter is not None:
                            flipped_base_positions = [i // 2 for (i, _, _) in base_indices]
                            position_error_counter.update(flipped_base_positions)
                        return int_array_to_dna(list(guess_bytes))
                    return None

                i, base, base_index = base_indices[index]
                possible_substitutions = []

                for target_index in range(4):
                    target_base = BASES[target_index]
                    if target_base == base:
                        continue
                    substitution_prob = TM_matrix[base_index][target_index]
                    possible_substitutions.append((target_base, substitution_prob))

                # Sort substitutions by likelihood (descending)
                possible_substitutions.sort(key=lambda x: -x[1])

                for alt_base, _ in possible_substitutions:
                    alt_bits = BASE_TO_BIT[alt_base]
                    new_bits = current_bits[:i] + alt_bits + current_bits[i+2:]
                    result = recursive_substitute(index + 1, new_bits)
                    if result:
                        return result

                return None

            result = recursive_substitute(0, bitstring)
            if result:
                return result

    return None