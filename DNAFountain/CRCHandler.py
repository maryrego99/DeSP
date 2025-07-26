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