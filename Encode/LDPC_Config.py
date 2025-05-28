import numpy as np
from pyldpc import make_ldpc, encode, decode

# ------------------------- CONFIG ------------------------- #
LDPC_N = 128   # Length of the codeword
LDPC_DV = 2    # Variable node degree
LDPC_DC = 4    # Check node degree
LDPC_SYSTEMATIC = True
LDPC_SEED = 42

# ------------------------- MATRIX GEN ------------------------- #
H, G = make_ldpc(LDPC_N, LDPC_DV, LDPC_DC, systematic=LDPC_SYSTEMATIC, seed=LDPC_SEED)
LDPC_K = G.shape[0]  # Message length

# ------------------------- FUNCTIONS ------------------------- #

def ldpc_encode(data_bytes):
    data_bits = np.unpackbits(np.frombuffer(data_bytes, dtype=np.uint8))
    # Pad or truncate to LDPC_K
    if len(data_bits) > LDPC_K:
        data_bits = data_bits[:LDPC_K]
    else:
        data_bits = np.pad(data_bits, (0, LDPC_K - len(data_bits)))

    codeword = encode(G, data_bits, snr=2)  # Adjust SNR if needed
    return np.packbits(codeword).tobytes()

def ldpc_decode(codeword_bytes):
    codeword_bits = np.unpackbits(np.frombuffer(codeword_bytes, dtype=np.uint8))
    decoded_bits = decode(H, codeword_bits, maxiter=50)
    return np.packbits(decoded_bits[:LDPC_K]).tobytes()