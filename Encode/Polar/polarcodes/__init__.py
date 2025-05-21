"""
Polar Codes in Python
=============================================
"""

from polarcodes.utils import *
from polarcodes.decoder_utils import *
from polarcodes.Construct import Construct
from polarcodes.Shorten import Shorten
from polarcodes.Encode import Encode
from polarcodes.Decode import Decode
from polarcodes.AWGN import AWGN
from polarcodes.Puncture import Puncture
from polarcodes.PolarCode import PolarCode
from polarcodes.GUI import GUI

from .PolarCode import PolarCode
from .Encode import Encode
from .Decode import Decode
import streamlit as st

def polar_encode_bits(input_bits, N=128, K=64):
    """
    Encode input bits using Polar Codes in chunks.

    Args:
        input_bits (List[int]): Full binary message
        N (int): Block length (must be power of 2)
        K (int): Message length

    Returns:
        List[int]: Full encoded bit stream
    """
    st.success(f"Polar ECC is active — Encoding in chunks (K={K}, N={N})")

    encoded_bits = []
    for i in range(0, len(input_bits), K):
        chunk = input_bits[i:i+K]
        # Pad if chunk is too short
        if len(chunk) < K:
            chunk += [0] * (K - len(chunk))

        pc = PolarCode(N, K)
        pc.message = chunk
        _ = Encode(pc, encoder_name='polar_encode')
        encoded_bits.extend(pc.x.tolist())

    return encoded_bits


def polar_decode_bits(encoded_bits, N=128, K=64):
    """
    Decode encoded bits using Polar Codes in chunks.

    Args:
        encoded_bits (List[int]): Full encoded stream
        N (int): Block length
        K (int): Message length

    Returns:
        List[int]: Recovered message
    """
    st.warning(f"Polar decoding engaged — Decoding in chunks (N={N}, K={K})")

    decoded_bits = []
    for i in range(0, len(encoded_bits), N):
        chunk = encoded_bits[i:i+N]
        if len(chunk) < N:
            chunk += [0] * (N - len(chunk))  # pad if necessary

        pc = PolarCode(N, K)
        pc.x_noisy = chunk
        _ = Decode(pc, decoder_name='scd')
        decoded_bits.extend(pc.message_received.tolist())

    return decoded_bits