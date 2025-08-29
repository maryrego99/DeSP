
import zlib
from .ECCDecoder import ECCDecoder
from Encode.Helper_Functions import dna_to_int_array
from DNAFountain.CRCHandler import grand_crc_repair  

class CRCGrandDecoder(ECCDecoder):
    def __init__(self, max_flips=2):
        self.max_flips = max_flips

    def decode(self, data, original_dna=None):
        if len(data) < 4:
            return -1, None, False, False

        received_crc = int.from_bytes(data[-4:], byteorder='big')
        data_wo_crc = data[:-4]
        calc_crc = zlib.crc32(bytes(data_wo_crc))

        if calc_crc == received_crc:
            return 0, data_wo_crc, True, False  # CRC passed

        # Try GRAND repair
        if original_dna is None:
            return -1, None, False, False

        print("using CRCGrandDecoder")
        repaired_dna = grand_crc_repair(original_dna, max_flips=self.max_flips)
        if repaired_dna:
            repaired_data = dna_to_int_array(repaired_dna)
            received_crc = int.from_bytes(repaired_data[-4:], byteorder='big')
            data_wo_crc = repaired_data[:-4]
            calc_crc = zlib.crc32(bytes(data_wo_crc))
            if calc_crc == received_crc:
                return 0, data_wo_crc, False, True  # GRAND passed

        return -1, None, False, False